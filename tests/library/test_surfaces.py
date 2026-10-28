"""Tests for loading surfaces from TOML files."""

import textwrap

import numpy as np
import pytest

from pulse_tube_cooldown.core.errors import ConfigurationError
from pulse_tube_cooldown.core.types import Surface
from pulse_tube_cooldown.library.surfaces.load import load_surface

# Emissivity rises linearly from 0.01 at 10 K to 0.03 at 100 K.
_TEST_SURFACE = textwrap.dedent("""
    name = "Test Finish"

    [emissivity]
    form = "table"
    table = [[10.0, 0.01], [100.0, 0.03]]
    range = [10.0, 100.0]
    source = "Emissivity reference"
""")


def _write(directory, filename, text):
    (directory / filename).write_text(text, encoding="utf-8")


def _replace(old, new):
    assert old in _TEST_SURFACE
    return _TEST_SURFACE.replace(old, new)


def test_loads_surface_by_name(tmp_path):
    _write(tmp_path, "test.toml", _TEST_SURFACE)

    surface = load_surface("Test Finish", directory=tmp_path)

    assert isinstance(surface, Surface)
    assert surface.name == "Test Finish"
    np.testing.assert_allclose(surface.emissivity(np.array([10.0, 55.0, 100.0])), [0.01, 0.02, 0.03])


def test_loads_analytic_emissivity(tmp_path):
    # log10(emissivity) = -3 + log10(T), so emissivity = T / 1000.
    _write(
        tmp_path,
        "test.toml",
        _replace(
            'form = "table"\ntable = [[10.0, 0.01], [100.0, 0.03]]', 'form = "log_poly"\ncoefficients = [-3.0, 1.0]'
        ),
    )

    surface = load_surface("Test Finish", directory=tmp_path)

    np.testing.assert_allclose(surface.emissivity(np.array([20.0, 50.0])), [0.02, 0.05])


def test_does_not_warn_within_range_including_its_ends(tmp_path, recwarn):
    _write(tmp_path, "test.toml", _TEST_SURFACE)

    load_surface("Test Finish", directory=tmp_path).emissivity(np.array([10.0, 100.0]))

    assert not recwarn.list


def test_emissivity_takes_value_at_nearest_range_end_with_warning(tmp_path):
    _write(tmp_path, "test.toml", _TEST_SURFACE)

    surface = load_surface("Test Finish", directory=tmp_path)

    with pytest.warns(RuntimeWarning, match=r"Test Finish.*emissivity.*below.*10.*100"):
        below = surface.emissivity(np.array([4.0]))
    with pytest.warns(RuntimeWarning, match=r"Test Finish.*emissivity.*above.*10.*100"):
        above = surface.emissivity(np.array([300.0]))

    np.testing.assert_allclose([below[0], above[0]], [0.01, 0.03])


def test_preserves_temperature_array_shape(tmp_path):
    _write(tmp_path, "test.toml", _TEST_SURFACE)

    surface = load_surface("Test Finish", directory=tmp_path)

    assert surface.emissivity(np.array([[10.0, 20.0], [30.0, 40.0]])).shape == (2, 2)


def test_finds_surface_by_name_field_from_string_directory(tmp_path):
    _write(tmp_path, "unrelated_filename.toml", _TEST_SURFACE)

    assert load_surface("Test Finish", directory=str(tmp_path)).name == "Test Finish"


def test_unknown_name_is_reported(tmp_path):
    _write(tmp_path, "test.toml", _TEST_SURFACE)

    with pytest.raises(ConfigurationError, match=r"0 surfaces named 'Black Paint'"):
        load_surface("Black Paint", directory=tmp_path)


def test_duplicate_names_are_reported(tmp_path):
    _write(tmp_path, "a.toml", _TEST_SURFACE)
    _write(tmp_path, "b.toml", _TEST_SURFACE)

    with pytest.raises(ConfigurationError, match=r"2 surfaces named 'Test Finish'"):
        load_surface("Test Finish", directory=tmp_path)


def test_malformed_emissivity_table_is_reported(tmp_path):
    _write(tmp_path, "test.toml", _replace("[[10.0, 0.01], [100.0, 0.03]]", "[[100.0, 0.01], [10.0, 0.03]]"))

    with pytest.raises(ValueError, match="strictly increasing"):
        load_surface("Test Finish", directory=tmp_path)


def test_missing_emissivity_is_reported(tmp_path):
    _write(tmp_path, "test.toml", 'name = "Test Finish"\n')

    with pytest.raises(KeyError, match="emissivity"):
        load_surface("Test Finish", directory=tmp_path)
