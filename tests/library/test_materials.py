"""Tests for loading materials from TOML files."""

import textwrap
import warnings

import numpy as np
import pytest

from pulse_tube_cooldown.core.errors import ConfigurationError
from pulse_tube_cooldown.core.types import Material
from pulse_tube_cooldown.library.materials.load import load_material

# Specific heat is log10(y) = log10(T), so y = T. Thermal conductivity is log10(y) = 1 + sqrt(T), so 1e3 at T = 4 and
# 1e4 at T = 9.
_POLY_MATERIAL = """
    name = "Poly Metal"

    [specific_heat]
    form = "log_poly"
    coefficients = [0.0, 1.0]
    range = [4.0, 300.0]
    source = "Specific heat reference"

    [thermal_conductivity]
    form = "rational_poly"
    coefficients = [1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    range = [4.0, 9.0]
    source = "Conductivity reference"
"""

_TABLE_MATERIAL = """
    name = "Table Metal"

    [specific_heat]
    form = "table"
    table = [[10.0, 1.0], [20.0, 3.0], [40.0, 7.0]]
    range = [10.0, 40.0]
    source = "Specific heat table"

    [thermal_conductivity]
    form = "table"
    table = [[10.0, 100.0], [40.0, 400.0]]
    range = [10.0, 40.0]
    source = "Conductivity table"
"""


def _write(directory, filename, text):
    (directory / filename).write_text(textwrap.dedent(text), encoding="utf-8")


def test_loads_analytic_fits_by_name(tmp_path):
    _write(tmp_path, "poly.toml", _POLY_MATERIAL)

    material = load_material("Poly Metal", directory=tmp_path)

    assert isinstance(material, Material)
    assert material.name == "Poly Metal"
    np.testing.assert_allclose(material.specific_heat(np.array([10.0, 100.0])), [10.0, 100.0], rtol=1e-12)
    np.testing.assert_allclose(material.thermal_conductivity(np.array([4.0, 9.0])), [1.0e3, 1.0e4], rtol=1e-12)


def test_loads_table_fits_by_name(tmp_path):
    _write(tmp_path, "table.toml", _TABLE_MATERIAL)

    material = load_material("Table Metal", directory=tmp_path)

    np.testing.assert_allclose(material.specific_heat(np.array([15.0, 30.0])), [2.0, 5.0])
    np.testing.assert_allclose(material.thermal_conductivity(np.array([25.0])), [250.0])


def test_evaluates_at_nearest_range_end_outside_range(tmp_path):
    _write(tmp_path, "poly.toml", _POLY_MATERIAL)

    material = load_material("Poly Metal", directory=tmp_path)

    # Specific heat equals T inside [4, 300] K; conductivity is 1e3 at 4 K and 1e4 at 9 K.
    with pytest.warns(RuntimeWarning):
        specific_heat = material.specific_heat(np.array([1.0, 4.0, 300.0, 1000.0]))
        thermal_conductivity = material.thermal_conductivity(np.array([2.0, 16.0]))

    np.testing.assert_allclose(specific_heat, [4.0, 4.0, 300.0, 300.0])
    np.testing.assert_allclose(thermal_conductivity, [1.0e3, 1.0e4], rtol=1e-12)


@pytest.mark.parametrize(
    ("temperature", "side"),
    [([1.0, 10.0], "below"), ([10.0, 1000.0], "above")],
)
def test_warns_naming_material_property_and_side_outside_range(tmp_path, temperature, side):
    _write(tmp_path, "poly.toml", _POLY_MATERIAL)

    material = load_material("Poly Metal", directory=tmp_path)

    with pytest.warns(RuntimeWarning, match=rf"Poly Metal.*specific_heat.*{side}.*4.*300"):
        material.specific_heat(np.array(temperature))


def test_does_not_warn_within_range_including_its_ends(tmp_path, recwarn):
    _write(tmp_path, "poly.toml", _POLY_MATERIAL)

    material = load_material("Poly Metal", directory=tmp_path)
    material.specific_heat(np.array([4.0, 100.0, 300.0]))
    material.thermal_conductivity(np.array([4.0, 9.0]))

    assert not recwarn.list


def test_repeated_evaluation_outside_range_warns_once_per_side_under_default_filter(tmp_path):
    _write(tmp_path, "poly.toml", _POLY_MATERIAL)
    material = load_material("Poly Metal", directory=tmp_path)

    # A solver evaluates the same property from one line at many temperatures.
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("default")
        for temperature in [1.0, 2.0, 3.0, 400.0, 500.0]:
            material.specific_heat(np.array([temperature]))

    assert len(caught) == 2


def test_accepts_directory_as_string(tmp_path):
    _write(tmp_path, "poly.toml", _POLY_MATERIAL)

    assert load_material("Poly Metal", directory=str(tmp_path)).name == "Poly Metal"


def test_accepts_integer_range_bounds(tmp_path):
    _write(tmp_path, "poly.toml", _POLY_MATERIAL.replace("range = [4.0, 300.0]", "range = [4, 300]"))

    material = load_material("Poly Metal", directory=tmp_path)

    with pytest.warns(RuntimeWarning):
        np.testing.assert_allclose(material.specific_heat(np.array([1.0, 1000.0])), [4.0, 300.0])


def test_preserves_temperature_array_shape(tmp_path):
    _write(tmp_path, "poly.toml", _POLY_MATERIAL)

    material = load_material("Poly Metal", directory=tmp_path)

    np.testing.assert_allclose(
        material.specific_heat(np.array([[4.0, 10.0], [100.0, 300.0]])), [[4.0, 10.0], [100.0, 300.0]]
    )


def test_finds_material_by_name_field_rather_than_filename(tmp_path):
    _write(tmp_path, "a.toml", _TABLE_MATERIAL)
    _write(tmp_path, "b.toml", _POLY_MATERIAL)

    assert load_material("Poly Metal", directory=tmp_path).name == "Poly Metal"
    assert load_material("Table Metal", directory=tmp_path).name == "Table Metal"


def test_unknown_name_lists_available_materials(tmp_path):
    _write(tmp_path, "poly.toml", _POLY_MATERIAL)
    _write(tmp_path, "table.toml", _TABLE_MATERIAL)

    with pytest.raises(ConfigurationError, match=r"Copper.*Poly Metal.*Table Metal"):
        load_material("Copper", directory=tmp_path)


def test_duplicate_names_identify_both_files(tmp_path):
    _write(tmp_path, "first.toml", _POLY_MATERIAL)
    _write(tmp_path, "second.toml", _POLY_MATERIAL)

    with pytest.raises(ConfigurationError, match=r"Poly Metal.*first\.toml.*second\.toml"):
        load_material("Poly Metal", directory=tmp_path)


@pytest.mark.parametrize(
    ("text", "message"),
    [
        ('name = "Broken\n', r"broken\.toml"),
        ('source = "No name"\n', r"broken\.toml.*name"),
        ("name = 3\n", r"broken\.toml.*name"),
    ],
)
def test_unreadable_file_in_directory_is_identified(tmp_path, text, message):
    _write(tmp_path, "poly.toml", _POLY_MATERIAL)
    _write(tmp_path, "broken.toml", text)

    with pytest.raises(ConfigurationError, match=message):
        load_material("Poly Metal", directory=tmp_path)


def _replace(old, new):
    assert old in _POLY_MATERIAL
    return _POLY_MATERIAL.replace(old, new)


@pytest.mark.parametrize(
    ("text", "message"),
    [
        (_POLY_MATERIAL.split("[thermal_conductivity]")[0], r"thermal_conductivity"),
        (_replace('name = "Poly Metal"', 'name = "Poly Metal"\ndensity = 8960.0'), r"density"),
        (_replace('form = "log_poly"\n', ""), r"specific_heat.*form"),
        (_replace('form = "log_poly"', 'form = "spline"'), r"specific_heat.*spline.*log_poly.*rational_poly.*table"),
        (_replace("range = [4.0, 300.0]\n", ""), r"specific_heat.*range"),
        (_replace("range = [4.0, 300.0]", "range = [4.0]"), r"specific_heat.*range"),
        (_replace("range = [4.0, 300.0]", "range = [300.0, 4.0]"), r"specific_heat.*range"),
        (_replace("range = [4.0, 300.0]", "range = [0.0, 300.0]"), r"specific_heat.*range"),
        (_replace("range = [4.0, 300.0]", 'range = ["low", "high"]'), r"specific_heat.*range"),
        (_replace('source = "Specific heat reference"\n', ""), r"specific_heat.*source"),
        (_replace('source = "Specific heat reference"', 'source = ""'), r"specific_heat.*source"),
        (
            _replace("coefficients = [1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]", "coefficients = [1.0]"),
            r"thermal_conductivity.*9 coefficients",
        ),
        (_replace("coefficients = [0.0, 1.0]", "table = [[4.0, 1.0], [300.0, 2.0]]"), r"specific_heat"),
    ],
    ids=[
        "missing property",
        "unknown key",
        "missing form",
        "unknown form",
        "missing range",
        "short range",
        "reversed range",
        "nonpositive range",
        "nonnumeric range",
        "missing source",
        "empty source",
        "invalid fit parameters",
        "unexpected fit parameter",
    ],
)
def test_malformed_material_reports_file_and_problem(tmp_path, text, message):
    _write(tmp_path, "poly.toml", text)

    with pytest.raises(ConfigurationError, match=rf"poly\.toml.*{message}"):
        load_material("Poly Metal", directory=tmp_path)
