"""Tests for loading gases from TOML files."""

import textwrap

import numpy as np
import pytest

from pulse_tube_cooldown.core.errors import ConfigurationError
from pulse_tube_cooldown.core.types import Gas
from pulse_tube_cooldown.library.gases.load import load_gas

# Conductivity rises linearly from 0.01 W/(m K) at 10 K to 0.1 at 100 K. Saturation pressure is
# 2e5 * exp(-4 * (5 / T - 1)) Pa up to the critical point at 5 K.
_TEST_GAS = textwrap.dedent("""
    name = "Test Gas"
    molar_mass = 4.0e-3
    heat_capacity_ratio = 1.5
    adsorption_energy = 1.0e3
    adsorption_constant = 1.0e-6
    adsorption_capacity = 10.0

    [accommodation_coeff]
    form = "table"
    table = [[10.0, 0.9], [100.0, 0.5]]
    range = [10.0, 100.0]
    source = "Accommodation reference"

    [thermal_conductivity]
    form = "table"
    table = [[10.0, 0.01], [100.0, 0.1]]
    range = [10.0, 100.0]
    source = "Conductivity reference"

    [saturation_pressure]
    form = "wagner"
    critical_temperature = 5.0
    critical_pressure = 2.0e5
    coefficients = [-4.0]
    exponents = [1.0]
    range = [2.0, inf]
    source = "Vapor pressure reference"
""")


def _write(directory, filename, text):
    (directory / filename).write_text(text, encoding="utf-8")


def _replace(old, new):
    assert old in _TEST_GAS
    return _TEST_GAS.replace(old, new)


def _saturation_pressure(temperature):
    return 2.0e5 * np.exp(-4.0 * (5.0 / temperature - 1.0))


def test_loads_gas_by_name(tmp_path):
    _write(tmp_path, "test.toml", _TEST_GAS)

    gas = load_gas("Test Gas", directory=tmp_path)

    assert isinstance(gas, Gas)
    assert gas.name == "Test Gas"
    assert gas.molar_mass == 4.0e-3
    assert gas.heat_capacity_ratio == 1.5
    np.testing.assert_allclose(gas.accommodation_coeff(np.array([10.0, 55.0, 100.0])), [0.9, 0.7, 0.5])
    assert gas.adsorption_energy == 1.0e3
    assert gas.adsorption_constant == 1.0e-6
    assert gas.adsorption_capacity == 10.0
    np.testing.assert_allclose(gas.thermal_conductivity(np.array([10.0, 55.0, 100.0])), [0.01, 0.055, 0.1])
    np.testing.assert_allclose(
        gas.saturation_pressure(np.array([2.5, 4.0, 5.0])), _saturation_pressure(np.array([2.5, 4.0, 5.0]))
    )


def test_does_not_warn_within_range_including_its_ends(tmp_path, recwarn):
    _write(tmp_path, "test.toml", _TEST_GAS)

    gas = load_gas("Test Gas", directory=tmp_path)
    gas.thermal_conductivity(np.array([10.0, 100.0]))
    gas.saturation_pressure(np.array([2.0, 5.0]))

    assert not recwarn.list


def test_thermal_conductivity_takes_value_at_nearest_range_end_with_warning(tmp_path):
    _write(tmp_path, "test.toml", _TEST_GAS)

    gas = load_gas("Test Gas", directory=tmp_path)

    with pytest.warns(RuntimeWarning, match=r"Test Gas.*thermal_conductivity.*below.*10.*100"):
        below = gas.thermal_conductivity(np.array([5.0]))
    with pytest.warns(RuntimeWarning, match=r"Test Gas.*thermal_conductivity.*above.*10.*100"):
        above = gas.thermal_conductivity(np.array([300.0]))

    np.testing.assert_allclose([below[0], above[0]], [0.01, 0.1])


def test_saturation_pressure_takes_value_at_low_range_end_with_warning(tmp_path):
    _write(tmp_path, "test.toml", _TEST_GAS)

    gas = load_gas("Test Gas", directory=tmp_path)

    with pytest.warns(RuntimeWarning, match=r"Test Gas.*saturation_pressure.*below.*2.*inf"):
        pressure = gas.saturation_pressure(np.array([1.0]))

    np.testing.assert_allclose(pressure, _saturation_pressure(np.array([2.0])))


def test_saturation_pressure_is_infinite_above_critical_point_without_warning(tmp_path, recwarn):
    _write(tmp_path, "test.toml", _TEST_GAS)

    gas = load_gas("Test Gas", directory=tmp_path)

    assert np.all(np.isposinf(gas.saturation_pressure(np.array([5.5, 300.0]))))
    assert not recwarn.list


def test_preserves_temperature_array_shape(tmp_path):
    _write(tmp_path, "test.toml", _TEST_GAS)

    gas = load_gas("Test Gas", directory=tmp_path)
    temperature = np.array([[10.0, 20.0], [3.0, 4.0]])

    assert gas.thermal_conductivity(temperature).shape == (2, 2)
    assert gas.saturation_pressure(temperature).shape == (2, 2)


def test_accepts_directory_as_string_and_integer_values(tmp_path):
    _write(tmp_path, "test.toml", _replace("adsorption_capacity = 10.0", "adsorption_capacity = 10"))

    capacity = load_gas("Test Gas", directory=str(tmp_path)).adsorption_capacity

    assert isinstance(capacity, float)
    assert capacity == 10.0


def test_unknown_name_is_reported(tmp_path):
    _write(tmp_path, "test.toml", _TEST_GAS)

    with pytest.raises(ConfigurationError, match=r"0 gases named 'Neon'"):
        load_gas("Neon", directory=tmp_path)
