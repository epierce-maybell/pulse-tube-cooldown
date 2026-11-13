"""Tests for He-4 properties at pulse-tube head pressures."""

import numpy as np
import pytest
from hepak.properties import density_TP, enthalpy_TP, entropy_TP, thermal_conductivity_TP
from scipy.optimize import brentq

from pulse_tube_cooldown.library.gases.helium import (
    density,
    enthalpy,
    heat_capacity_ratio,
    isentropic_density,
    isobaric_expansion_coefficient,
    isobaric_specific_heat,
    properties,
    thermal_conductivity,
)

_GAS_CONSTANT = 8.314462618  # J/(mol K)
_MOLAR_MASS = 4.002602e-3  # kg/mol
_HEAD_PRESSURES = [0.8e6, 1.6e6, 2.4e6]  # Pa, spanning return to supply
_SECOND_VIRIAL_AT_300_K = 11.8e-6  # m^3/mol


@pytest.mark.parametrize("pressure", _HEAD_PRESSURES)
def test_near_ideal_at_room_temperature(pressure):
    temperature = 300.0
    molar_volume = _GAS_CONSTANT * temperature / pressure + _SECOND_VIRIAL_AT_300_K

    assert density(temperature, pressure) == pytest.approx(_MOLAR_MASS / molar_volume, rel=1e-3)
    assert isobaric_specific_heat(temperature, pressure) == pytest.approx(2.5 * _GAS_CONSTANT / _MOLAR_MASS, rel=1e-3)
    assert isobaric_expansion_coefficient(temperature, pressure) == pytest.approx(1.0 / temperature, rel=2e-2)


@pytest.mark.parametrize(
    ("pressure", "expected_density", "expected_specific_heat"),
    # NIST Chemistry WebBook isotherm at 4.0 K
    [(0.8e6, 148.654, 2866.90), (1.6e6, 160.184, 2461.96), (2.4e6, 168.467, 2240.26)],
)
def test_real_gas_at_four_kelvin(pressure, expected_density, expected_specific_heat):
    temperature = 4.0

    assert density(temperature, pressure) == pytest.approx(expected_density, rel=1e-2)
    assert isobaric_specific_heat(temperature, pressure) == pytest.approx(expected_specific_heat, rel=7e-2)
    assert isobaric_expansion_coefficient(temperature, pressure) * temperature < 0.5


@pytest.mark.parametrize("temperature", [2.2, 4.0, 10.0, 50.0, 300.0])
@pytest.mark.parametrize("pressure", _HEAD_PRESSURES)
def test_consistent_with_enthalpy_and_density(temperature, pressure):
    step = 1e-4 * temperature
    enthalpy_rise = enthalpy_TP(temperature + step, pressure, basis="mass") - enthalpy_TP(
        temperature - step, pressure, basis="mass"
    )
    log_density_fall = np.log(density(temperature - step, pressure)) - np.log(density(temperature + step, pressure))

    assert isobaric_specific_heat(temperature, pressure) == pytest.approx(enthalpy_rise / (2 * step), rel=1e-3)
    assert isobaric_expansion_coefficient(temperature, pressure) == pytest.approx(
        log_density_fall / (2 * step), rel=1e-3
    )


@pytest.mark.parametrize("pressure", [0.8e6, 2.4e6])
def test_isentropic_density_near_ideal_at_room_temperature(pressure):
    reference_pressure = 1.6e6
    expected = density(300.0, reference_pressure) * (pressure / reference_pressure) ** (3 / 5)

    assert isentropic_density(300.0, reference_pressure, pressure) == pytest.approx(expected, rel=5e-3)
    assert heat_capacity_ratio(300.0, pressure) == pytest.approx(5 / 3, rel=2e-3)


@pytest.mark.parametrize("temperature", [4.0, 10.0, 50.0])
@pytest.mark.parametrize("pressure", [0.8e6, 2.4e6])
def test_isentropic_density_keeps_entropy(temperature, pressure):
    reference_pressure = 1.6e6
    entropy = entropy_TP(temperature, reference_pressure, basis="mass")
    final_temperature = brentq(
        lambda t: entropy_TP(t, pressure, basis="mass") - entropy, 0.5 * temperature, 2 * temperature, xtol=1e-12
    )

    assert isentropic_density(temperature, reference_pressure, pressure) == pytest.approx(
        density_TP(final_temperature, pressure, basis="mass"), rel=1e-6
    )


def test_isentropic_density_at_reference_pressure_is_density():
    assert isentropic_density(20.0, 1.6e6, 1.6e6) == pytest.approx(density(20.0, 1.6e6), rel=1e-9)


@pytest.mark.parametrize("pressure", _HEAD_PRESSURES)
def test_enthalpy_rise_near_ideal_at_room_temperature(pressure):
    # h(300 K) - h(250 K) = 5/2 R / M * 50 K for a monatomic ideal gas.
    rise = enthalpy(300.0, pressure) - enthalpy(250.0, pressure)

    assert rise == pytest.approx(2.5 * _GAS_CONSTANT / _MOLAR_MASS * 50.0, rel=1e-3)


@pytest.mark.parametrize("temperature", [4.0, 10.0, 50.0])
@pytest.mark.parametrize("pressure", _HEAD_PRESSURES)
def test_enthalpy_rises_with_isobaric_specific_heat(temperature, pressure):
    step = 1e-3 * temperature
    rise = enthalpy(temperature + step, pressure) - enthalpy(temperature - step, pressure)

    assert rise / (2 * step) == pytest.approx(isobaric_specific_heat(temperature, pressure), rel=1e-3)


@pytest.mark.parametrize(
    "function",
    [
        density,
        enthalpy,
        isobaric_specific_heat,
        isobaric_expansion_coefficient,
        heat_capacity_ratio,
        lambda temperature, pressure: isentropic_density(temperature, 1.6e6, pressure),
    ],
)
def test_broadcasts_over_arrays(function):
    temperature = np.array([[2.0], [40.0], [300.0]])
    pressure = np.array([0.8e6, 2.4e6])

    values = function(temperature, pressure)

    assert values.shape == (3, 2)
    assert values[1, 0] == pytest.approx(function(40.0, 0.8e6))
    assert np.all(np.isfinite(values))


@pytest.mark.parametrize("pressure", [0.8e6, 1.7e6, 2.5e6])
@pytest.mark.parametrize("temperature", [1.5, 2.5, 3.0])
def test_thermal_conductivity_below_3p5_kelvin_holds_its_3p5_kelvin_value(temperature, pressure):
    """Verify that conductivity below 3.5 K is positive and equals the value at 3.5 K at the same pressure."""
    held = thermal_conductivity(3.5, pressure)

    assert held > 0.0
    assert thermal_conductivity(temperature, pressure) == pytest.approx(held, rel=1e-12)
    assert properties(temperature, pressure).thermal_conductivity == pytest.approx(held, rel=1e-12)


@pytest.mark.parametrize("pressure", [0.8e6, 2.5e6])
def test_thermal_conductivity_above_3p5_kelvin_follows_hepak(pressure):
    """Verify that conductivity above 3.5 K is HEPAK's at 10 K and 20 K."""
    assert thermal_conductivity(10.0, pressure) == pytest.approx(thermal_conductivity_TP(10.0, pressure), rel=1e-12)
    assert properties(20.0, pressure).thermal_conductivity == pytest.approx(
        thermal_conductivity_TP(20.0, pressure), rel=1e-12
    )
