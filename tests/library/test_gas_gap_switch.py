"""Gas-gap switch requirements from DESIGN's adsorption and series-resistance equations."""

from dataclasses import replace

import numpy as np
import pytest
from scipy.constants import R

from pulse_tube_cooldown.core.types import Gas, Material
from pulse_tube_cooldown.library.links.gas_gap_switch import GasGapSwitch, SwitchBracket, calc_gap_conductance


def _gas():
    # At 10 K: K = 1 / Pa, and the molecular conductance per area equals pressure numerically.
    return Gas(
        name="Worked example",
        molar_mass=R / (80 * np.pi),
        heat_capacity_ratio=3.0,
        accommodation_coeff=lambda temperature: np.full_like(temperature, 0.5),
        thermal_conductivity=lambda temperature: np.full_like(temperature, 0.02),
        saturation_pressure=lambda temperature: np.full_like(temperature, np.inf),
        adsorption_energy=10 * R * np.log(2),
        adsorption_constant=0.5,
        adsorption_capacity=3.0,
    )


def _switch():
    bracket = Material(
        name="Worked bracket",
        density=1.0,
        specific_heat=lambda temperature: np.ones_like(temperature),
        thermal_conductivity=lambda temperature: 0.5 * temperature,
        conductivity_integral=lambda temperature: 0.25 * temperature**2,
    )
    return GasGapSwitch(
        name="Switch",
        nodes=("Cold", "Warm"),
        gas=_gas(),
        fill_pressure=18.0,
        fill_temperature=30.0,
        gas_volume=5 * R,
        gap=0.01,
        area=3.0,
        charcoal_mass=1.0,
        brackets=(
            SwitchBracket(material=bracket, area=2.0, length=2.0),
            SwitchBracket(material=bracket, area=2.0, length=3.0),
        ),
    )


def test_pressure_matches_conserved_free_and_adsorbed_gas():
    # n = 3 mol; at 10 K and 2 Pa, free gas = 1 mol and adsorbed gas = 3 * 2 / 3 = 2 mol.
    switch = _switch()

    assert switch.calc_pressure(np.array([8.0, 12.0])) == pytest.approx(2.0)


def test_without_charcoal_pressure_follows_ideal_gas_cooling():
    switch = replace(_switch(), charcoal_mass=0.0)

    assert switch.calc_pressure(np.array([8.0, 12.0])) == pytest.approx(6.0)


def test_pressure_is_capped_at_colder_nodes_vapor_pressure():
    switch = _switch()
    switch = replace(switch, gas=replace(switch.gas, saturation_pressure=lambda temperature: temperature / 8))

    assert switch.calc_pressure(np.array([8.0, 12.0])) == pytest.approx(1.0)
    assert switch.calc_pressure(np.array([12.0, 8.0])) == pytest.approx(1.0)


@pytest.mark.parametrize(("pressure", "expected"), [(0.0, 0.0), (1e-12, 3e-12), (1e12, 6.0)])
def test_gap_conductance_reaches_vacuum_molecular_and_viscous_limits(pressure, expected):
    # Molecular limit: area * pressure. Viscous limit: area * k / gap = 6 W/K.
    conductance = calc_gap_conductance(_gas(), 10.0, pressure, area=3.0, gap=0.01)

    assert conductance == pytest.approx(expected, rel=1e-10, abs=0.0)


@pytest.mark.parametrize(("temperature", "expected"), [([8.0, 12.0], 4.8), ([12.0, 8.0], -4.8)])
def test_heat_flow_combines_gap_and_bracket_resistances_at_mean_temperature(temperature, expected):
    # At 10 K and 2 Pa, both gap coefficients are 2 W/(m^2 K); gap resistance = 1/3 K/W.
    # Bracket k = 5 W/(m K), resistance = 1/2 K/W. Total flow = 4 / (1/3 + 1/2) = 4.8 W.
    flow = _switch().calc_heat_flow(np.array(temperature))

    np.testing.assert_allclose(flow, [expected, -expected], rtol=1e-12)


def test_isothermal_switch_has_zero_heat_flow():
    np.testing.assert_array_equal(_switch().calc_heat_flow(np.array([10.0, 10.0])), [0.0, 0.0])


@pytest.mark.parametrize(("charcoal_mass", "expected"), [(2.0, 0.0), (0.5, 0.0015)])
def test_strong_adsorption_stays_finite_with_sufficient_or_exhausted_charcoal(charcoal_mass, expected):
    # At 0.005 K, inverse affinity underflows to zero. With 3 mol total and 1000 mol/Pa of free gas,
    # 6 mol of adsorption capacity traps everything; 1.5 mol capacity leaves 1.5 / 1000 Pa free.
    switch = replace(_switch(), charcoal_mass=charcoal_mass)

    with np.errstate(over="raise", invalid="raise", divide="raise"):
        pressure = switch.calc_pressure(np.array([0.004, 0.006]))
        flow = switch.calc_heat_flow(np.array([0.004, 0.006]))

    assert pressure == pytest.approx(expected, rel=1e-12, abs=0.0)
    assert np.all(np.isfinite(flow))
    if expected == 0.0:
        np.testing.assert_array_equal(flow, [0.0, 0.0])
