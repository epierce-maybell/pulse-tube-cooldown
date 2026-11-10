"""Test steady passive behavior and gas-transport limits for gas-gap switches."""

from dataclasses import replace

import numpy as np
import pytest

from pulse_tube_cooldown.core.types import Gas
from pulse_tube_cooldown.library.links.gas_gap_switch import GasGapSwitch


@pytest.fixture
def gas():
    """Return a synthetic gas with constant conductivity of 0.02 W/(m K)."""
    return Gas(
        name="Synthetic gas",
        molar_mass=0.004,
        heat_capacity_ratio=5.0 / 3.0,
        accommodation_coeff=lambda temperature: np.full_like(temperature, 0.5),
        thermal_conductivity=lambda temperature: np.full_like(temperature, 0.02),
        saturation_pressure=lambda temperature: np.full_like(temperature, np.inf),
        adsorption_energy=1000.0,
        adsorption_constant=1.0e-8,
        adsorption_capacity=10.0,
    )


@pytest.fixture
def link(gas):
    """Return a sealed switch with facing area 0.004 m^2 and gap width 0.0005 m."""
    return GasGapSwitch(
        name="Switch",
        nodes=("Plate", "Bath"),
        gas=gas,
        fill_pressure=100000.0,
        fill_temperature=300.0,
        gas_volume=3.0e-6,
        gap=0.0005,
        area=0.004,
        charcoal_mass=0.0001,
    )


@pytest.mark.parametrize("temperature", [(4.0, 300.0), (300.0, 4.0), (50.0, 77.0), (77.0, 77.0)])
def test_steady_switch_conserves_energy(link, temperature):
    """Verify that heat entering one endpoint equals heat leaving the other at steady state."""
    terms = link.calc_energy_balance_terms(np.array(temperature))

    assert terms.heat_flow.sum() == pytest.approx(0.0, abs=1e-12)


@pytest.mark.parametrize("temperature", [4.0, 77.0, 300.0])
def test_equal_temperatures_carry_no_heat(link, temperature):
    """Verify that equal endpoint temperatures produce zero heat flow."""
    terms = link.calc_energy_balance_terms(np.array([temperature, temperature]))

    np.testing.assert_array_equal(terms.heat_flow, [0.0, 0.0])


@pytest.mark.parametrize("temperature", [(4.0, 300.0), (300.0, 4.0), (50.0, 77.0), (77.0, 50.0)])
def test_heat_flows_from_hot_to_cold(link, temperature):
    """Verify that the colder endpoint gains heat and the hotter endpoint loses heat."""
    temperature = np.array(temperature)

    terms = link.calc_energy_balance_terms(temperature)

    assert terms.heat_flow[np.argmin(temperature)] > 0.0
    assert terms.heat_flow[np.argmax(temperature)] < 0.0


@pytest.mark.parametrize("mean_temperature", [50.0, 300.0])
@pytest.mark.parametrize("direction", [-1.0, 1.0])
def test_heat_flow_tends_continuously_to_zero_as_temperatures_converge(link, mean_temperature, direction):
    """Verify diminishing nonzero heat flow as the endpoint difference shrinks toward zero."""
    differences = direction * np.array([0.1, 0.001, 0.00001])
    flows = np.array(
        [
            link.calc_energy_balance_terms(
                np.array([mean_temperature - difference / 2, mean_temperature + difference / 2])
            ).heat_flow[0]
            for difference in differences
        ]
    )

    assert np.all(direction * flows > 0.0)
    assert np.all(np.abs(flows[1:]) < 0.02 * np.abs(flows[:-1]))


@pytest.mark.parametrize("temperature", [10.0, 77.0, 300.0])
def test_vacuum_carries_no_gas_heat(link, temperature):
    """Verify zero gas-gap conductance at zero pressure."""
    assert link.calc_gap_conductance(temperature, 0.0) == 0.0


@pytest.mark.parametrize("temperature", [10.0, 77.0, 300.0])
def test_lowering_pressure_weakens_coupling_below_continuum_limit(link, temperature):
    """Verify increasing conductance with pressure, bounded by continuum conduction."""
    pressures = [0.0, 1.0e-10, 1.0e-5, 1.0, 1.0e5, 1.0e10]

    conductances = np.array([link.calc_gap_conductance(temperature, pressure) for pressure in pressures])

    assert np.all(np.diff(conductances) > 0.0)
    # Fourier conduction: 0.02 W/(m K) * 0.004 m^2 / 0.0005 m = 0.16 W/K.
    assert np.all(conductances <= 0.16 * (1.0 + 1e-12))


@pytest.mark.parametrize("temperature", [10.0, 77.0, 300.0])
@pytest.mark.parametrize("area, scale", [(0.008, 2.0), (0.002, 0.5)])
@pytest.mark.parametrize("pressure", [1.0e-10, 1.0, 1.0e10])
def test_gap_conductance_scales_with_facing_area(link, temperature, area, scale, pressure):
    """Verify that scaling facing area scales gas-gap conductance by the same factor."""
    resized_link = replace(link, area=area)

    original = link.calc_gap_conductance(temperature, pressure)
    resized = resized_link.calc_gap_conductance(temperature, pressure)

    assert resized == pytest.approx(scale * original, rel=1e-12, abs=0.0)


@pytest.mark.parametrize("temperature", [10.0, 77.0, 300.0])
@pytest.mark.parametrize("pressure, scale", [(1.0e-9, 10.0), (1.0e-11, 0.1)])
def test_molecular_conductance_scales_with_pressure(link, temperature, pressure, scale):
    """Verify pressure-proportional conductance in the low-pressure molecular limit."""
    original = link.calc_gap_conductance(temperature, 1.0e-10)
    adjusted = link.calc_gap_conductance(temperature, pressure)

    assert adjusted == pytest.approx(scale * original, rel=1e-6, abs=0.0)


@pytest.mark.parametrize("temperature", [10.0, 77.0, 300.0])
@pytest.mark.parametrize("gap", [0.00025, 0.001])
def test_molecular_conductance_is_independent_of_gap_width(link, temperature, gap):
    """Verify gap-width independence in the low-pressure molecular limit."""
    resized_link = replace(link, gap=gap)

    original = link.calc_gap_conductance(temperature, 1.0e-10)
    resized = resized_link.calc_gap_conductance(temperature, 1.0e-10)

    assert resized == pytest.approx(original, rel=1e-6, abs=0.0)


@pytest.mark.parametrize("temperature", [10.0, 77.0, 300.0])
def test_high_pressure_conductance_approaches_fourier_limit(link, temperature):
    """Verify continuum conductance of 0.16 W/K from conductivity times area over gap."""
    conductance = link.calc_gap_conductance(temperature, 1.0e10)

    assert conductance == pytest.approx(0.16, rel=1e-6, abs=0.0)


@pytest.mark.parametrize("temperature", [10.0, 77.0, 300.0])
@pytest.mark.parametrize("gap, scale", [(0.001, 0.5), (0.00025, 2.0)])
def test_continuum_conductance_scales_inversely_with_gap_width(link, temperature, gap, scale):
    """Verify inverse gap-width scaling in the high-pressure continuum limit."""
    resized_link = replace(link, gap=gap)

    original = link.calc_gap_conductance(temperature, 1.0e10)
    resized = resized_link.calc_gap_conductance(temperature, 1.0e10)

    assert resized == pytest.approx(scale * original, rel=1e-6, abs=0.0)
