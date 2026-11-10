"""Test passive coupling and shared-pressure behavior for exchange gas."""

from dataclasses import replace

import numpy as np
import pytest

from pulse_tube_cooldown.core.types import Gas
from pulse_tube_cooldown.library.links.exchange_gas import ExchangeGas, FacingSurfaces


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
        adsorption_energy=0.0,
        adsorption_constant=0.0,
        adsorption_capacity=0.0,
    )


@pytest.fixture
def link(gas):
    """Return two facing pairs at a shared pressure of 1 Pa with a separate pressure-control node."""
    return ExchangeGas(
        name="Exchange gas",
        nodes=("Can", "Pressure", "Middle", "Inner"),
        gas=gas,
        pressure_node="Pressure",
        pressure_coefficients=(0.0, 0.0, 0.0, 1.0),
        surfaces=(
            FacingSurfaces(nodes=("Can", "Middle"), area=0.004, gap=0.0005),
            FacingSurfaces(nodes=("Middle", "Inner"), area=0.006, gap=0.002),
        ),
    )


@pytest.fixture(params=[1.0e-10, 1.0, 1.0e10])
def coupled_link(request, link):
    """Return coupled surfaces at molecular, intermediate, or continuum pressure in Pa."""
    return replace(link, pressure_coefficients=(0.0, 0.0, 0.0, request.param))


@pytest.mark.parametrize(
    "temperature",
    [(300.0, 200.0, 150.0, 50.0), (50.0, 200.0, 150.0, 300.0), (150.0, 200.0, 300.0, 50.0)],
)
def test_exchange_gas_conserves_energy(coupled_link, temperature):
    """Verify that heat flows sum to zero across all nodes sharing the gas."""
    terms = coupled_link.calc_energy_balance_terms(np.array(temperature))

    assert terms.heat_flow.sum() == pytest.approx(0.0, abs=1e-12)


@pytest.mark.parametrize("temperature", [50.0, 77.0, 300.0])
def test_equal_surface_temperatures_carry_no_heat(coupled_link, temperature):
    """Verify zero flow when all facing surfaces are isothermal and the pressure node differs."""
    terms = coupled_link.calc_energy_balance_terms(np.array([temperature, 200.0, temperature, temperature]))

    np.testing.assert_array_equal(terms.heat_flow, [0.0, 0.0, 0.0, 0.0])


@pytest.mark.parametrize("pair_index, endpoint_indices", [(0, (0, 2)), (1, (2, 3))])
@pytest.mark.parametrize("temperature", [(300.0, 200.0, 150.0, 50.0), (50.0, 200.0, 150.0, 300.0)])
def test_each_facing_pair_transfers_heat_from_hot_to_cold(coupled_link, pair_index, endpoint_indices, temperature):
    """Verify hot-to-cold transfer through each isolated facing pair."""
    pair_link = replace(coupled_link, surfaces=(coupled_link.surfaces[pair_index],))
    temperature = np.array(temperature)
    first, second = endpoint_indices
    cold, hot = (first, second) if temperature[first] < temperature[second] else (second, first)

    terms = pair_link.calc_energy_balance_terms(temperature)

    assert terms.heat_flow[cold] > 0.0
    assert terms.heat_flow[hot] < 0.0


@pytest.mark.parametrize(
    "temperature",
    [(300.0, 200.0, 150.0, 50.0), (50.0, 200.0, 150.0, 300.0), (150.0, 200.0, 300.0, 50.0)],
)
def test_hottest_surface_loses_heat_and_coldest_gains_heat(coupled_link, temperature):
    """Verify passive heat transfer at the temperature extremes of the connected surfaces."""
    temperature = np.array(temperature)

    terms = coupled_link.calc_energy_balance_terms(temperature)

    surface_flows = terms.heat_flow[[0, 2, 3]]
    surface_temperatures = temperature[[0, 2, 3]]
    assert surface_flows[np.argmax(surface_temperatures)] < 0.0
    assert surface_flows[np.argmin(surface_temperatures)] > 0.0


@pytest.mark.parametrize(
    "temperature",
    [(300.0, 200.0, 150.0, 50.0), (300.0, 200.0, 50.0, 150.0), (50.0, 200.0, 300.0, 150.0)],
)
def test_combined_pairs_sum_their_separate_heat_flows(coupled_link, temperature):
    """Verify accumulation at the shared surface for opposing and reinforcing heat flows."""
    first_pair = replace(coupled_link, surfaces=(coupled_link.surfaces[0],))
    second_pair = replace(coupled_link, surfaces=(coupled_link.surfaces[1],))
    temperature = np.array(temperature)

    combined = coupled_link.calc_energy_balance_terms(temperature)
    first = first_pair.calc_energy_balance_terms(temperature)
    second = second_pair.calc_energy_balance_terms(temperature)

    np.testing.assert_allclose(combined.heat_flow, first.heat_flow + second.heat_flow, rtol=1e-12, atol=0.0)


def test_partitioning_facing_area_preserves_heat_flow(coupled_link):
    """Verify that splitting one facing area into a quarter and three quarters preserves all flows."""
    first, second = coupled_link.surfaces
    partitioned_link = replace(
        coupled_link,
        surfaces=(replace(first, area=0.001), replace(first, area=0.003), second),
    )
    temperature = np.array([300.0, 200.0, 150.0, 50.0])

    original = coupled_link.calc_energy_balance_terms(temperature)
    partitioned = partitioned_link.calc_energy_balance_terms(temperature)

    np.testing.assert_allclose(partitioned.heat_flow, original.heat_flow, rtol=1e-12, atol=0.0)


@pytest.mark.parametrize("pressure_temperature", [150.0, 200.0, 300.0])
def test_pressure_control_alone_transfers_no_heat_to_its_node(link, pressure_temperature):
    """Verify that a node controlling positive gas pressure receives heat only through facing pairs."""
    controlled_link = replace(link, pressure_coefficients=(0.0, 0.0, 0.01, -1.0))

    terms = controlled_link.calc_energy_balance_terms(np.array([300.0, pressure_temperature, 150.0, 50.0]))

    assert terms.heat_flow[1] == 0.0
    assert np.any(terms.heat_flow[[0, 2, 3]] != 0.0)


@pytest.mark.parametrize("pressure_temperature, pressure", [(100.0, 0.5), (200.0, 5.5)])
def test_pressure_node_sets_the_same_pressure_for_every_pair(link, pressure_temperature, pressure):
    """Verify that the control temperature produces the equivalent prescribed pressure across both pairs."""
    # Pressure in Pa is (temperature / 100)^3 - (temperature / 100)^2 + temperature / 100 - 0.5.
    controlled_link = replace(link, pressure_coefficients=(1.0e-6, -1.0e-4, 0.01, -0.5))
    prescribed_link = replace(link, pressure_coefficients=(0.0, 0.0, 0.0, pressure))
    temperature = np.array([300.0, pressure_temperature, 150.0, 50.0])

    controlled = controlled_link.calc_energy_balance_terms(temperature)
    prescribed = prescribed_link.calc_energy_balance_terms(temperature)

    np.testing.assert_allclose(controlled.heat_flow, prescribed.heat_flow, rtol=1e-12, atol=0.0)


@pytest.mark.parametrize("pressure_temperature", [100.0, 99.0, 50.0])
def test_zero_or_negative_pressure_relation_carries_no_heat(link, pressure_temperature):
    """Verify zero coupling when the configured pressure relation reaches or falls below zero."""
    controlled_link = replace(link, pressure_coefficients=(0.0, 0.0, 0.01, -1.0))

    terms = controlled_link.calc_energy_balance_terms(np.array([300.0, pressure_temperature, 150.0, 50.0]))

    np.testing.assert_array_equal(terms.heat_flow, [0.0, 0.0, 0.0, 0.0])


def test_all_heat_flows_tend_continuously_to_zero_at_pressure_cutoff(link):
    """Verify diminishing nonzero coupling as the pressure control approaches its zero crossing."""
    controlled_link = replace(link, pressure_coefficients=(0.0, 0.0, 0.01, -1.0))

    flows = np.array(
        [
            controlled_link.calc_energy_balance_terms(np.array([300.0, control, 150.0, 50.0])).heat_flow[[0, 2, 3]]
            for control in [101.0, 100.01, 100.0001]
        ]
    )

    assert np.all(np.abs(flows) > 0.0)
    assert np.all(np.abs(flows[1:]) < 0.02 * np.abs(flows[:-1]))


@pytest.mark.parametrize(
    "temperature, expected",
    [((50.0, 200.0, 150.0, 77.0), [8.0, 0.0, -8.0, 0.0]), ((300.0, 200.0, 100.0, 77.0), [-32.0, 0.0, 32.0, 0.0])],
)
@pytest.mark.parametrize("area, gap, scale", [(0.004, 0.0005, 1.0), (0.006, 0.002, 0.375)])
def test_single_pair_approaches_fourier_conduction(link, temperature, expected, area, gap, scale):
    """Verify the continuum Fourier limit with temperature-dependent conductivity and distinct geometries."""
    gas = replace(link.gas, thermal_conductivity=lambda temperature: temperature / 10000.0)
    pair_link = replace(
        link,
        gas=gas,
        pressure_coefficients=(0.0, 0.0, 0.0, 1.0e10),
        surfaces=(FacingSurfaces(nodes=("Can", "Middle"), area=area, gap=gap),),
    )

    terms = pair_link.calc_energy_balance_terms(np.array(temperature))

    # Fourier conduction through 0.004 m^2 over 0.0005 m gives 8 W and 32 W for the two temperature pairs.
    np.testing.assert_allclose(terms.heat_flow, scale * np.array(expected), rtol=1e-6, atol=0.0)
