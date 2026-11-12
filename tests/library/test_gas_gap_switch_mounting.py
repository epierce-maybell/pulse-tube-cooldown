"""Worked inventory and heat-flow cases for mounted charcoal and connected gaps."""

from dataclasses import replace

import numpy as np
import pytest
from scipy.constants import R

from pulse_tube_cooldown.core.errors import ConfigurationError
from pulse_tube_cooldown.core.types import Gas, Material
from pulse_tube_cooldown.library.links.gas_gap_switch import (
    GasGapSwitch,
    SwitchBracket,
    load_gas_gap_switch,
)


@pytest.fixture
def single_switch():
    gas = Gas(
        name="Worked gas",
        molar_mass=R / (80 * np.pi),
        heat_capacity_ratio=3.0,
        accommodation_coeff=lambda temperature: np.full_like(temperature, 0.5),
        thermal_conductivity=lambda temperature: np.full_like(temperature, 0.02),
        saturation_pressure=lambda temperature: np.full_like(temperature, np.inf),
        adsorption_energy=10 * R * np.log(2),
        adsorption_constant=0.5,
        adsorption_capacity=3.0,
    )
    return GasGapSwitch(
        name="Worked switch",
        nodes=("Cold", "Warm"),
        gas=gas,
        fill_pressure=18.0,
        fill_temperature=30.0,
        gas_volume=5 * R,
        gap=0.01,
        area=3.0,
        charcoal_mass=1.0,
    )


@pytest.fixture
def double_switch(single_switch):
    return replace(
        single_switch,
        nodes=("Cold", "Mount", "Warm"),
        charcoal_node="Mount",
        gas=replace(single_switch.gas, adsorption_energy=12 * R * np.log(2)),
        fill_pressure=16.875,
    )


@pytest.mark.parametrize(
    ("nodes", "temperature"),
    [(("Cold", "Warm"), [5.0, 15.0]), (("Warm", "Cold"), [15.0, 5.0])],
)
def test_named_charcoal_uses_mount_temperature_after_reordering(single_switch, nodes, temperature):
    switch = replace(single_switch, nodes=nodes, charcoal_node="Cold")

    # At 5 K, affinity is 2 / Pa. At 1.5 Pa, free gas is 0.75 mol and adsorption is 2.25 mol.
    assert switch.calc_pressure(np.array(temperature)) == pytest.approx(1.5, rel=1e-12)


def test_omitting_charcoal_mount_preserves_two_node_mean_temperature(single_switch):
    # At the 10 K mean, affinity is 1 / Pa: 1 mol free + 2 mol adsorbed = 3 mol.
    assert single_switch.calc_pressure(np.array([5.0, 15.0])) == pytest.approx(2.0, rel=1e-12)


def test_connected_gaps_share_inventory_with_their_own_gas_temperatures(double_switch):
    # Fill contains 2.8125 mol. At 2 Pa, the 10 K and 16 K gaps hold 0.5 + 0.3125 mol;
    # charcoal at 12 K has affinity 1 / Pa and holds 2 mol.
    assert double_switch.calc_pressure(np.array([8.0, 12.0, 20.0])) == pytest.approx(2.0, rel=1e-12)


def test_total_volume_is_divided_between_connected_gaps(double_switch):
    switch = replace(double_switch, charcoal_mass=0.0)

    # 2.8125 mol / ((2.5 / 10 + 2.5 / 16) mol/Pa) = 90 / 13 Pa.
    assert switch.calc_pressure(np.array([8.0, 12.0, 20.0])) == pytest.approx(90 / 13, rel=1e-12)


def test_scaling_total_volume_and_charcoal_preserves_pressure(double_switch):
    switch = replace(
        double_switch,
        gas_volume=7 * double_switch.gas_volume,
        charcoal_mass=7 * double_switch.charcoal_mass,
    )

    assert switch.calc_pressure(np.array([8.0, 12.0, 20.0])) == pytest.approx(2.0, rel=1e-12)


def test_each_gap_uses_its_mean_temperature_and_balances_middle_node(double_switch):
    bracket = Material(
        name="Worked bracket",
        density=1.0,
        specific_heat=lambda temperature: np.ones_like(temperature),
        specific_heat_integral=lambda temperature: temperature,
        thermal_conductivity=lambda temperature: 0.5 * temperature,
        conductivity_integral=lambda temperature: 0.25 * temperature**2,
    )
    switch = replace(
        double_switch,
        brackets=(
            SwitchBracket(material=bracket, area=2.0, length=2.0),
            SwitchBracket(material=bracket, area=2.0, length=3.0),
        ),
    )

    terms = switch.calc_energy_balance_terms(np.array([8.0, 12.0, 20.0]))

    # First gap: molecular + viscous + bracket resistances = 1/6 + 1/6 + 1/2 K/W.
    # Second gap: 1/(3 sqrt(2.5)) + 1/6 + 5/16 K/W, from its 16 K mean.
    second_flow = 8 / (1 / (3 * np.sqrt(2.5)) + 1 / 6 + 5 / 16)
    np.testing.assert_allclose(terms.heat_flow, [4.8, second_flow - 4.8, -second_flow], rtol=1e-12)
    assert np.sum(terms.heat_flow) == pytest.approx(0.0, abs=1e-12)
    np.testing.assert_array_equal(terms.heat_capacity, 0.0)


@pytest.mark.parametrize("temperature", [[8.0, 12.0, 20.0], [12.0, 8.0, 20.0], [20.0, 12.0, 8.0]])
def test_connected_pressure_is_capped_by_coldest_surface(double_switch, temperature):
    switch = replace(
        double_switch,
        gas=replace(double_switch.gas, saturation_pressure=lambda temperature: temperature / 8),
    )

    assert switch.calc_pressure(np.array(temperature)) == pytest.approx(1.0, rel=1e-12)


@pytest.mark.parametrize("nodes", [(), ("Cold",), ("Cold", "Mount", "Warm", "Extra")])
def test_invalid_switch_node_count_raises_configuration_error(single_switch, nodes):
    with pytest.raises(ConfigurationError):
        replace(single_switch, nodes=nodes)


def test_charcoal_mount_must_belong_to_switch(single_switch):
    with pytest.raises(ConfigurationError):
        replace(single_switch, charcoal_node="Outside")


def test_three_node_switch_requires_explicit_charcoal_mount(single_switch):
    with pytest.raises(ConfigurationError):
        replace(single_switch, nodes=("Cold", "Mount", "Warm"))


@pytest.mark.parametrize("include_mount", [True, False])
def test_loader_handles_connected_switch_mounting(include_mount):
    config = {
        "name": "Connected switch",
        "nodes": ["Cold", "Mount", "Warm"],
        "gas": "He-4",
        "fill_pressure": 16.875,
        "fill_temperature": 30.0,
        "gas_volume": 5 * R,
        "gap": 0.01,
        "area": 3.0,
        "charcoal_mass": 0.0,
    }
    if include_mount:
        config["charcoal_node"] = "Mount"
        switch = load_gas_gap_switch(config)

        assert switch.nodes == ("Cold", "Mount", "Warm")
        assert switch.charcoal_node == "Mount"
        assert switch.calc_pressure(np.array([8.0, 12.0, 20.0])) == pytest.approx(90 / 13, rel=1e-12)
    else:
        with pytest.raises(ConfigurationError):
            load_gas_gap_switch(config)
