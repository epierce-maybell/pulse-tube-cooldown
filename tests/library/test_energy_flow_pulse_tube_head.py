"""Test the energy-flow pulse-tube head: loading, acoustic power, conduction limits, balances, and manifold."""

from dataclasses import replace

import numpy as np
import pytest

from pulse_tube_cooldown.config import load_link
from pulse_tube_cooldown.core.errors import ConfigurationError
from pulse_tube_cooldown.core.types import Material
from pulse_tube_cooldown.library.gases import helium
from pulse_tube_cooldown.library.links.pulse_tube.design import load_head_design
from pulse_tube_cooldown.library.links.pulse_tube.energy_flow.manifold import Manifold


@pytest.fixture
def manifold_table():
    """Return a manifold table for one F-70 compressor on one head."""
    return {
        "compressor": "F-70",
        "head_count": 1,
        "compressor_count": 1,
        "fill_pressure": 1.7e6,
        "supply_volume": 12.8e-3,
        "return_volume": 15.8e-3,
    }


@pytest.fixture
def head_table(manifold_table):
    """Return an unclamped PT1.5 link table."""
    return {
        "name": "Head",
        "type": "energy_flow_pulse_tube_head",
        "design": "PT1.5",
        "stages": ["PRP", "CFP"],
        "warm_flange_temperature": 295.0,
        "valve_frequency": 1.1,
        "manifold": manifold_table,
    }


@pytest.fixture
def clamp_table():
    """Return the clamp keys of a link table, with a contact resistance of 2 K/W."""
    return {"clamp": "RGP", "contact_resistance": 2.0}


@pytest.fixture
def clamped_head(head_table, clamp_table):
    """Return the PT1.5 head loaded from a link table with a clamp and a contact resistance of 2 K/W."""
    return load_link({**head_table, **clamp_table})


@pytest.fixture
def unclamped_head(head_table):
    """Return the PT1.5 head loaded from a link table with no clamp."""
    return load_link(head_table)


@pytest.fixture
def constant_material():
    """Return a material with constant conductivity of 2 W/(m K)."""
    return Material(
        name="Constant conductivity",
        density=1000.0,
        specific_heat=lambda temperature: np.full_like(temperature, 100.0),
        specific_heat_integral=lambda temperature: 100.0 * temperature,
        thermal_conductivity=lambda temperature: np.full_like(temperature, 2.0),
        conductivity_integral=lambda temperature: 2.0 * temperature,
    )


@pytest.fixture
def conduction_head(clamped_head, constant_material):
    """Return a clamped head with no acoustic power, no matrix gas, and walls of constant conductivity 2 W/(m K)."""
    return replace(
        clamped_head,
        wall_material=constant_material,
        effective_impedance=1.0e30,
        section_porosity=np.full_like(clamped_head.section_porosity, 1.0e-9),
        axial_conductivity_coeff=0.0,
    )


@pytest.fixture
def screen_head(clamped_head, constant_material):
    """Return a clamped head with walls of constant conductivity 2 W/(m K) and simple correlations.

    The Nusselt number is 1 + Pe^0.5, and the axial conductivity ratio is 1 + 0.5 Pe^0.5 / porosity.
    """
    return replace(
        clamped_head,
        wall_material=constant_material,
        nusselt_coeff=1.0,
        nusselt_peclet_exp=0.5,
        nusselt_porosity_exp=0.0,
        axial_conductivity_coeff=0.5,
        axial_conductivity_peclet_exp=0.5,
        axial_conductivity_porosity_exp=-1.0,
    )


def section_counts(head):
    """Return the numbers of sections in regenerator 1 and in regenerator 2 above and below the clamp."""
    return len(head.regenerator_1), len(head.regenerator_2_upper), len(head.regenerator_2_lower)


def linear_temperatures(head, warm, stage_1, stage_2):
    """Return the section temperatures of a linear profile, warm to stage 1 to stage 2, at section centers."""
    r1, _, _ = section_counts(head)
    lengths = head.section_length
    centers = np.cumsum(lengths) - lengths / 2
    regenerator_1 = warm + (stage_1 - warm) * centers[:r1] / lengths[:r1].sum()
    regenerator_2 = stage_1 + (stage_2 - stage_1) * (centers[r1:] - lengths[:r1].sum()) / lengths[r1:].sum()
    return np.concatenate([regenerator_1, regenerator_2])


def realistic_temperatures(head, stage_1=60.0, stage_2=5.0, plate=20.0):
    """Return node temperatures for a clamped or unclamped head with geometrically spaced sections."""
    sections = np.geomspace(250.0, 6.0, len(head.internal_nodes))
    clamp = [] if head.regenerator_plate is None else [plate]
    return np.array([stage_1, stage_2, *clamp, *sections])


def test_loading_clamped_head_orders_nodes_stages_plate_sections(clamped_head):
    """Verify that a clamped head lists stage 1, stage 2, the plate, then sections warm to cold."""
    sections = tuple(node.name for node in clamped_head.internal_nodes)

    assert clamped_head.nodes[:3] == ("PRP", "CFP", "RGP")
    assert clamped_head.nodes[3:] == sections
    r1, upper, lower = section_counts(clamped_head)
    design = load_head_design("PT1.5")
    assert len(sections) == r1 + upper + lower
    assert r1 >= len(design.regenerator_1)
    assert upper + lower >= len(design.regenerator_2)
    assert upper >= 1
    assert lower >= 1
    assert len(set(clamped_head.nodes)) == len(clamped_head.nodes)


def test_loading_unclamped_head_has_no_plate_and_no_contact(unclamped_head):
    """Verify that a head loaded without a clamp has no plate node and zero contact conductance."""
    assert unclamped_head.nodes[:2] == ("PRP", "CFP")
    assert len(unclamped_head.nodes) == 2 + sum(section_counts(unclamped_head))
    assert unclamped_head.clamp_contact_conductance == 0.0


def test_loading_clamp_converts_resistance_to_conductance(clamped_head):
    """Verify that a contact resistance of 2 K/W gives a contact conductance of 0.5 W/K."""
    assert clamped_head.clamp_contact_conductance == pytest.approx(0.5)


@pytest.mark.parametrize(
    "extra", [{"clamp": "RGP"}, {"contact_resistance": 2.0}, {"clamp": "RGP", "contact_resistance": 0.0}]
)
def test_loading_rejects_incomplete_or_nonpositive_clamp(head_table, extra):
    """Verify that a clamp needs both a node and a positive resistance."""
    with pytest.raises(ConfigurationError):
        load_link({**head_table, **extra})


@pytest.mark.parametrize("clamped", [True, False])
def test_balance_returns_one_heat_flow_per_node(clamped, clamped_head, unclamped_head):
    """Verify that the energy balance has an entry for each node the head lists."""
    head = clamped_head if clamped else unclamped_head

    terms = head.calc_energy_balance_terms(realistic_temperatures(head))

    assert terms.heat_flow.shape == (len(head.nodes),)
    assert np.all(np.isfinite(terms.heat_flow))


def test_chain_acoustic_power_scales_with_temperature_for_ideal_gas(unclamped_head):
    """Verify that unit transmission carries acoustic power in proportion to temperature."""
    # Points at 300, 200, 100, 50 K carry 12 W at the warm end: 12 * T / 300 = 12, 8, 4, 2 W.
    powers = unclamped_head.chain_acoustic_powers(12.0, 300.0, np.array([200.0, 100.0]), 50.0, np.array([1.0, 1.0]))

    np.testing.assert_allclose(powers, [12.0, 8.0, 4.0, 2.0], rtol=1e-12)


def test_chain_acoustic_power_steps_with_mean_transmission(unclamped_head):
    """Verify that each step uses the mean of the transmissions of its two neighboring sections."""
    # Intervals have transmission 1, 0.5, 0 from the end sections and their mean.
    powers = unclamped_head.chain_acoustic_powers(8.0, 400.0, np.array([200.0, 100.0]), 50.0, np.array([1.0, 0.0]))

    # 8 * (200/400)^1 = 4; 4 * (100/200)^0.5 = 2 sqrt(2); the last step has exponent 0.
    np.testing.assert_allclose(powers, [8.0, 4.0, 2.0 * np.sqrt(2.0), 2.0 * np.sqrt(2.0)], rtol=1e-12)


def test_acoustic_powers_split_at_stage_1(unclamped_head):
    """Verify the input power p_a^2 / (2 Z), the stage-1 split, and the powers delivered to the pulse tubes."""
    # p_a^2 / (2 Z) = (1e5)^2 / (2 * 5e8) = 10 W enters regenerator 1.
    head = replace(unclamped_head, effective_impedance=5.0e8, stage_1_split=0.25)
    sections = np.array([200.0, 120.0, 90.0, 80.0, 70.0, 30.0, 10.0])

    section_powers, pulse_tube_1, pulse_tube_2 = head.acoustic_powers(60.0, sections, 5.0, np.ones(7), 1.0e5)

    # Ideal gas: power is 10 W * T / 295 K in regenerator 1, then 0.25 of the stage-1 power * T / 60 K.
    stage_1_power = 10.0 * 60.0 / 295.0
    expected_sections = np.concatenate([10.0 * sections[:2] / 295.0, 0.25 * stage_1_power * sections[2:] / 60.0])
    np.testing.assert_allclose(section_powers, expected_sections, rtol=1e-12)
    assert pulse_tube_1 == pytest.approx(0.75 * stage_1_power, rel=1e-12)
    assert pulse_tube_2 == pytest.approx(0.25 * stage_1_power * 5.0 / 60.0, rel=1e-12)


def test_half_section_transport_wall_conduction_without_gas(conduction_head):
    """Verify that without acoustic power or matrix gas a half section's resistance is (L / 2) / (k A_wall)."""
    n = len(conduction_head.internal_nodes)
    gas = helium.properties(np.full(n, 50.0), 1.7e6)

    resistance, offset = conduction_head.half_section_transport(
        np.full(n, 50.0),
        np.zeros(n),
        1.0e5,
        gas=gas,
        matrix_fraction=np.full(n, 0.9),
        transmission=np.full(n, 0.5),
    )

    expected = 0.5 * conduction_head.section_length / (2.0 * conduction_head.section_wall_area)
    np.testing.assert_allclose(resistance, expected, rtol=1e-8)
    np.testing.assert_allclose(offset, 0.0, atol=1e-30)


def test_axial_conductivity_ratio_is_one_at_zero_peclet_number(screen_head):
    """Verify that the gas keeps its molecular conductivity with no flow and gains 0.5 Pe^0.5 / porosity with flow."""
    ratio = screen_head.axial_conductivity_ratio(np.array([0.0, 4.0]), np.array([0.5, 0.5]))

    # 1 + 0.5 * sqrt(4) / 0.5 = 3.
    np.testing.assert_allclose(ratio, [1.0, 3.0], rtol=1e-12)


def test_half_section_transport_gas_conducts_axially_without_power(screen_head):
    """Verify that with no acoustic power K = k_gas porosity A_bore + k_wall A_wall."""
    head = screen_head
    n = len(head.internal_nodes)
    gas = replace(helium.properties(np.full(n, 50.0), 1.7e6), thermal_conductivity=np.full(n, 0.5))

    resistance, _ = head.half_section_transport(
        np.full(n, 50.0),
        np.zeros(n),
        1.0e5,
        gas=gas,
        matrix_fraction=np.full(n, 0.5),
        transmission=np.full(n, 0.5),
    )

    conductance = 0.5 * head.section_porosity * head.section_flow_area + 2.0 * head.section_wall_area
    np.testing.assert_allclose(resistance, 0.5 * head.section_length / conductance, rtol=1e-12)


def test_half_section_transport_adds_gas_terms_to_wall_conductance(screen_head):
    """Verify K = alpha^2 c_p^2 mean(mdot^2) / hA' + N_k k_gas porosity A_bore + k_wall A_wall, and the offset."""
    head = screen_head
    n = len(head.internal_nodes)
    # Density 1, specific heat 2, conductivity 0.5: mdot amplitude = 1 * 2 / 4 = 0.5, mean(mdot^2) = 0.5.
    gas = replace(
        helium.properties(np.full(n, 50.0), 1.7e6),
        density=np.full(n, 1.0),
        isobaric_specific_heat=np.full(n, 2.0),
        thermal_conductivity=np.full(n, 0.5),
    )

    resistance, offset = head.half_section_transport(
        np.full(n, 50.0),
        np.full(n, 2.0),
        4.0,
        gas=gas,
        matrix_fraction=np.full(n, 0.5),
        transmission=np.full(n, 0.4),
    )

    porosity, bore, diameter = head.section_porosity[0], head.section_flow_area[0], head.section_hydraulic_diameter[0]
    void = porosity * bore
    # Pe = sqrt(0.5) / void * d * c_p / k_gas; (alpha c_p)^2 = 1.
    peclet = np.sqrt(0.5) / void * diameter * 2.0 / 0.5
    exchange = (1.0 + np.sqrt(peclet)) * 0.5 * 4.0 * void / diameter**2
    axial_gas = (1.0 + 0.5 * np.sqrt(peclet) / porosity) * 0.5 * void
    conductance = 0.5 / exchange + axial_gas + 2.0 * head.section_wall_area[0]
    assert resistance[0] == pytest.approx(0.5 * head.section_length[0] / conductance, rel=1e-12)
    assert offset[0] == pytest.approx(resistance[0] * 0.6 * 2.0, rel=1e-12)


def test_pulse_tube_energy_flow_subtracts_wall_conduction_from_acoustic_power(conduction_head):
    """Verify H = W - (A_wall / L) (I(T0) - I(T_stage)) with I(T) = 2 T."""
    head = conduction_head
    flow_1, flow_2 = head.pulse_tube_energy_flows(80.0, 4.0, 3.0, 0.5)

    conduction = head.pulse_tube_wall_area / head.pulse_tube_length * 2.0 * np.array([295.0 - 80.0, 295.0 - 4.0])
    assert flow_1 == pytest.approx(3.0 - conduction[0], rel=1e-12)
    assert flow_2 == pytest.approx(0.5 - conduction[1], rel=1e-12)


def test_pulse_tube_energy_flow_equals_acoustic_power_at_warm_temperature(conduction_head):
    """Verify that a stage at the warm-flange temperature has no wall conduction."""
    flows = conduction_head.pulse_tube_energy_flows(295.0, 295.0, 3.0, 0.5)

    assert flows == pytest.approx((3.0, 0.5), rel=1e-12)


def test_wall_conduction_limit_gives_stage_balances(conduction_head):
    """Verify that without acoustic power the stage heat flows are the wall conduction in the tubes and regenerators."""
    head = conduction_head
    r1, upper, _ = section_counts(head)
    warm, stage_1, stage_2 = 295.0, 50.0, 4.0
    sections = linear_temperatures(head, warm, stage_1, stage_2)
    junction = (
        stage_1 + (stage_2 - stage_1) * head.section_length[r1 : r1 + upper].sum() / head.section_length[r1:].sum()
    )
    temperature = np.array([stage_1, stage_2, junction, *sections])

    heat_flow = head.calc_energy_balance_terms(temperature).heat_flow

    # Conduction with k = 2 W/(m K), positive toward the cold end.
    regenerator_1 = 2.0 * head.section_wall_area[0] * (warm - stage_1) / head.section_length[:r1].sum()
    regenerator_2 = 2.0 * head.section_wall_area[-1] * (stage_1 - stage_2) / head.section_length[r1:].sum()
    tube_1, tube_2 = head.pulse_tube_wall_area / head.pulse_tube_length * 2.0 * (warm - np.array([stage_1, stage_2]))
    stage_1_cooling = -tube_1 + regenerator_2 - regenerator_1
    stage_2_cooling = -tube_2 - regenerator_2
    # Linear profile: no heat at the plate and none in any section.
    expected = np.array([-stage_1_cooling, -stage_2_cooling, 0.0, *np.zeros(len(sections))])
    np.testing.assert_allclose(heat_flow, expected, atol=1e-9)


@pytest.mark.parametrize("plate_offset", [-3.0, 0.0, 5.0])
def test_clamp_heat_flow_follows_contact_conductance(conduction_head, plate_offset):
    """Verify that the plate heat flow is G B dT / (G + B), with dT the plate offset and B the section conductance."""
    head = conduction_head
    r1, upper, _ = section_counts(head)
    warm, stage_1, stage_2 = 295.0, 50.0, 4.0
    sections = linear_temperatures(head, warm, stage_1, stage_2)
    upper_length, lower_length = head.section_length[r1 + upper - 1], head.section_length[r1 + upper]
    junction = (
        stage_1 + (stage_2 - stage_1) * head.section_length[r1 : r1 + upper].sum() / head.section_length[r1:].sum()
    )
    temperature = np.array([stage_1, stage_2, junction + plate_offset, *sections])

    heat_flow = head.calc_energy_balance_terms(temperature).heat_flow

    # Both half sections beside the clamp join it in parallel: B = 2 k A (1 / l_upper + 1 / l_lower).
    series = 2.0 * 2.0 * head.section_wall_area[-1] * (1.0 / upper_length + 1.0 / lower_length)
    contact = 0.5
    plate_cooling = contact * series * plate_offset / (contact + series)
    assert heat_flow[2] == pytest.approx(-plate_cooling, abs=1e-9)


def test_wall_conduction_heat_flows_sum_to_warm_end_balance(conduction_head):
    """Verify that all heat flows sum to the regenerator-1 warm-end flow minus the two pulse-tube flows."""
    head = conduction_head
    r1, upper, _ = section_counts(head)
    warm, stage_1, stage_2 = 295.0, 50.0, 4.0
    sections = linear_temperatures(head, warm, stage_1, stage_2)
    temperature = np.array([stage_1, stage_2, 17.0, *sections])

    heat_flow = head.calc_energy_balance_terms(temperature).heat_flow

    regenerator_1 = 2.0 * head.section_wall_area[0] * (warm - stage_1) / head.section_length[:r1].sum()
    tube_1, tube_2 = head.pulse_tube_wall_area / head.pulse_tube_length * 2.0 * (warm - np.array([stage_1, stage_2]))
    assert heat_flow.sum() == pytest.approx(regenerator_1 + tube_1 + tube_2, rel=1e-7)


def test_vanishing_contact_conductance_matches_unclamped_head(clamped_head, unclamped_head):
    """Verify that a clamp with no contact leaves the stage and section flows of the unclamped head unchanged."""
    isolated = replace(clamped_head, clamp_contact_conductance=1.0e-12)

    clamped_flows = isolated.calc_energy_balance_terms(realistic_temperatures(isolated)).heat_flow
    unclamped_flows = unclamped_head.calc_energy_balance_terms(realistic_temperatures(unclamped_head)).heat_flow

    assert clamped_flows[2] == pytest.approx(0.0, abs=1e-9)
    np.testing.assert_allclose(np.delete(clamped_flows, 2), unclamped_flows, rtol=1e-8)


def test_warmer_plate_is_cooled_by_the_clamp(clamped_head):
    """Verify that a plate warmer than the regenerator at the clamp loses heat and a colder one gains it."""
    flows = [
        clamped_head.calc_energy_balance_terms(realistic_temperatures(clamped_head, plate=plate)).heat_flow[2]
        for plate in (5.0, 40.0, 150.0)
    ]

    assert flows[0] > flows[1] > flows[2]
    assert flows[0] > 0.0 > flows[2]


def test_colder_second_stage_has_less_cooling_power(clamped_head):
    """Verify that lowering stage 2 at fixed sections reduces its cooling power."""
    cooling = [
        -clamped_head.calc_energy_balance_terms(realistic_temperatures(clamped_head, stage_2=t)).heat_flow[1]
        for t in (4.0, 8.0, 12.0)
    ]

    assert cooling[0] < cooling[1] < cooling[2]


def test_drive_is_capped_by_relief_amplitude(clamped_head):
    """Verify that two compressors on one head are limited to the relief amplitude, and one is not."""
    temperature = realistic_temperatures(clamped_head)
    sections = temperature[3:]
    one = clamped_head.drive(60.0, sections, 5.0)
    manifold = replace(clamped_head.manifold, compressor_count=2)
    two = replace(clamped_head, manifold=manifold).drive(60.0, sections, 5.0)

    assert one[1] < clamped_head.manifold.relief_amplitude
    assert two[1] == clamped_head.manifold.relief_amplitude


def test_manifold_amplitude_balances_delivery_and_charge_is_conserved():
    """Verify that the fill charge is conserved and the amplitude is the delivered flow over the admittance."""
    manifold = Manifold(
        head_count=2,
        compressor_count=1,
        displacement_rate=5.0e-3,
        volumetric_efficiency=0.9,
        relief_amplitude=1.0e9,
        fill_pressure=1.7e6,
        supply_volume=10.0e-3,
        return_volume=16.0e-3,
    )
    warm = helium.properties(295.0, 1.7e6)
    rho, slope = float(warm.density), float(warm.density * warm.isothermal_compressibility)
    head_volume, head_mass, head_mass_slope, admittance = 4.0e-3, 2.0 * rho * 4.0e-3 * 0.8, 3.0e-9, 4.0e-6

    mean_pressure, amplitude = manifold.pressures(
        warm_temperature=295.0,
        head_volume=head_volume,
        head_mass=head_mass,
        head_mass_slope=head_mass_slope,
        head_admittance=admittance,
    )

    shift = mean_pressure - 1.7e6
    charge = (
        10.0e-3 * (rho + slope * (shift + amplitude))
        + 16.0e-3 * (rho + slope * (shift - amplitude))
        + 2.0 * (head_mass + head_mass_slope * shift)
    )
    assert charge == pytest.approx(rho * (26.0e-3 + 2.0 * head_volume), rel=1e-9)
    delivery = 0.9 * 5.0e-3 * (rho + slope * (shift - amplitude))
    assert amplitude == pytest.approx(np.pi * delivery / (2.0 * admittance), rel=1e-9)


def test_manifold_symmetric_unloaded_charge_leaves_mean_pressure_at_fill():
    """Verify that equal line volumes and a head at the fill density keep the mean pressure at the fill pressure."""
    manifold = Manifold(
        head_count=1,
        compressor_count=1,
        displacement_rate=5.0e-3,
        volumetric_efficiency=0.9,
        relief_amplitude=1.0e9,
        fill_pressure=1.7e6,
        supply_volume=10.0e-3,
        return_volume=10.0e-3,
    )
    rho = float(helium.properties(295.0, 1.7e6).density)

    mean_pressure, _ = manifold.pressures(
        warm_temperature=295.0,
        head_volume=4.0e-3,
        head_mass=rho * 4.0e-3,
        head_mass_slope=1.0e-9,
        head_admittance=4.0e-6,
    )

    assert mean_pressure == pytest.approx(1.7e6, rel=1e-12)


def test_manifold_relief_cap_sets_amplitude_and_conserves_charge():
    """Verify that a delivery above the relief limit gives the relief amplitude, with the charge still conserved."""
    manifold = Manifold(
        head_count=1,
        compressor_count=2,
        displacement_rate=5.0e-3,
        volumetric_efficiency=0.9,
        relief_amplitude=1.0e5,
        fill_pressure=1.7e6,
        supply_volume=10.0e-3,
        return_volume=16.0e-3,
    )
    warm = helium.properties(295.0, 1.7e6)
    rho, slope = float(warm.density), float(warm.density * warm.isothermal_compressibility)
    head_volume, head_mass, head_mass_slope = 4.0e-3, rho * 4.0e-3 * 0.8, 3.0e-9

    mean_pressure, amplitude = manifold.pressures(
        warm_temperature=295.0,
        head_volume=head_volume,
        head_mass=head_mass,
        head_mass_slope=head_mass_slope,
        head_admittance=1.0e-7,
    )

    shift = mean_pressure - 1.7e6
    charge = (
        10.0e-3 * (rho + slope * (shift + amplitude))
        + 16.0e-3 * (rho + slope * (shift - amplitude))
        + head_mass
        + head_mass_slope * shift
    )
    assert amplitude == 1.0e5
    assert charge == pytest.approx(rho * (26.0e-3 + head_volume), rel=1e-9)


def test_manifold_more_compressors_raise_amplitude_and_more_heads_lower_it():
    """Verify that uncapped amplitude grows with compressor count and falls as heads share the supply."""

    def amplitude(heads, compressors):
        manifold = Manifold(
            head_count=heads,
            compressor_count=compressors,
            displacement_rate=5.0e-3,
            volumetric_efficiency=0.9,
            relief_amplitude=1.0e9,
            fill_pressure=1.7e6,
            supply_volume=10.0e-3,
            return_volume=10.0e-3,
        )
        rho = float(helium.properties(295.0, 1.7e6).density)
        return manifold.pressures(
            warm_temperature=295.0,
            head_volume=4.0e-3,
            head_mass=rho * 4.0e-3,
            head_mass_slope=1.0e-9,
            head_admittance=4.0e-6,
        )[1]

    assert amplitude(1, 2) > amplitude(1, 1) > amplitude(2, 1)
