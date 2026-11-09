"""Requirements for one head's gas inventory, mass swing, acoustic power, breakthrough loss, and heat flows."""

import dataclasses
import math

import numpy as np
import pytest
from scipy.integrate import quad

from pulse_tube_cooldown.core.types import Material
from pulse_tube_cooldown.library.gases import helium
from pulse_tube_cooldown.library.heads import HeadModel, MatrixLayer, PulseTube, RegeneratorSection, load_head_model
from pulse_tube_cooldown.library.links.cooler.head import (
    Head,
    calc_acoustic_power,
    calc_breakthrough_loss,
    calc_gas_mass,
    calc_heat_flow,
    calc_mass_swing,
    calc_pore_heat_capacity,
)

_MOLAR_MASS = 4.002602e-3  # kg/mol
_GAS_CONSTANT = 8.314462618  # J/(mol K)
_SECOND_VIRIAL_AT_300_K = 11.8e-6  # m^3/mol
_MATRIX = Material(
    name="Matrix",
    density=1.0,
    specific_heat=lambda t: t,
    specific_heat_integral=lambda t: t**2 / 2.0,
    thermal_conductivity=None,
    conductivity_integral=None,
)
_WALL_CONDUCTIVITY = 2.0  # W/(m K)
_WALL = Material(
    name="Wall",
    density=1.0,
    specific_heat=lambda t: t,
    specific_heat_integral=lambda t: t**2 / 2.0,
    thermal_conductivity=lambda t: np.full_like(t, _WALL_CONDUCTIVITY),
    conductivity_integral=lambda t: _WALL_CONDUCTIVITY * np.asarray(t),
)


def _model(reservoir_volume=5.0e-3, heat_exchanger_voids=(0.0, 0.0), matrix=_MATRIX):
    """A head whose stage-2 section's second layer splits into two nodes, giving four nodes in all."""

    def section(*layers):
        return RegeneratorSection(
            length=0.2,
            wall_material=_WALL,
            wall_area=1e-4,
            layers=tuple(
                MatrixLayer(material=matrix, mass=1.0, pore_volume=volume, length_fraction=fraction, node_count=count)
                for volume, fraction, count in layers
            ),
        )

    def tube(volume):
        return PulseTube(
            length=0.3, volume=volume, wall_material=_WALL, wall_area=1e-5, reservoir_volume=reservoir_volume
        )

    return HeadModel(
        name="Test",
        regenerators=(section((1.0e-4, 1.0, 1)), section((3.0e-5, 0.5, 1), (4.0e-5, 0.5, 2))),
        pulse_tubes=(tube(6.0e-5), tube(2.0e-4)),
        clamp_fraction=0.5,
        dead_volume=4.5e-5,
        waveform_factor=math.sqrt(2),
        buffer_fraction=0.5,
        stroke_spread=0.1,
        heat_exchanger_voids=heat_exchanger_voids,
    )


# Chain order: warm flange, stage-1 section node, stage 1, stage-2 section's three nodes, stage 2.
_COLD = np.array([300.0, 60.0, 40.0, 20.0, 10.0, 6.0, 4.0])
_PORE = {1: 1.0e-4, 3: 3.0e-5, 4: 2.0e-5, 5: 2.0e-5}  # m^3 at each node's chain index
_TUBE = {2: 6.0e-5, 6: 2.0e-4}  # m^3 at each stage's chain index
_NO_VOIDS = (0.0, 0.0)
_VOIDS = (4.0e-5, 1.0e-5)  # m^3 of the stage-1 and stage-2 heat exchangers
_FREQUENCY = 1.2
_STIFF = (1e10, 3e10)  # Pa s/m^3, holding each reservoir near the mean pressure
_OPEN = (2e7, 6e7)  # Pa s/m^3, letting each reservoir swing


def _mean_over_tube(function, warm, cold):
    """Mean of `function(T)` over a linear temperature profile from `warm` to `cold`, by adaptive quadrature."""
    return quad(function, cold, warm, epsrel=1e-10, limit=200)[0] / (warm - cold)


@pytest.mark.parametrize("voids", [_NO_VOIDS, _VOIDS])
def test_gas_mass_at_uniform_temperature_fills_every_volume(voids):
    # HEPAK's density at 300 K and 1.7 MPa times the pores, tubes, reservoirs, dead volume, and heat-exchanger voids.
    volume = 1.7e-4 + 2.6e-4 + 1.0e-2 + 4.5e-5 + sum(voids)

    assert calc_gas_mass(_model(heat_exchanger_voids=voids), np.full(7, 300.0), 1.7e6) == pytest.approx(
        2.7069548683 * volume, rel=1e-8
    )


@pytest.mark.parametrize("voids", [_NO_VOIDS, _VOIDS])
def test_gas_mass_places_each_volume_at_its_temperature(voids):
    pressure = 1.4e6
    pores = sum(volume * float(helium.density(_COLD[index], pressure)) for index, volume in _PORE.items())
    tubes = sum(
        volume * _mean_over_tube(lambda t: float(helium.density(t, pressure)), _COLD[0], _COLD[index])
        for index, volume in _TUBE.items()
    )
    heat_exchangers = sum(
        volume * float(helium.density(_COLD[index], pressure)) for index, volume in zip(_TUBE, voids, strict=True)
    )
    warm = (1.0e-2 + 4.5e-5) * float(helium.density(_COLD[0], pressure))

    assert calc_gas_mass(_model(heat_exchanger_voids=voids), _COLD, pressure) == pytest.approx(
        pores + tubes + heat_exchangers + warm, rel=1e-5
    )


def test_pore_heat_capacity_holds_pores_and_voids_at_their_temperatures():
    # Pore gas at each node and each heat-exchanger void at its stage, V rho c_p, with none at the warm flange.
    pressure = 1.4e6
    volume = np.zeros(7)
    volume[list(_PORE)] = list(_PORE.values())
    volume[list(_TUBE)] = _VOIDS
    expected = volume * helium.density(_COLD, pressure) * helium.isobaric_specific_heat(_COLD, pressure)

    heat_capacity = calc_pore_heat_capacity(_model(heat_exchanger_voids=_VOIDS), _COLD, pressure)

    np.testing.assert_allclose(heat_capacity, expected, rtol=1e-12)
    assert heat_capacity[0] == 0.0


@pytest.mark.parametrize("voids", [_NO_VOIDS, _VOIDS])
def test_mass_swing_in_ideal_gas_limit(voids):
    # Uniform 300 K and 0.4 MPa, where helium is ideal to 0.3%. Pores and heat-exchanger voids swing isothermally with
    # the second virial coefficient; the dead volume, tubes, and reservoirs swing adiabatically from the mean state as
    # an ideal gas with gamma = 5/3, each reservoir by the swing reduced by sqrt(1 + (omega R B)^2) for
    # B = V_r / (gamma P).
    mean_pressure, swing = 0.4e6, 0.4e6
    resistance = (1.8e7, 3.6e7)  # Pa s/m^3, omega R B near 1 and 2
    mean_density = 0.6407114641  # HEPAK at 300 K and 0.4 MPa

    def isothermal(volume):
        def virial_density(pressure):
            return _MOLAR_MASS / (_GAS_CONSTANT * 300.0 / pressure + _SECOND_VIRIAL_AT_300_K)

        return volume * (virial_density(mean_pressure + swing / 2) - virial_density(mean_pressure - swing / 2))

    def adiabatic(volume, swing=swing):
        upper, lower = (1 + swing / (2 * mean_pressure)), (1 - swing / (2 * mean_pressure))
        return volume * mean_density * (upper**0.6 - lower**0.6)

    response = 2 * math.pi * _FREQUENCY * np.array(resistance) * 5.0e-3 / (5 / 3 * mean_pressure)
    reservoir = [adiabatic(5.0e-3, swing / math.sqrt(1 + r**2)) for r in response]
    expected = [
        adiabatic(4.5e-5),
        isothermal(1.0e-4),
        adiabatic(6.0e-5) + reservoir[0] + isothermal(voids[0]),
        isothermal(3.0e-5),
        isothermal(2.0e-5),
        isothermal(2.0e-5),
        adiabatic(2.0e-4) + reservoir[1] + isothermal(voids[1]),
    ]

    head = Head(
        model=_model(heat_exchanger_voids=voids), packing=(0.0, 0.0), resistance=resistance, frequency=_FREQUENCY
    )

    mass_swing = calc_mass_swing(head, np.full(7, 300.0), mean_pressure, swing)

    np.testing.assert_allclose(mass_swing, expected, rtol=3e-3)


@pytest.mark.parametrize("voids", [_NO_VOIDS, _VOIDS])
def test_mass_swing_with_cold_pores_tubes_and_voids(voids):
    # Without reservoirs: pores swing at their node temperatures and heat-exchanger voids at their stage temperatures,
    # and each tube's gas adiabatically from the mean state along a linear profile from the warm flange to its stage,
    # averaged by the default quadrature to within 0.1% of adaptive quadrature.
    mean_pressure, swing = 1.4e6, 1.2e6
    upper, lower = mean_pressure + swing / 2, mean_pressure - swing / 2

    def adiabatic(t):
        return float(
            helium.isentropic_density(t, mean_pressure, upper) - helium.isentropic_density(t, mean_pressure, lower)
        )

    expected = np.zeros(7)
    expected[0] = 4.5e-5 * adiabatic(300.0)
    for index, volume in _PORE.items():
        expected[index] = volume * float(helium.density(_COLD[index], upper) - helium.density(_COLD[index], lower))
    for (index, volume), void in zip(_TUBE.items(), voids, strict=True):
        expected[index] = volume * _mean_over_tube(adiabatic, _COLD[0], _COLD[index]) + void * float(
            helium.density(_COLD[index], upper) - helium.density(_COLD[index], lower)
        )

    head = Head(
        model=_model(reservoir_volume=0.0, heat_exchanger_voids=voids),
        packing=(0.0, 0.0),
        resistance=(1e9, 1e9),
        frequency=_FREQUENCY,
    )

    mass_swing = calc_mass_swing(head, _COLD, mean_pressure, swing)

    np.testing.assert_allclose(mass_swing, expected, rtol=1e-3)


def test_acoustic_power_follows_reservoir_impedance():
    # E = kappa^2 dP^2 / 4 * omega^2 R B^2 / (2 [1 + (omega R B)^2]) for B = V_r / (gamma P), with gamma from HEPAK
    # at the warm flange's temperature and the mean pressure.
    mean_pressure, swing, warm = 1.4e6, 1.2e6, 300.0
    resistance = (2.0e7, 6.0e7)
    head = Head(model=_model(), packing=(0.0, 0.0), resistance=resistance, frequency=_FREQUENCY)
    omega = 2 * math.pi * _FREQUENCY
    response = 5.0e-3 / (float(helium.heat_capacity_ratio(warm, mean_pressure)) * mean_pressure)
    amplitude = math.sqrt(2) * swing / 2
    expected = [amplitude**2 * omega**2 * r * response**2 / (2 * (1 + (omega * r * response) ** 2)) for r in resistance]

    np.testing.assert_allclose(calc_acoustic_power(head, mean_pressure, swing, warm), expected, rtol=1e-10)


def test_acoustic_power_with_stiff_reservoir_is_amplitude_squared_over_twice_resistance():
    # With omega R B near 4000, E = P_a^2 / (2 R) for P_a = sqrt(2) * 1.2 MPa / 2, to 1 part in 10^7.
    head = Head(model=_model(), packing=(0.0, 0.0), resistance=(1e11, 4e11), frequency=_FREQUENCY)

    power = calc_acoustic_power(head, 1.4e6, 1.2e6, 300.0)

    np.testing.assert_allclose(power, [0.72e12 / 2e11, 0.72e12 / 8e11], rtol=1e-6)


def test_heat_flow_without_swing_is_wall_conduction():
    # Walls of k = 2 W/(m K). Regenerator walls (A = 1e-4 m^2) conduct along a piecewise-linear profile through the
    # node centres: stage-1 section centre 0.1 m from either end; stage-2 section centres 0.05, 0.125, and 0.175 m
    # from stage 1 in its 0.2 m. Each pulse tube (A = 1e-5 m^2, 0.3 m) conducts from the warm flange to its stage.
    #   stage-1 section: 2e-4 * (300 - 60) / 0.1 = 0.48 W in, 2e-4 * (60 - 40) / 0.1 = 0.04 W out
    #   stage-2 section: 2e-4 * [(40 - 20) / 0.05, (20 - 10) / 0.075, (10 - 6) / 0.05, (6 - 4) / 0.025]
    #                    = 0.08, 0.02667, 0.016, 0.016 W
    #   tubes:           2e-5 / 0.3 * [260, 296] = 0.017333, 0.019733 W
    # A stroke spread of 1% of each tube's gas leaves no breakthrough while the gas stands still.
    head = Head(
        model=dataclasses.replace(_model(), stroke_spread=0.01),
        packing=(0.1, 0.1),
        resistance=(2e7, 6e7),
        frequency=_FREQUENCY,
    )
    expected = [
        0.0,
        0.48 - 0.04,
        0.04 + 0.0173333 - 0.08,
        0.08 - 0.0266667,
        0.0266667 - 0.016,
        0.016 - 0.016,
        0.016 + 0.0197333,
    ]

    np.testing.assert_allclose(calc_heat_flow(head, _COLD, 1.4e6, 0.0), expected, atol=1e-6)


def _capacity_ratio(temperature, pressure, pore_volume, matrix_mass):
    """Pore gas over matrix heat capacity, for the test matrix's specific heat of T J/(kg K)."""
    gas = pore_volume * float(
        helium.density(temperature, pressure) * helium.isobaric_specific_heat(temperature, pressure)
    )
    return gas / (matrix_mass * temperature)


def test_heat_flow_at_uniform_temperature_is_pressure_transport():
    # Without gradients only F = (1 - beta T / (1 + epsilon)) E flows, with E unchanged along each section. Each plate
    # boundary takes the epsilon of its node. Stage 2 nets -beta T E_2 / (1 + epsilon), stage 1 the stage-1 section's
    # cold-end flow of E_1 + E_2 less E_1 and the stage-2 section's warm-end flow, and each section's nodes together
    # the difference between its ends' flows.
    temperature, mean_pressure, swing = 10.0, 1.4e6, 1.2e6
    head = Head(model=_model(), packing=(0.1, 0.1), resistance=(2e7, 6e7), frequency=_FREQUENCY)
    first, second = calc_acoustic_power(head, mean_pressure, swing, temperature)
    expansion = temperature * float(helium.isobaric_expansion_coefficient(temperature, mean_pressure))
    stage_1_section = _capacity_ratio(temperature, mean_pressure, 1.0e-4, 1.0)
    stage_2_warm = _capacity_ratio(temperature, mean_pressure, 3.0e-5, 1.0)
    stage_2_cold = _capacity_ratio(temperature, mean_pressure, 2.0e-5, 0.5)

    def flow(capacity_ratio, power):
        return (1 - expansion / (1 + capacity_ratio)) * power

    heat = calc_heat_flow(head, np.full(7, temperature), mean_pressure, swing)

    assert heat[0] == 0.0
    assert heat[1] == pytest.approx(0.0, abs=1e-9 * first)
    assert heat[2] == pytest.approx(
        flow(stage_1_section, first + second) - first - flow(stage_2_warm, second), rel=1e-6
    )
    assert heat[3:6].sum() == pytest.approx(flow(stage_2_warm, second) - flow(stage_2_cold, second), rel=1e-6)
    assert heat[6] == pytest.approx(-expansion * second / (1 + stage_2_cold), rel=1e-6)


def _isothermal_swing(volume, temperature, mean_pressure, swing):
    """Mass in kg a volume at `temperature` takes in over the swing."""
    return volume * float(
        helium.density(temperature, mean_pressure + swing / 2) - helium.density(temperature, mean_pressure - swing / 2)
    )


@pytest.mark.parametrize("voids", [_NO_VOIDS, _VOIDS])
def test_regeneration_at_stage_2_follows_packing(voids):
    # The packing coefficient adds a (m c_p)^2 / (f C'_m) dT/dx at the stage-2 boundary, with m the mass flow of tube
    # 2, reservoir 2, and the stage-2 heat exchanger's void at 4 K, and c_p at T_2. C'_m is half the coldest node,
    # 0.25 kg, at the test matrix's mean specific heat between its centre's 6 K and stage 2's 4 K,
    # (6^2 - 4^2) / 2 / (6 - 4) = 5 J/(kg K), over the 0.025 m between them.
    mean_pressure, swing, packing = 1.4e6, 1.2e6, 0.3
    bare = Head(model=_model(), packing=(0.1, 0.0), resistance=(2e7, 6e7), frequency=_FREQUENCY)
    loose = dataclasses.replace(bare, model=_model(heat_exchanger_voids=voids))
    packed = dataclasses.replace(loose, packing=(0.1, packing))
    mass_flow = _FREQUENCY * (
        calc_mass_swing(bare, _COLD, mean_pressure, swing)[6] + _isothermal_swing(voids[1], 4.0, mean_pressure, swing)
    )
    capacity_rate = mass_flow * float(helium.isobaric_specific_heat(4.0, mean_pressure))
    expected = packing * capacity_rate**2 / (_FREQUENCY * 0.25 * 5.0 / 0.025) * (6.0 - 4.0) / 0.025

    change = calc_heat_flow(packed, _COLD, mean_pressure, swing) - calc_heat_flow(loose, _COLD, mean_pressure, swing)

    assert change[6] == pytest.approx(expected, rel=1e-9)


@pytest.mark.parametrize("voids", [_NO_VOIDS, _VOIDS])
def test_regeneration_at_stage_1_carries_every_cold_side_volume(voids):
    # At the stage-1 boundary, m is the mass flow of everything from stage 1 down: the void-free head's swing there
    # and beyond, plus both heat-exchanger voids at their stage temperatures. c_p is at T_1 = 40 K. C'_m is half the
    # stage-1 section's node, 0.5 kg, at the mean specific heat between its centre's 60 K and stage 1's 40 K,
    # (60^2 - 40^2) / 2 / (60 - 40) = 50 J/(kg K), over the 0.1 m between them.
    mean_pressure, swing, packing = 1.4e6, 1.2e6, 0.3
    bare = Head(model=_model(), packing=(0.0, 0.1), resistance=(2e7, 6e7), frequency=_FREQUENCY)
    loose = dataclasses.replace(bare, model=_model(heat_exchanger_voids=voids))
    packed = dataclasses.replace(loose, packing=(packing, 0.1))
    mass_flow = _FREQUENCY * (
        calc_mass_swing(bare, _COLD, mean_pressure, swing)[2:].sum()
        + _isothermal_swing(voids[0], 40.0, mean_pressure, swing)
        + _isothermal_swing(voids[1], 4.0, mean_pressure, swing)
    )
    capacity_rate = mass_flow * float(helium.isobaric_specific_heat(40.0, mean_pressure))
    expected = packing * capacity_rate**2 / (_FREQUENCY * 0.5 * 50.0 / 0.1) * (60.0 - 40.0) / 0.1

    change = calc_heat_flow(packed, _COLD, mean_pressure, swing) - calc_heat_flow(loose, _COLD, mean_pressure, swing)

    assert change[2] == pytest.approx(expected, rel=1e-9)


def test_regeneration_between_layers_takes_each_half_nodes_mean_specific_heat():
    # The stage-2 section's boundary 0.1 m from stage 1 lies between the first layer's node, 1 kg centred 0.05 m in at
    # 20 K, and the second layer's first node, 0.5 kg centred at 0.125 m at 10 K. C'_m is both half-nodes,
    # 0.5 + 0.25 kg, at the mean specific heat between 20 K and 10 K, (20^2 - 10^2) / 2 / (20 - 10) = 15 J/(kg K),
    # over the 0.075 m between centres. m is the mass flow of everything on the cold side, and c_p is at the boundary's
    # 20 - 10 * 0.05 / 0.075 K. Packing only the stage-2 section leaves stage 1's flow from its own section unchanged,
    # so stage 1 and the 20 K node together lose this boundary's flow.
    mean_pressure, swing, packing = 1.4e6, 1.2e6, 0.3
    loose = Head(model=_model(), packing=(0.1, 0.0), resistance=(2e7, 6e7), frequency=_FREQUENCY)
    packed = dataclasses.replace(loose, packing=(0.1, packing))
    mass_flow = _FREQUENCY * calc_mass_swing(loose, _COLD, mean_pressure, swing)[4:].sum()
    capacity_rate = mass_flow * float(helium.isobaric_specific_heat(20.0 - 10.0 * 0.05 / 0.075, mean_pressure))
    expected = packing * capacity_rate**2 / (_FREQUENCY * 0.75 * 15.0 / 0.075) * (20.0 - 10.0) / 0.075

    change = calc_heat_flow(packed, _COLD, mean_pressure, swing) - calc_heat_flow(loose, _COLD, mean_pressure, swing)

    assert change[2] + change[3] == pytest.approx(-expected, rel=1e-9)


def test_pressure_transport_at_stage_2_takes_the_mean_specific_heat():
    # Without packing, stage 2 takes -beta T E_2 / (1 + epsilon) from pressure transport at its boundary, and nothing
    # else there depends on the matrix. Against a matrix a trillion times heavier, where epsilon vanishes, it gains
    # beta T E_2 epsilon / (1 + epsilon), with beta T at 4 K. epsilon is half the coldest node's pore gas, 1e-5 m^3 at
    # 6 K, over half its matrix, 0.25 kg at the mean specific heat of 5 J/(kg K) between 6 K and 4 K.
    mean_pressure, swing = 1.4e6, 1.2e6
    head = Head(model=_model(), packing=(0.0, 0.0), resistance=(2e7, 6e7), frequency=_FREQUENCY)
    heavy_matrix = dataclasses.replace(
        _MATRIX, specific_heat=lambda t: 1e12 * t, specific_heat_integral=lambda t: 5e11 * t**2
    )
    heavy = dataclasses.replace(head, model=_model(matrix=heavy_matrix))
    _, power = calc_acoustic_power(head, mean_pressure, swing, _COLD[0])
    expansion = 4.0 * float(helium.isobaric_expansion_coefficient(4.0, mean_pressure))
    gas = 1.0e-5 * float(helium.density(6.0, mean_pressure) * helium.isobaric_specific_heat(6.0, mean_pressure))
    ratio = gas / (0.25 * 5.0)
    expected = expansion * power * ratio / (1 + ratio)

    change = calc_heat_flow(head, _COLD, mean_pressure, swing) - calc_heat_flow(heavy, _COLD, mean_pressure, swing)

    assert change[6] == pytest.approx(expected, rel=1e-9)


def test_heat_flow_with_constant_specific_heat_is_unchanged_by_boundary_averaging():
    # A matrix of constant specific heat, 100 J/(kg K), has the same capacity at a boundary whether its specific heat
    # is taken at the node centres or averaged between them. Heat flows recorded from the point-sampled implementation.
    constant = dataclasses.replace(
        _MATRIX,
        specific_heat=lambda t: np.full_like(np.asarray(t, dtype=float), 100.0),
        specific_heat_integral=lambda t: 100.0 * np.asarray(t),
    )
    head = Head(model=_model(matrix=constant), packing=(0.1, 0.3), resistance=(2e7, 6e7), frequency=_FREQUENCY)
    expected = [
        0.0,
        10894.92186635623,
        7094.392734657332,
        3103.710953063086,
        -2112.114870413167,
        -1921.017924767319,
        5443.247966596633,
    ]

    np.testing.assert_allclose(calc_heat_flow(head, _COLD, 1.4e6, 1.2e6), expected, rtol=1e-9)


def _tube_gas_mass(index, pressure):
    """Gas mass in kg of the pulse tube at a stage's chain index, by adaptive quadrature over its profile."""
    return _TUBE[index] * _mean_over_tube(lambda t: float(helium.density(t, pressure)), _COLD[0], _COLD[index])


def _enthalpy_drop(index, pressure):
    """Specific enthalpy in J/kg of helium at the warm flange less that at the stage at a chain index."""
    return float(helium.enthalpy(_COLD[0], pressure) - helium.enthalpy(_COLD[index], pressure))


def test_breakthrough_loss_vanishes_while_stroke_stays_within_buffer():
    # A 0.2 MPa swing with stiff reservoirs moves 2% and 0.1% of each tube's gas through its warm end, at least four
    # spreads of 10% short of its buffer of half the tube's gas. Each loss is below 1e-4 of the enthalpy one spread of
    # gas carries per second, f (h(T_0) - h(T_i)) sigma m_i.
    mean_pressure = 1.4e6
    head = Head(model=_model(), packing=(0.0, 0.0), resistance=_STIFF, frequency=_FREQUENCY)

    loss = calc_breakthrough_loss(head, _COLD, mean_pressure, 0.2e6)

    for value, index in zip(loss, _TUBE, strict=True):
        spread_enthalpy = _FREQUENCY * _enthalpy_drop(index, mean_pressure) * 0.1 * _tube_gas_mass(index, mean_pressure)
        assert 0.0 <= value < 1e-4 * spread_enthalpy


def test_breakthrough_loss_ignores_the_tubes_own_compression():
    # Without reservoirs, a 2.0 MPa swing compresses 89% and 70% of each tube's gas through its cold end, but none
    # passes the warm end, so nothing breaks through: each loss is below 1e-4 of f (h(T_0) - h(T_i)) sigma m_i.
    mean_pressure = 1.4e6
    head = Head(model=_model(reservoir_volume=0.0), packing=(0.0, 0.0), resistance=_OPEN, frequency=_FREQUENCY)

    loss = calc_breakthrough_loss(head, _COLD, mean_pressure, 2.0e6)

    for value, index in zip(loss, _TUBE, strict=True):
        spread_enthalpy = _FREQUENCY * _enthalpy_drop(index, mean_pressure) * 0.1 * _tube_gas_mass(index, mean_pressure)
        assert 0.0 <= value < 1e-4 * spread_enthalpy


def test_breakthrough_loss_far_past_buffer_carries_excess_enthalpy():
    # With open valves a 1.2 MPa swing moves 18 and 2.0 times each tube's gas through its warm end, many spreads past
    # its buffer. Each cycle the excess over half the tube's gas carries the enthalpy difference between the warm
    # flange and the stage: L_i = f (h(T_0) - h(T_i)) (dm_i - m_i / 2), with dm_i the reservoir's swing, the head's
    # swing at the stage less that of the same head without reservoirs, and m_i averaged by the default quadrature to
    # within 1e-4 of adaptive quadrature.
    mean_pressure, swing = 1.4e6, 1.2e6
    head = Head(model=_model(), packing=(0.0, 0.0), resistance=_OPEN, frequency=_FREQUENCY)
    closed = dataclasses.replace(head, model=_model(reservoir_volume=0.0))
    reservoir_swing = calc_mass_swing(head, _COLD, mean_pressure, swing) - calc_mass_swing(
        closed, _COLD, mean_pressure, swing
    )
    expected = [
        _FREQUENCY
        * _enthalpy_drop(index, mean_pressure)
        * (reservoir_swing[index] - 0.5 * _tube_gas_mass(index, mean_pressure))
        for index in _TUBE
    ]

    np.testing.assert_allclose(calc_breakthrough_loss(head, _COLD, mean_pressure, swing), expected, rtol=1e-4)


def test_breakthrough_loss_rises_smoothly_past_buffer():
    # Swings from 0.2 to 2.0 MPa carry each tube's stroke through its warm end from 10% and 9% of its gas, within its
    # buffer, to 96% and 93%, past it. The loss rises throughout with continuous derivatives, as in the heat-flow
    # roughness test below.
    head = Head(model=_model(), packing=(0.0, 0.0), resistance=(2e9, 3e8), frequency=_FREQUENCY)

    def losses(points):
        return np.array(
            [calc_breakthrough_loss(head, _COLD, 1.4e6, swing) for swing in np.linspace(0.2e6, 2.0e6, points)]
        )

    def roughness(loss):
        return np.max(np.abs(np.diff(loss, 2, axis=0)), axis=0) / np.max(np.abs(np.diff(loss, axis=0)), axis=0)

    coarse, fine = losses(41), losses(161)

    assert np.all(np.diff(fine, axis=0) > 0)
    assert np.all(fine[0] < 1e-3 * fine[-1])
    assert np.all(roughness(fine) <= 0.5 * roughness(coarse))


def test_heat_flow_takes_breakthrough_loss_at_each_stage():
    # Against a head whose buffers hold a hundred times each tube's gas, the heat flow gains L_1 at stage 1 and L_2 at
    # stage 2, and nothing elsewhere.
    mean_pressure, swing = 1.4e6, 1.2e6
    head = Head(model=_model(), packing=(0.1, 0.1), resistance=_OPEN, frequency=_FREQUENCY)
    buffered = dataclasses.replace(head, model=dataclasses.replace(head.model, buffer_fraction=100.0))
    expected = np.zeros(7)
    expected[list(_TUBE)] = calc_breakthrough_loss(head, _COLD, mean_pressure, swing)

    change = calc_heat_flow(head, _COLD, mean_pressure, swing) - calc_heat_flow(buffered, _COLD, mean_pressure, swing)

    np.testing.assert_allclose(change, expected, rtol=1e-9, atol=1e-12 * np.max(expected))


def test_stage_coupling_loads_stage_1_with_stage_2s_gross_cooling():
    # Against the uncoupled head, stage 1 gains c E_2 x (1 - x) and nothing else changes. x = beta T_2 / (1 + epsilon)
    # at T_2 = 4 K, with epsilon the coldest node's: pore gas of 2e-5 m^3 at 6 K over 0.5 kg of matrix at a constant
    # 5 J/(kg K).
    mean_pressure, swing, coupling = 1.4e6, 1.2e6, 2.0
    constant = dataclasses.replace(
        _MATRIX, specific_heat=lambda t: np.full(np.shape(t), 5.0), specific_heat_integral=lambda t: 5.0 * np.asarray(t)
    )
    model = _model()
    model = dataclasses.replace(
        model,
        regenerators=tuple(
            dataclasses.replace(
                section, layers=tuple(dataclasses.replace(layer, material=constant) for layer in section.layers)
            )
            for section in model.regenerators
        ),
    )
    uncoupled = Head(model=model, packing=(0.1, 0.1), resistance=(2e7, 6e7), frequency=_FREQUENCY)
    coupled = dataclasses.replace(uncoupled, model=dataclasses.replace(model, stage_coupling=coupling))
    _, stage_2_power = calc_acoustic_power(uncoupled, mean_pressure, swing, _COLD[0])
    expansion = 4.0 * float(helium.isobaric_expansion_coefficient(4.0, mean_pressure))
    gas = 2.0e-5 * float(helium.density(6.0, mean_pressure) * helium.isobaric_specific_heat(6.0, mean_pressure))
    fraction = expansion / (1 + gas / (0.5 * 5.0))
    expected = np.zeros(7)
    expected[2] = coupling * stage_2_power * fraction * (1 - fraction)

    change = calc_heat_flow(coupled, _COLD, mean_pressure, swing) - calc_heat_flow(
        uncoupled, _COLD, mean_pressure, swing
    )

    np.testing.assert_allclose(change, expected, rtol=1e-9, atol=1e-12 * abs(expected[2]))


# PT1.5 with placeholder packing coefficients and valve resistances, its matrix given a smooth
# specific heat of 0.1 T^3 + 20 J/(kg K).
_SMOOTH_MATRIX = dataclasses.replace(
    _MATRIX,
    specific_heat=lambda t: 0.1 * np.asarray(t) ** 3 + 20.0,
    specific_heat_integral=lambda t: 0.025 * np.asarray(t) ** 4 + 20.0 * np.asarray(t),
)
_PT1_5 = load_head_model("PT1.5")
_SMOOTH_PT1_5 = Head(
    model=dataclasses.replace(
        _PT1_5,
        regenerators=tuple(
            dataclasses.replace(
                section,
                layers=tuple(dataclasses.replace(layer, material=_SMOOTH_MATRIX) for layer in section.layers),
            )
            for section in _PT1_5.regenerators
        ),
    ),
    packing=(0.15, 0.05),
    resistance=(1e10, 3e10),
    frequency=1.4,
)


def test_heat_flow_has_continuous_derivatives_through_a_cooldown():
    # Stages cooling geometrically from 295 K to 40 K and 4 K at a fixed drive. Along a path with continuous
    # derivatives, second differences shrink against first differences in proportion to the step; at a kink they do not.
    # The warm flange carries no heat flow.
    def roughness(points):
        heat = []
        for progress in np.linspace(0.0, 1.0, points):
            stages = (295.0 * (40.0 / 295.0) ** progress, 295.0 * (4.0 / 295.0) ** progress)
            temperature = _linear_chain(_PT1_5, stages, warm_temperature=295.0)
            heat.append(calc_heat_flow(_SMOOTH_PT1_5, temperature, 1.3e6, 1.2e6)[1:])
        heat = np.array(heat)
        return np.max(np.abs(np.diff(heat, 2, axis=0)), axis=0) / np.max(np.abs(np.diff(heat, axis=0)), axis=0)

    coarse, fine = roughness(41), roughness(161)

    assert np.all(fine <= 0.5 * coarse)


def _linear_chain(model, stage_temperature, *, warm_temperature):
    """Chain temperatures with each section's nodes linear in position between its ends' temperatures."""
    ends = (warm_temperature, *stage_temperature)
    sections = []
    for section, warm, cold in zip(model.regenerators, ends[:-1], ends[1:], strict=True):
        fractions = np.array(
            [layer.length_fraction / layer.node_count for layer in section.layers for _ in range(layer.node_count)]
        )
        sections.append(warm + (cold - warm) * (np.cumsum(fractions) - fractions / 2))
    return np.concatenate(([warm_temperature], sections[0], [ends[1]], sections[1], [ends[2]]))
