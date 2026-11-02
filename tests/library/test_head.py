"""Tests for the network equations of one pulse-tube head."""

import numpy as np
import pytest
from scipy.constants import gas_constant

from pulse_tube_cooldown.library.links.cooler.head import (
    calc_acoustic_power,
    calc_capacity_ratio,
    calc_enthalpy_flow,
    calc_gas_mass,
    calc_heat_flow,
    calc_mass_flow,
    calc_utilization,
)

_MOLAR_MASS = 4.002602e-3  # kg/mol

# A head with one stage-1 layer and two stage-2 layers, in chain order: warm flange, stage-1 column, stage 1,
# stage-2 column, stage 2.
_TEMPERATURE = np.array([300.0, 60.0, 40.0, 20.0, 8.0, 4.0])
_PORE_VOLUME = (np.array([1.0e-4]), np.array([3.0e-5, 2.0e-5]))
_PULSE_TUBE_VOLUME = (6.0e-5, 2.0e-4)
_HEAT_CAPACITY = (np.array([800.0]), np.array([40.0, 15.0]))
_SUPPLY_PRESSURE = 2.0e6
_RETURN_PRESSURE = 0.8e6
_MEAN_PRESSURE = (_SUPPLY_PRESSURE + _RETURN_PRESSURE) / 2
_FREQUENCY = 1.2
_PACKING_CONSTANT = (0.8, 1.3)
_ENTHALPY_FLOW = (30.0, 1.5)
# Pore-gas to matrix heat-capacity ratios of the three nodes, from `test_capacity_ratio_matches_hand_calculation`.
_CAPACITY_RATIO = (np.array([7.163822e-3]), np.array([0.1480389, 0.9654703]))

# Hand computation from HEPAK at these states and the ADR's network form. Mass swing rho(T, P_h) - rho(T, P_l) times
# volume, in kg: pores 9.058225e-4 (60 K), 8.335998e-4 (20 K), 1.131456e-3 (8 K); pulse tubes by adaptive quadrature
# over a linear profile from 300 K, 2.581565e-4 (to 40 K) and 1.457110e-3 (to 4 K). Boundary flows are f times the
# swing beyond each boundary.
_MASS_FLOW = np.array([5.503374e-3, 4.416387e-3, 4.106599e-3, 3.106279e-3, 1.748532e-3])
# Boundary conduction 2 a (m c_p)^2 (T_w - T_c) / (f (C_w + C_c)), c_p at 1.4 MPa and the boundary's mean temperature:
# 326.9320, 18.23414, 553.1576, 202.8832, 30.42731 W. Stage-2 boundaries at 40 K, 14 K, and 4 K, with HEPAK densities
# 16.24194, 51.76106, 158.0951 kg/m^3 and beta T of 0.9933046, 1.349434, 0.1623257 at 1.4 MPa. The 20 K node spans
# 40 K to 14 K and the 8 K node 14 K to 4 K, so the acoustic power at those boundaries is H_2 times
# (158.0951 / 51.76106)^(1 / 1.9654703) times (51.76106 / 16.24194)^(1 / 1.1480389), H_2 times the first factor, and
# H_2: 7.265492, 2.647332, 1.5 W. With eps of 0.1480389, 0.5567546, and 0.9654703 at the boundaries, the flows
# (1 - beta T / (1 + eps)) E_2 toward stage 2 are 0.9792531, 0.3525575, 1.376117 W. Stage sinks H_1 and H_2 land on the
# warm flange.
_ACOUSTIC_POWER = [7.265492, 2.647332, 1.5]
_HEAT_FLOW = [-295.4320, 308.6979, -565.9027, 350.9011, 171.4323, 30.30343]


def test_mass_flow_matches_hand_calculation():
    mass_flow = calc_mass_flow(
        _TEMPERATURE,
        pore_volume=_PORE_VOLUME,
        pulse_tube_volume=_PULSE_TUBE_VOLUME,
        supply_pressure=_SUPPLY_PRESSURE,
        return_pressure=_RETURN_PRESSURE,
        frequency=_FREQUENCY,
    )

    np.testing.assert_allclose(mass_flow, _MASS_FLOW, rtol=1e-6)


def test_mass_flow_of_head_at_one_temperature():
    # HEPAK's density swing at 300 K between 2.0 and 0.8 MPa is 1.901227 kg/m^3, the same for pores and pulse tubes.
    # Gas volume beyond each boundary of the test head, in m^3.
    volume_beyond = np.array([4.1e-4, 3.1e-4, 2.5e-4, 2.2e-4, 2.0e-4])

    mass_flow = calc_mass_flow(
        np.full(6, 300.0),
        pore_volume=_PORE_VOLUME,
        pulse_tube_volume=_PULSE_TUBE_VOLUME,
        supply_pressure=_SUPPLY_PRESSURE,
        return_pressure=_RETURN_PRESSURE,
        frequency=_FREQUENCY,
    )

    np.testing.assert_allclose(mass_flow, _FREQUENCY * 1.901227 * volume_beyond, rtol=1e-6)


def test_gas_mass_matches_hand_calculation():
    # HEPAK density at 1.4 MPa times volume for the pores, and by adaptive quadrature over a linear profile from 300 K
    # for the pulse tubes.
    gas_mass = calc_gas_mass(
        _TEMPERATURE, pore_volume=_PORE_VOLUME, pulse_tube_volume=_PULSE_TUBE_VOLUME, pressure=_MEAN_PRESSURE
    )

    assert gas_mass == pytest.approx(6.7574036e-3, rel=1e-6)


def test_gas_mass_of_head_at_one_temperature():
    # HEPAK's density at 300 K and 1.4 MPa, 2.2322996 kg/m^3, times the head's 4.1e-4 m^3 of gas volume.
    gas_mass = calc_gas_mass(
        np.full(6, 300.0), pore_volume=_PORE_VOLUME, pulse_tube_volume=_PULSE_TUBE_VOLUME, pressure=_MEAN_PRESSURE
    )

    assert gas_mass == pytest.approx(2.2322996 * 4.1e-4, rel=1e-6)


def test_acoustic_power_matches_hand_calculation():
    acoustic_power = calc_acoustic_power(
        _TEMPERATURE, capacity_ratio=_CAPACITY_RATIO, enthalpy_flow=_ENTHALPY_FLOW[1], mean_pressure=_MEAN_PRESSURE
    )

    np.testing.assert_allclose(acoustic_power, _ACOUSTIC_POWER, rtol=1e-6)


def test_acoustic_power_without_pore_gas_scales_with_density():
    # With eps = 0 the product of density ratios telescopes to rho(T_2) / rho(T_b): 158.0951 / 16.24194 at 40 K and
    # 158.0951 / 51.76106 at 14 K.
    acoustic_power = calc_acoustic_power(
        _TEMPERATURE,
        capacity_ratio=(np.zeros(1), np.zeros(2)),
        enthalpy_flow=_ENTHALPY_FLOW[1],
        mean_pressure=_MEAN_PRESSURE,
    )

    np.testing.assert_allclose(acoustic_power, [1.5 * 9.733758, 1.5 * 3.054327, 1.5], rtol=1e-6)


def test_heat_flow_matches_hand_calculation():
    heat_flow = calc_heat_flow(
        _TEMPERATURE,
        heat_capacity=_HEAT_CAPACITY,
        capacity_ratio=_CAPACITY_RATIO,
        mass_flow=_MASS_FLOW,
        enthalpy_flow=_ENTHALPY_FLOW,
        packing_constant=_PACKING_CONSTANT,
        mean_pressure=_MEAN_PRESSURE,
        frequency=_FREQUENCY,
    )

    np.testing.assert_allclose(heat_flow, _HEAT_FLOW, rtol=1e-6)
    assert heat_flow.sum() == pytest.approx(0.0, abs=1e-9)


def test_utilization_matches_hand_calculation():
    # Cold-end flow times c_p at the stage temperature and 1.4 MPa, over f times the column's summed heat capacity.
    utilization = calc_utilization(
        _TEMPERATURE,
        heat_capacity=_HEAT_CAPACITY,
        mass_flow=_MASS_FLOW,
        mean_pressure=_MEAN_PRESSURE,
        frequency=_FREQUENCY,
    )

    np.testing.assert_allclose(utilization, [0.02466503, 0.07090502], rtol=1e-6)


def test_capacity_ratio_matches_hand_calculation():
    # Pore volume times HEPAK density and c_p at 1.4 MPa, over the node's matrix heat capacity. At 60 K, 20 K, and
    # 8 K: densities 10.89474, 33.43007, 114.8121 kg/m^3 and c_p 5260.389, 5904.419, 6306.848 J/(kg K).
    first, second = calc_capacity_ratio(
        _TEMPERATURE, heat_capacity=_HEAT_CAPACITY, pore_volume=_PORE_VOLUME, mean_pressure=_MEAN_PRESSURE
    )

    np.testing.assert_allclose(first, [7.163822e-3], rtol=1e-6)
    np.testing.assert_allclose(second, [0.1480389, 0.9654703], rtol=1e-6)


def test_capacity_ratio_in_ideal_gas_limit():
    # Near 300 K helium's volumetric heat capacity rho c_p approaches the ideal monatomic 5 P / (2 T).
    temperature = np.array([300.0, 290.0, 280.0, 270.0, 260.0, 250.0])
    pressure = 1.0e6

    first, second = calc_capacity_ratio(
        temperature, heat_capacity=_HEAT_CAPACITY, pore_volume=_PORE_VOLUME, mean_pressure=pressure
    )

    ideal = 2.5 * pressure / temperature
    np.testing.assert_allclose(first, _PORE_VOLUME[0] * ideal[[1]] / _HEAT_CAPACITY[0], rtol=1e-2)
    np.testing.assert_allclose(second, _PORE_VOLUME[1] * ideal[[3, 4]] / _HEAT_CAPACITY[1], rtol=1e-2)


def test_enthalpy_flow_from_pressure_swing():
    # H = C (kappa dP / 2)^2 / 2 = 1e-10 * (0.9 * 1e6 / 2)^2 / 2.
    assert calc_enthalpy_flow(1.0e-10, 1.0e6, 0.9) == pytest.approx(10.125)


def test_ideal_gas_limit_reproduces_section_equations():
    # Warm enough for helium to be near ideal, with no pore gas, so each column carries one mass flow. With each
    # column's single node at its midpoint, the node is in balance and each column carries the section loss
    # L = a (m c_p)^2 dT / (f M c_m).
    warm, stage_1, stage_2 = 300.0, 200.0, 150.0
    temperature = np.array([warm, 250.0, stage_1, 175.0, stage_2])
    supply_pressure, return_pressure, frequency = 1.2e6, 0.6e6, 1.0
    tube_1, tube_2 = 6.0e-5, 2.0e-4
    matrix_1, matrix_2 = 300.0, 100.0
    packing_constant = (0.8, 1.3)
    enthalpy_flow = (30.0, 0.5)

    def tube_swing(volume, cold):
        # Ideal-gas density integrated over a linear profile from the warm flange.
        density_swing_times_temperature = _MOLAR_MASS * (supply_pressure - return_pressure) / gas_constant
        return volume * density_swing_times_temperature * np.log(warm / cold) / (warm - cold)

    flow_2 = frequency * tube_swing(tube_2, stage_2)
    flow_1 = flow_2 + frequency * tube_swing(tube_1, stage_1)
    specific_heat = 2.5 * gas_constant / _MOLAR_MASS
    loss_1 = packing_constant[0] * (flow_1 * specific_heat) ** 2 * (warm - stage_1) / (frequency * matrix_1)
    loss_2 = packing_constant[1] * (flow_2 * specific_heat) ** 2 * (stage_1 - stage_2) / (frequency * matrix_2)

    mass_flow = calc_mass_flow(
        temperature,
        pore_volume=(np.zeros(1), np.zeros(1)),
        pulse_tube_volume=(tube_1, tube_2),
        supply_pressure=supply_pressure,
        return_pressure=return_pressure,
        frequency=frequency,
    )
    heat_flow = calc_heat_flow(
        temperature,
        heat_capacity=(np.array([matrix_1]), np.array([matrix_2])),
        capacity_ratio=(np.zeros(1), np.zeros(1)),
        mass_flow=mass_flow,
        enthalpy_flow=enthalpy_flow,
        packing_constant=packing_constant,
        mean_pressure=(supply_pressure + return_pressure) / 2,
        frequency=frequency,
    )

    np.testing.assert_allclose(mass_flow, [flow_1, flow_1, flow_2, flow_2], rtol=2e-2)
    expected = [
        -loss_1 + sum(enthalpy_flow),
        0.0,
        loss_1 - loss_2 - enthalpy_flow[0],
        0.0,
        loss_2 - enthalpy_flow[1],
    ]
    np.testing.assert_allclose(heat_flow, expected, atol=3e-2 * max(loss_1, loss_2))


def test_ideal_gas_limit_of_finite_capacity_terms():
    # Warm enough for helium to be near ideal, so beta T = 1 and the real-gas flow vanishes, at a low pressure to keep
    # it so. With no packing loss, each stage-2 boundary carries only the finite-capacity flow eps / (1 + eps) E_2
    # toward stage 2, where E_2 scales from H_2 by (T_warm / T_cold)^(1 / (1 + eps)) across each node and eps is
    # V (5 P / 2 T) / C at a node.
    temperature = np.array([300.0, 250.0, 200.0, 180.0, 160.0, 150.0])
    supply_pressure, return_pressure = 0.4e6, 0.2e6
    mean_pressure = (supply_pressure + return_pressure) / 2
    heat_capacity = (np.array([300.0]), np.array([10.0, 5.0]))
    pore_volume = (np.zeros(1), np.array([1.0e-3, 1.0e-3]))
    enthalpy_flow = (30.0, 20.0)
    warm, middle = 2.5 * mean_pressure * pore_volume[1] / (temperature[[3, 4]] * heat_capacity[1])
    boundary_ratio = np.array([warm, (warm + middle) / 2, middle])
    acoustic_power = enthalpy_flow[1] * np.array(
        [(170 / 150) ** (1 / (1 + middle)) * (200 / 170) ** (1 / (1 + warm)), (170 / 150) ** (1 / (1 + middle)), 1.0]
    )
    flow = boundary_ratio / (1 + boundary_ratio) * acoustic_power

    capacity_ratio = calc_capacity_ratio(
        temperature, heat_capacity=heat_capacity, pore_volume=pore_volume, mean_pressure=mean_pressure
    )
    heat_flow = calc_heat_flow(
        temperature,
        heat_capacity=heat_capacity,
        capacity_ratio=capacity_ratio,
        mass_flow=calc_mass_flow(
            temperature,
            pore_volume=pore_volume,
            pulse_tube_volume=(6.0e-5, 2.0e-4),
            supply_pressure=supply_pressure,
            return_pressure=return_pressure,
            frequency=1.0,
        ),
        enthalpy_flow=enthalpy_flow,
        packing_constant=(0.0, 0.0),
        mean_pressure=mean_pressure,
        frequency=1.0,
    )

    np.testing.assert_allclose(
        calc_acoustic_power(temperature, capacity_ratio, enthalpy_flow[1], mean_pressure), acoustic_power, rtol=1e-2
    )
    expected = [
        sum(enthalpy_flow),
        0.0,
        -enthalpy_flow[0] - flow[0],
        flow[0] - flow[1],
        flow[1] - flow[2],
        flow[2] - enthalpy_flow[1],
    ]
    np.testing.assert_allclose(heat_flow, expected, atol=1e-2 * enthalpy_flow[1])
