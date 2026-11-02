"""Tests for the network equations of one pulse-tube head."""

import numpy as np
import pytest
from scipy.constants import gas_constant

from pulse_tube_cooldown.library.links.head import (
    calc_enthalpy_flow,
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

# Hand computation from HEPAK at these states and the ADR's network form. Mass swing rho(T, P_h) - rho(T, P_l) times
# volume, in kg: pores 9.058225e-4 (60 K), 8.335998e-4 (20 K), 1.131456e-3 (8 K); pulse tubes by adaptive quadrature
# over a linear profile from 300 K, 2.581565e-4 (to 40 K) and 1.457110e-3 (to 4 K). Boundary flows are f times the
# swing beyond each boundary.
_MASS_FLOW = np.array([5.503374e-3, 4.416387e-3, 4.106599e-3, 3.106279e-3, 1.748532e-3])
# Boundary conduction 2 a (m c_p)^2 (T_w - T_c) / (f (C_w + C_c)), c_p at 1.4 MPa and the boundary's mean temperature:
# 326.9320, 18.23414, 553.1576, 202.8832, 30.42731 W. Real-gas flow (1 - beta T) H_2 at 1.4 MPa toward stage 2 at
# 40 K, 14 K, and 4 K: 0.01004304, -0.5241515, 1.256511 W. Stage sinks H_1 and H_2 land on the warm flange.
_HEAT_FLOW = [-295.4320, 308.6979, -564.9335, 350.8086, 170.6753, 30.18382]


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


def test_heat_flow_matches_hand_calculation():
    heat_flow = calc_heat_flow(
        _TEMPERATURE,
        heat_capacity=_HEAT_CAPACITY,
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
