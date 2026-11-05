"""Tests for heat flow across a pressed joint."""

import numpy as np

from pulse_tube_cooldown.library.links.joint_contact import JointContact


def _joint(force=1.0e4, conductance_per_force=1.9e-3, reference_temperature=4.2, temperature_exponent=1.1):
    return JointContact(
        name="Stage 2 flange",
        nodes=("PT stage 2", "CFP"),
        force=force,
        conductance_per_force=conductance_per_force,
        reference_temperature=reference_temperature,
        temperature_exponent=temperature_exponent,
        source="Dhuley 2019",
    )


def test_matches_hand_calculation_at_reference_temperature():
    # Mean 4.2 K: G = 1.9e-3 W/(K N) * 1e4 N = 19 W/K; 19 W/K * 0.4 K = 7.6 W into the colder node.
    flow = _joint().calc_energy_balance_terms(np.array([4.0, 4.4])).heat_flow

    np.testing.assert_allclose(flow, [7.6, -7.6], rtol=1e-12)


def test_conductance_follows_power_law_at_mean_temperature():
    # Mean 20 K: G = 1e-3 W/(K N) * 1e3 N * (20 / 10)^2 = 4 W/K; 4 W/K * 10 K = 40 W into the colder node.
    joint = _joint(force=1.0e3, conductance_per_force=1.0e-3, reference_temperature=10.0, temperature_exponent=2.0)

    flow = joint.calc_energy_balance_terms(np.array([15.0, 25.0])).heat_flow

    np.testing.assert_allclose(flow, [40.0, -40.0], rtol=1e-12)


def test_warmer_first_node_loses_heat():
    joint = _joint()

    forward = joint.calc_energy_balance_terms(np.array([40.0, 60.0])).heat_flow
    reverse = joint.calc_energy_balance_terms(np.array([60.0, 40.0])).heat_flow

    assert reverse[0] < 0.0 < reverse[1]
    np.testing.assert_allclose(reverse, -forward, rtol=1e-12)
