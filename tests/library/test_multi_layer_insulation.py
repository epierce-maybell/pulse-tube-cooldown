"""Tests for heat flow through a multi-layer insulation blanket."""

import numpy as np
import pytest

from pulse_tube_cooldown.library.links.multi_layer_insulation import MultiLayerInsulation


def _blanket(area=0.5, layers=10, layer_density=20.0) -> MultiLayerInsulation:
    return MultiLayerInsulation(
        name="Blanket",
        nodes=("PRP", "Can"),
        area=area,
        layers=layers,
        layer_density=layer_density,
        spacer_coefficient=8.95e-8,
        layer_density_exponent=2.56,
        radiation_coefficient=5.39e-10,
        emissivity=0.03,
        source="Lockheed 1974",
    )


def test_matches_hand_calculation():
    # Spacers: 8.95e-8 * 20^2.56 * 175 * (300 - 50) / 11 = 0.762158 W/m^2.
    # Radiation: 5.39e-10 * 0.03 * (300^4.67 - 50^4.67) / 10 = 0.598088 W/m^2.
    # Over 0.5 m^2: 0.680123 W.
    flow = _blanket().calc_heat_flow(np.array([50.0, 300.0]))

    np.testing.assert_allclose(flow, [0.680123, -0.680123], rtol=1e-5)


def test_equal_temperatures_carry_no_heat():
    flow = _blanket().calc_heat_flow(np.array([77.0, 77.0]))

    np.testing.assert_array_equal(flow, [0.0, 0.0])


def test_warmer_inner_node_loses_heat_to_outer_node():
    blanket = _blanket()

    forward = blanket.calc_heat_flow(np.array([40.0, 250.0]))
    reverse = blanket.calc_heat_flow(np.array([250.0, 40.0]))

    assert reverse[0] < 0.0 < reverse[1]
    assert reverse.sum() == pytest.approx(0.0, abs=1e-12 * abs(reverse[0]))
    np.testing.assert_allclose(reverse, -forward, rtol=1e-12)


def test_heat_scales_with_area():
    temperature = np.array([4.0, 300.0])

    small = _blanket(area=0.1).calc_heat_flow(temperature)
    large = _blanket(area=0.3).calc_heat_flow(temperature)

    np.testing.assert_allclose(large, 3.0 * small, rtol=1e-12)
