"""Tests for heat flow through a multi-layer insulation blanket."""

import numpy as np

from pulse_tube_cooldown.library.links.multi_layer_insulation import MultiLayerInsulation


def _blanket(area=0.5, layers=10, layer_density=2000.0) -> MultiLayerInsulation:
    return MultiLayerInsulation(
        name="Blanket",
        nodes=("PRP", "Can"),
        area=area,
        layers=layers,
        layer_density=layer_density,
        spacer_coefficient=8.95e-8 * 100.0**-2.56,
        layer_density_exponent=2.56,
        radiation_coefficient=5.39e-10,
        emissivity=0.03,
        source="Lockheed 1974",
    )


def test_matches_hand_calculation():
    # Spacers in published per-cm units, 20 layers/cm: 8.95e-8 * 20^2.56 * 175 * (300 - 50) / 11 = 0.762158 W/m^2.
    # Radiation: 5.39e-10 * 0.03 * (300^4.67 - 50^4.67) / 10 = 0.598088 W/m^2.
    # Over 0.5 m^2: 0.680123 W.
    flow = _blanket().calc_energy_balance_terms(np.array([50.0, 300.0])).heat_flow

    np.testing.assert_allclose(flow, [0.680123, -0.680123], rtol=1e-5)
