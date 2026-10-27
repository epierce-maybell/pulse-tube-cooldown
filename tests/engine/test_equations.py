"""Tests for the engine's thermal equations."""

import numpy as np

from pulse_tube_cooldown.core.types import Material, Node
from pulse_tube_cooldown.engine.equations import calc_heat_capacity


def _material(name: str, specific_heat) -> Material:
    return Material(name=name, specific_heat=specific_heat, thermal_conductivity=lambda t: 0.0 * t)


LINEAR = _material("Linear", lambda t: 2.0 * t)
QUADRATIC = _material("Quadratic", lambda t: 0.5 * t**2)
CONSTANT = _material("Constant", lambda t: 400.0 + 0.0 * t)


def test_single_material_node_capacity_follows_temperature():
    node = Node("PRP", ((LINEAR, 3.0),))

    # C = m * c(T) = 3 kg * 2 T.
    np.testing.assert_allclose(calc_heat_capacity([node], np.array([4.0])), [24.0])
    np.testing.assert_allclose(calc_heat_capacity([node], np.array([300.0])), [1800.0])


def test_mixed_material_node_sums_mass_weighted_specific_heats():
    node = Node("RGP", ((LINEAR, 6.35), (QUADRATIC, 5.0), (CONSTANT, 0.25)))

    # At 10 K: 6.35 * 20 + 5.0 * 50 + 0.25 * 400 = 127 + 250 + 100.
    np.testing.assert_allclose(calc_heat_capacity([node], np.array([10.0])), [477.0])


def test_each_node_is_evaluated_at_its_own_temperature():
    nodes = [
        Node("PRP", ((LINEAR, 1.0),)),
        Node("RGP", ((QUADRATIC, 2.0), (LINEAR, 1.0))),
        Node("CFP", ((CONSTANT, 0.5),)),
    ]

    # PRP at 300 K: 1 * 600. RGP at 4 K: 2 * 8 + 1 * 8. CFP at 50 K: 0.5 * 400.
    np.testing.assert_allclose(calc_heat_capacity(nodes, np.array([300.0, 4.0, 50.0])), [600.0, 24.0, 200.0])
