"""Tests for the core component types."""

import dataclasses

import numpy as np
import pytest

from pulse_tube_cooldown.core.errors import ConfigurationError
from pulse_tube_cooldown.core.types import Material, Node


def _linear_material() -> Material:
    return Material(
        name="Toy",
        density=1.0,
        specific_heat=lambda temperature: 2.0 * temperature,
        specific_heat_integral=lambda temperature: temperature**2,
        thermal_conductivity=lambda temperature: 3.0 + 0.0 * temperature,
        conductivity_integral=lambda temperature: 3.0 * temperature,
    )


def test_node_without_materials_is_rejected():
    with pytest.raises(ConfigurationError, match="PRP"):
        Node("PRP", ())


@pytest.mark.parametrize(("temperature", "capacity"), [(4.0, 24.0), (300.0, 1800.0)])
def test_node_evaluates_heat_capacity_at_its_temperature(temperature, capacity):
    node = Node("PRP", ((_linear_material(), 3.0),))

    # Matrix mass 3 kg times specific heat 2 T J/(kg K).
    assert node.calc_heat_capacity(temperature) == pytest.approx(capacity)


def test_node_heat_capacity_sums_its_material_contributions():
    linear = _linear_material()
    quadratic = dataclasses.replace(linear, name="Quadratic", specific_heat=lambda temperature: 0.5 * temperature**2)
    constant = dataclasses.replace(
        linear, name="Constant", specific_heat=lambda temperature: np.full_like(temperature, 400.0)
    )
    node = Node("RGP", ((linear, 6.35), (quadratic, 5.0), (constant, 0.25)))

    # At 10 K: 6.35 kg * 20 J/(kg K) + 5 kg * 50 J/(kg K) + 0.25 kg * 400 J/(kg K).
    assert node.calc_heat_capacity(10.0) == pytest.approx(477.0)


@pytest.mark.parametrize("mass", [0.0, -1.0])
def test_node_with_nonpositive_mass_is_rejected(mass):
    steel = Material(
        name="Steel",
        density=1.0,
        specific_heat=lambda t: t,
        specific_heat_integral=lambda t: t**2 / 2.0,
        thermal_conductivity=lambda t: t,
        conductivity_integral=lambda t: t**2 / 2.0,
    )

    with pytest.raises(ConfigurationError, match=r"RGP(.|\n)*Steel"):
        Node("RGP", ((_linear_material(), 6.35), (steel, mass)))
