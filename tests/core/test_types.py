"""Tests for the core component types."""

import dataclasses

import numpy as np
import pytest
from numpy.typing import NDArray

from pulse_tube_cooldown.core.errors import ConfigurationError
from pulse_tube_cooldown.core.types import Gas, Link, Material, Node


def _linear_material() -> Material:
    return Material(
        name="Toy",
        specific_heat=lambda temperature: 2.0 * temperature,
        thermal_conductivity=lambda temperature: 3.0 + 0.0 * temperature,
    )


@dataclasses.dataclass(frozen=True, slots=True)
class _ConstantConductance:
    """Two-node link with constant conductance, following the Link convention."""

    name: str
    nodes: tuple[str, ...]
    conductance: float

    def calc_heat_flow(self, temperature: NDArray[np.float64]) -> NDArray[np.float64]:
        into_first = self.conductance * (temperature[1] - temperature[0])
        return np.array([into_first, -into_first])


@dataclasses.dataclass(frozen=True, slots=True)
class _ConstantSink:
    """Single-node link that removes a fixed power."""

    name: str
    nodes: tuple[str, ...]
    power: float

    def calc_heat_flow(self, temperature: NDArray[np.float64]) -> NDArray[np.float64]:
        return np.full(temperature.shape, -self.power)


def test_material_evaluates_properties_over_temperature_arrays():
    material = _linear_material()
    temperature = np.array([4.0, 77.0, 300.0])

    np.testing.assert_allclose(material.specific_heat(temperature), [8.0, 154.0, 600.0])
    np.testing.assert_allclose(material.thermal_conductivity(temperature), [3.0, 3.0, 3.0])


def test_material_and_gas_are_keyword_only():
    with pytest.raises(TypeError):
        Material("Toy", lambda t: t, lambda t: t)  # type: ignore[misc]
    with pytest.raises(TypeError):
        Gas("He-4", 4.0026e-3, 5.0 / 3.0)  # type: ignore[misc]


def test_gas_stores_si_properties():
    helium = Gas(name="He-4", molar_mass=4.0026e-3, heat_capacity_ratio=5.0 / 3.0)

    assert helium.name == "He-4"
    assert helium.molar_mass == pytest.approx(4.0026e-3)
    assert helium.heat_capacity_ratio == pytest.approx(5.0 / 3.0)


def test_node_composition_pairs_materials_with_masses():
    copper = _linear_material()
    steel = Material(name="Steel", specific_heat=lambda t: t, thermal_conductivity=lambda t: t)
    node = Node("RGP", ((copper, 6.35), (steel, 5.0)))

    assert node.name == "RGP"
    assert [(material.name, mass) for material, mass in node.composition] == [("Toy", 6.35), ("Steel", 5.0)]


def test_node_without_materials_is_rejected():
    with pytest.raises(ConfigurationError, match="PRP"):
        Node("PRP", ())


@pytest.mark.parametrize("mass", [0.0, -1.0])
def test_node_with_nonpositive_mass_is_rejected(mass):
    steel = Material(name="Steel", specific_heat=lambda t: t, thermal_conductivity=lambda t: t)

    with pytest.raises(ConfigurationError, match=r"RGP(.|\n)*Steel"):
        Node("RGP", ((_linear_material(), 6.35), (steel, mass)))


@pytest.mark.parametrize(
    "component",
    [
        _linear_material(),
        Gas(name="He-4", molar_mass=4.0026e-3, heat_capacity_ratio=5.0 / 3.0),
        Node("PRP", ((_linear_material(), 1.0),)),
    ],
)
def test_components_are_immutable(component):
    with pytest.raises(dataclasses.FrozenInstanceError):
        component.name = "changed"


def test_two_node_link_returns_signed_flows_in_declared_node_order():
    link: Link = _ConstantConductance(name="strut", nodes=("PRP", "RGP"), conductance=0.5)

    flow = link.calc_heat_flow(np.array([300.0, 100.0]))

    np.testing.assert_allclose(flow, [-100.0, 100.0])
    assert flow.sum() == pytest.approx(0.0)


def test_reversing_declared_node_order_reverses_flow_order():
    forward = _ConstantConductance(name="strut", nodes=("PRP", "RGP"), conductance=0.5)
    reverse = _ConstantConductance(name="strut", nodes=("RGP", "PRP"), conductance=0.5)

    np.testing.assert_allclose(
        reverse.calc_heat_flow(np.array([100.0, 300.0])),
        forward.calc_heat_flow(np.array([300.0, 100.0]))[::-1],
    )


def test_single_node_sink_returns_negative_flow():
    sink = _ConstantSink(name="cooler", nodes=("PRP",), power=40.0)

    np.testing.assert_allclose(sink.calc_heat_flow(np.array([50.0])), [-40.0])
