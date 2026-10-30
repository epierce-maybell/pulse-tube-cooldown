"""Tests for the engine's thermal equations."""

import dataclasses

import numpy as np
from numpy.typing import NDArray

from pulse_tube_cooldown.core.types import Graph, Material, Node
from pulse_tube_cooldown.engine.equations import (
    calc_heat_capacity,
    calc_link_heat_flows,
    calc_net_heat_flow,
    calc_temperature_rate_of_change,
)


def _material(name: str, specific_heat) -> Material:
    return Material(
        name=name,
        density=1.0,
        specific_heat=specific_heat,
        thermal_conductivity=lambda t: 0.0 * t,
        conductivity_integral=lambda t: 0.0 * t,
    )


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


@dataclasses.dataclass(frozen=True, slots=True)
class _Conductance:
    """Constant conductance G from the first node to the second."""

    name: str
    nodes: tuple[str, str]
    conductance: float

    def calc_heat_flow(self, temperature: NDArray[np.float64]) -> NDArray[np.float64]:
        flow = self.conductance * (temperature[1] - temperature[0])
        return np.array([flow, -flow])


@dataclasses.dataclass(frozen=True, slots=True)
class _SharedBath:
    """Relaxes each of its nodes toward their mean temperature with conductance G."""

    name: str
    nodes: tuple[str, ...]
    conductance: float

    def calc_heat_flow(self, temperature: NDArray[np.float64]) -> NDArray[np.float64]:
        return self.conductance * (temperature.mean() - temperature)


@dataclasses.dataclass(frozen=True, slots=True)
class _Source:
    """Constant power into a single node; negative power is a sink."""

    name: str
    nodes: tuple[str]
    power: float

    def calc_heat_flow(self, temperature: NDArray[np.float64]) -> NDArray[np.float64]:
        return np.full_like(temperature, self.power)


def _network() -> Graph:
    nodes = tuple(Node(name, ((CONSTANT, 1.0),)) for name in ("A", "B", "C", "D"))
    links = (
        _Conductance("AB", ("A", "B"), 2.0),
        _Conductance("BA", ("B", "A"), 0.5),
        _SharedBath("gas", ("D", "B", "C"), 1.0),
        _Source("heater", ("C",), 5.0),
        _Source("cooler", ("D",), -3.0),
    )
    return Graph(nodes=nodes, links=links)


# Node temperatures A, B, C, D in K.
_TEMPERATURE = np.array([300.0, 100.0, 40.0, 4.0])


def test_link_heat_flows_place_each_link_output_at_its_nodes():
    flows = calc_link_heat_flows(_network(), _TEMPERATURE)

    # Columns A, B, C, D; NaN where a link does not touch a node.
    # AB: 2 * (100 - 300) into A. BA: 0.5 * (300 - 100) into B.
    # gas: mean of D, B, C is (4 + 100 + 40) / 3 = 48, so D +44, B -52, C +8.
    nan = np.nan
    expected = [
        [-400.0, 400.0, nan, nan],
        [-100.0, 100.0, nan, nan],
        [nan, -52.0, 8.0, 44.0],
        [nan, nan, 5.0, nan],
        [nan, nan, nan, -3.0],
    ]
    np.testing.assert_allclose(flows, expected)


def test_net_heat_flow_sums_signed_flows_into_each_node():
    net = calc_net_heat_flow(_network(), _TEMPERATURE)

    # A: -400 - 100. B: 400 + 100 - 52. C: 8 + 5. D: 44 - 3.
    np.testing.assert_allclose(net, [-500.0, 448.0, 13.0, 41.0])


def test_net_heat_flow_propagates_nan_from_a_link_only_to_its_nodes():
    nodes = tuple(Node(name, ((CONSTANT, 1.0),)) for name in ("A", "B", "C"))
    graph = Graph(nodes=nodes, links=(_Conductance("AB", ("A", "B"), 1.0), _Source("broken", ("C",), np.nan)))

    # A: 1 * (20 - 10) = 10. B: -10. C: the broken link's NaN.
    np.testing.assert_allclose(calc_net_heat_flow(graph, np.array([10.0, 20.0, 30.0])), [10.0, -10.0, np.nan])


def test_graph_without_links_has_no_heat_flow():
    graph = Graph(nodes=(Node("A", ((CONSTANT, 1.0),)), Node("B", ((CONSTANT, 1.0),))), links=())

    assert calc_link_heat_flows(graph, np.array([10.0, 20.0])).shape == (0, 2)
    np.testing.assert_allclose(calc_net_heat_flow(graph, np.array([10.0, 20.0])), [0.0, 0.0])


def _boundary_network() -> Graph:
    nodes = (
        Node("PT", ((LINEAR, 3.0),)),
        Node("RGP", ((QUADRATIC, 2.0),)),
        Node("CFP", ((CONSTANT, 0.5),)),
    )
    links = (
        _Conductance("PT-RGP", ("PT", "RGP"), 2.0),
        _Conductance("RGP-CFP", ("RGP", "CFP"), 1.0),
        _Source("heater", ("CFP",), 5.0),
    )
    return Graph(nodes=nodes, links=links)


# Node temperatures PT, RGP, CFP in K.
# Net heat flow: PT 2 * (10 - 300) = -580. RGP 580 + 1 * (4 - 10) = 574. CFP 6 + 5 = 11.
# Heat capacity: PT 3 * 2 * 300 = 1800. RGP 2 * 0.5 * 10^2 = 100. CFP 0.5 * 400 = 200.
_BOUNDARY_TEMPERATURE = np.array([300.0, 10.0, 4.0])


def test_temperature_rate_of_change_divides_net_heat_flow_by_capacity():
    rate = calc_temperature_rate_of_change(_boundary_network(), (), 0.0, _BOUNDARY_TEMPERATURE)

    np.testing.assert_allclose(rate, [-580.0 / 1800.0, 5.74, 0.055])


def test_fixed_boundary_node_holds_while_its_neighbors_exchange_heat_with_it():
    rate = calc_temperature_rate_of_change(_boundary_network(), (0,), 0.0, _BOUNDARY_TEMPERATURE)

    # RGP still receives 580 W from PT.
    np.testing.assert_allclose(rate, [0.0, 5.74, 0.055])


def test_every_fixed_boundary_node_holds():
    rate = calc_temperature_rate_of_change(_boundary_network(), (2, 0), 0.0, _BOUNDARY_TEMPERATURE)

    np.testing.assert_allclose(rate, [0.0, 5.74, 0.0])
