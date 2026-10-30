"""Tests for graph connectivity and name lookups."""

import dataclasses

import numpy as np
import pytest
from numpy.typing import NDArray

from pulse_tube_cooldown.core.errors import ConfigurationError, PTCError
from pulse_tube_cooldown.core.types import Graph, Material, Node

_TOY = Material(
    name="Toy",
    density=1.0,
    specific_heat=lambda t: t,
    thermal_conductivity=lambda t: t,
    conductivity_integral=lambda t: t**2 / 2.0,
)


@dataclasses.dataclass(frozen=True, slots=True)
class _StubLink:
    """Link whose flow is irrelevant to connectivity."""

    name: str
    nodes: tuple[str, ...]

    def calc_heat_flow(self, temperature: NDArray[np.float64]) -> NDArray[np.float64]:
        return np.zeros_like(temperature)


def _nodes(*names: str) -> tuple[Node, ...]:
    return tuple(Node(name, ((_TOY, 1.0),)) for name in names)


def test_connectivity_maps_link_local_order_to_graph_indices():
    graph = Graph(
        nodes=_nodes("PRP", "RGP", "STILL", "CP"),
        links=(
            _StubLink("forward", ("PRP", "RGP")),
            _StubLink("reverse", ("RGP", "PRP")),
            _StubLink("parallel", ("PRP", "RGP")),
            _StubLink("shared gas", ("CP", "PRP", "STILL")),
            _StubLink("sink", ("STILL",)),
        ),
    )

    assert graph.connectivity == ((0, 1), (1, 0), (0, 1), (3, 0, 2), (2,))


def test_lookups_return_graph_indices():
    graph = Graph(
        nodes=_nodes("PRP", "RGP", "STILL"),
        links=(_StubLink("a", ("RGP", "STILL")), _StubLink("b", ("PRP", "RGP"))),
    )

    assert [graph.get_node_index(name) for name in ("PRP", "RGP", "STILL")] == [0, 1, 2]
    assert [graph.get_link_index(name) for name in ("a", "b")] == [0, 1]


def test_graph_without_links_has_empty_connectivity():
    graph = Graph(nodes=_nodes("PRP"), links=())

    assert graph.connectivity == ()


@pytest.mark.parametrize(
    ("lookup", "name"),
    [("get_node_index", "Bath"), ("get_link_index", "Rod C")],
)
def test_unknown_lookup_names_raise_configuration_errors_naming_them(lookup, name):
    graph = Graph(nodes=_nodes("Plate", "Sink"), links=(_StubLink("Rod A", ("Plate", "Sink")),))

    with pytest.raises(ConfigurationError, match=repr(name)) as error:
        getattr(graph, lookup)(name)

    assert isinstance(error.value, ValueError)
    assert isinstance(error.value, PTCError)


@pytest.mark.parametrize(
    ("nodes", "links", "mentions"),
    [
        (_nodes("Plate", "Bath"), (_StubLink("Rod A", ("Plate", "Bth")),), ["Rod A", "Bth"]),
        (_nodes("Plate", "Plate"), (), ["Plate"]),
        (_nodes("Plate", "Bath"), (_StubLink("Rod", ("Plate", "Bath")), _StubLink("Rod", ("Bath",))), ["Rod"]),
    ],
    ids=["unknown node reference", "duplicate node", "duplicate link"],
)
def test_invalid_structure_raises_value_error_naming_offenders(nodes, links, mentions):
    with pytest.raises(ConfigurationError) as error:
        Graph(nodes=nodes, links=links)

    assert isinstance(error.value, ValueError)
    assert isinstance(error.value, PTCError)
    for mention in mentions:
        assert repr(mention) in str(error.value)


def test_graph_is_immutable():
    graph = Graph(nodes=_nodes("Plate"), links=())

    with pytest.raises(dataclasses.FrozenInstanceError):
        graph.nodes = ()  # type: ignore[misc]
