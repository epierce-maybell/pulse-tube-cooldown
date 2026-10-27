"""Tests for loading configuration files."""

import tomllib

import pytest

from pulse_tube_cooldown.config import load_fridge
from pulse_tube_cooldown.core.errors import ConfigurationError
from pulse_tube_cooldown.core.types import Graph
from pulse_tube_cooldown.library.links.conduction import Conduction

_FRIDGE = tomllib.loads(
    """
    [[nodes]]
    name = "Plate"
    composition = [["OFHC Cu RRR100", 1.5]]

    [[nodes]]
    name = "Mixed"
    composition = [["OFHC Cu RRR100", 2.0], ["304 SS", 0.5]]

    [[nodes]]
    name = "Bath"
    composition = [["304 SS", 3.0]]

    [[links]]
    name = "Rod A"
    type = "conduction"
    nodes = ["Plate", "Bath"]
    material = "304 SS"
    area = 1.0e-4
    length = 0.1

    [[links]]
    name = "Rod B"
    type = "conduction"
    nodes = ["Bath", "Mixed"]
    material = "OFHC Cu RRR100"
    area = 2.0e-5
    length = 0.3
    """
)


def _fridge(**changes):
    return _FRIDGE | changes


def test_loads_nodes_in_order_with_compositions():
    graph = load_fridge(_FRIDGE)

    assert isinstance(graph, Graph)
    assert [node.name for node in graph.nodes] == ["Plate", "Mixed", "Bath"]
    assert [[(material.name, mass) for material, mass in node.composition] for node in graph.nodes] == [
        [("OFHC Cu RRR100", 1.5)],
        [("OFHC Cu RRR100", 2.0), ("304 SS", 0.5)],
        [("304 SS", 3.0)],
    ]


def test_loads_links_with_parameters_and_connectivity():
    graph = load_fridge(_FRIDGE)

    assert all(isinstance(link, Conduction) for link in graph.links)
    assert [(link.name, link.nodes, link.material.name, link.area, link.length) for link in graph.links] == [
        ("Rod A", ("Plate", "Bath"), "304 SS", 1.0e-4, 0.1),
        ("Rod B", ("Bath", "Mixed"), "OFHC Cu RRR100", 2.0e-5, 0.3),
    ]
    assert graph.connectivity == ((0, 2), (2, 1))


def test_fridge_without_links_has_only_nodes():
    graph = load_fridge({"nodes": _FRIDGE["nodes"]})

    assert len(graph.nodes) == 3
    assert graph.links == ()


def test_link_to_unknown_node_names_link_and_node():
    links = [_FRIDGE["links"][0] | {"nodes": ["Plate", "Tub"]}]

    with pytest.raises(ConfigurationError, match=r"(?s)Rod A.*Tub"):
        load_fridge(_fridge(links=links))


def test_unknown_link_type_names_the_link():
    links = [_FRIDGE["links"][0] | {"type": "strut"}]

    with pytest.raises(ConfigurationError, match=r"(?s)Rod A.*strut"):
        load_fridge(_fridge(links=links))


def test_unknown_node_material_names_the_material():
    nodes = [{"name": "Plate", "composition": [["Unobtainium", 1.0]]}]

    with pytest.raises(ConfigurationError, match="Unobtainium"):
        load_fridge({"nodes": nodes})


def test_duplicate_node_names_are_rejected():
    nodes = _FRIDGE["nodes"] + [_FRIDGE["nodes"][0]]

    with pytest.raises(ConfigurationError, match="Plate"):
        load_fridge(_fridge(nodes=nodes))


@pytest.mark.parametrize("node", [{"composition": [["304 SS", 1.0]]}, {"name": "Plate"}])
def test_node_missing_a_key_is_rejected(node):
    with pytest.raises(ConfigurationError, match="name|composition"):
        load_fridge({"nodes": [node]})
