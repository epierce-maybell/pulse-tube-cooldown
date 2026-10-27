"""Building models from TOML configuration files."""

from collections.abc import Mapping
from typing import Any

from pulse_tube_cooldown.core.errors import ConfigurationError
from pulse_tube_cooldown.core.types import Graph, Node
from pulse_tube_cooldown.library.links.load import load_link
from pulse_tube_cooldown.library.materials.load import load_material

__all__ = ["load_fridge"]


def load_fridge(fridge: Mapping[str, Any]) -> Graph:
    """Builds the full graph, with every link active, from the `[fridge]` table.

    Args:
        fridge: An array of `nodes` tables, each with a `name` and a `composition` of [material name, mass in kg]
            pairs, and an optional array of `links` tables, each accepted by `load_link`.

    Returns:
        The graph of the nodes and links, in configuration order.

    Raises:
        ConfigurationError: A node or link is missing a key, or a name refers to an unknown material, node, or link
            type.
    """
    nodes = tuple(_load_node(node) for node in fridge["nodes"])
    links = tuple(load_link(link) for link in fridge.get("links", []))

    return Graph(nodes=nodes, links=links)


def _load_node(config: Mapping[str, Any]) -> Node:
    """Builds a node from its table, loading each material in its composition."""
    try:
        composition = tuple((load_material(material), mass) for material, mass in config["composition"])
        return Node(config["name"], composition)
    except KeyError as error:
        raise ConfigurationError(f"Node {config.get('name')!r} is missing {error}.") from error
