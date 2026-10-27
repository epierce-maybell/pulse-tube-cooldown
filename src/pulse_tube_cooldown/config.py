"""Building models from TOML configuration files."""

from collections.abc import Mapping
from typing import Any

from pulse_tube_cooldown.core.errors import ConfigurationError
from pulse_tube_cooldown.core.types import Graph, Node, Phase, Procedure
from pulse_tube_cooldown.library.conditions import load_condition
from pulse_tube_cooldown.library.links.load import load_link
from pulse_tube_cooldown.library.materials.load import load_material

__all__ = ["load_fridge", "load_procedure"]


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


def load_procedure(procedure: Mapping[str, Any], graph: Graph) -> Procedure:
    """Builds the procedure from the `[procedure]` table, activating the full graph's links phase by phase.

    Args:
        procedure: An array of `phases` tables and an optional `until` condition that stops the run. Each phase has a
            `name`, optional `disable` and `enable` arrays of link names, applied to the previous phase's links, and
            an `until` condition unless it is the last phase. Conditions are tables accepted by `load_condition`.
        graph: The full graph, with every link active, whose nodes the conditions refer to.

    Returns:
        The procedure, with each phase's links in graph order.

    Raises:
        ConfigurationError: A phase disables a link that is not active or enables one that is not inactive, a
            condition is invalid, or the phases do not form a valid procedure.
    """
    links = {link.name for link in graph.links}
    active = set(links)
    phases = []

    for config in procedure["phases"]:
        disable, enable = set(config.get("disable", [])), set(config.get("enable", []))
        invalid = (disable - active) | (enable - (links - active))
        if invalid:
            raise ConfigurationError(
                f"Phase {config.get('name')!r} can only disable active links and enable inactive ones, not "
                f"{sorted(invalid)}; active links are {sorted(active)}."
            )

        active = (active - disable) | enable
        until = load_condition(config["until"], graph) if "until" in config else None
        phases.append(Phase(config["name"], tuple(link for link in graph.links if link.name in active), until))

    stop_condition = load_condition(procedure["until"], graph) if "until" in procedure else None

    return Procedure(tuple(phases), stop_condition)


def _load_node(config: Mapping[str, Any]) -> Node:
    """Builds a node from its table, loading each material in its composition."""
    try:
        composition = tuple((load_material(material), mass) for material, mass in config["composition"])
        return Node(config["name"], composition)
    except KeyError as error:
        raise ConfigurationError(f"Node {config.get('name')!r} is missing {error}.") from error
