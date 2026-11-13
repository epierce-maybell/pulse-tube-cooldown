"""Building models from TOML configuration files."""

import os
import tomllib
from collections.abc import Callable, Mapping
from typing import Any

import numpy as np

from pulse_tube_cooldown.core.errors import ConfigurationError
from pulse_tube_cooldown.core.types import Graph, Link, Model, Node, Phase, Procedure
from pulse_tube_cooldown.library.conditions import load_condition
from pulse_tube_cooldown.library.links.conduction import load_conduction
from pulse_tube_cooldown.library.links.exchange_gas import load_exchange_gas
from pulse_tube_cooldown.library.links.gas_gap_switch import load_gas_gap_switch
from pulse_tube_cooldown.library.links.gas_pushing import load_gas_pushing
from pulse_tube_cooldown.library.links.joint_contact import load_joint_contact
from pulse_tube_cooldown.library.links.multi_layer_insulation import load_multi_layer_insulation
from pulse_tube_cooldown.library.links.pulse_tube.empirical import load_pulse_tube_head
from pulse_tube_cooldown.library.links.pulse_tube.energy_flow import load_energy_flow_pulse_tube_head
from pulse_tube_cooldown.library.links.radiation import load_radiation
from pulse_tube_cooldown.library.materials.load import load_material

__all__ = ["LINK_TYPES", "configure_model", "load_fridge", "load_link", "load_procedure"]

LINK_TYPES: Mapping[str, Callable[[Mapping[str, Any]], Link]] = {
    "conduction": load_conduction,
    "exchange_gas": load_exchange_gas,
    "gas_gap_switch": load_gas_gap_switch,
    "gas_pushing": load_gas_pushing,
    "joint_contact": load_joint_contact,
    "multi_layer_insulation": load_multi_layer_insulation,
    "radiation": load_radiation,
    "pulse_tube_head": lambda config: load_pulse_tube_head(config, load_graph=load_fridge),
    "energy_flow_pulse_tube_head": load_energy_flow_pulse_tube_head,
}


def configure_model(path_to_config: str | os.PathLike[str]) -> Model:
    """Builds a model from a TOML configuration file.

    Args:
        path_to_config: Path to the TOML configuration file.

    Returns:
        The model, ready to simulate.

    Raises:
        ConfigurationError: A table is invalid, `initial_temperature` does not match the nodes, or a fixed boundary
            node is unknown.
    """
    with open(path_to_config, "rb") as file:
        config = tomllib.load(file)

    graph = load_fridge(config["fridge"])
    procedure = load_procedure(config["procedure"], graph)
    solver_options = dict(config["simulation"])
    names = [node.name for node in graph.nodes]

    initial_temperature = np.asarray(solver_options.pop("initial_temperature"), dtype=np.float64)
    if initial_temperature.ndim == 0:
        initial_temperature = np.full(len(names), initial_temperature)
    elif initial_temperature.shape != (len(names),):
        raise ConfigurationError(
            f"`initial_temperature` has {initial_temperature.size} values but needs one for each node: {names}."
        )

    fixed_boundary_nodes = tuple(graph.get_node_index(name) for name in solver_options.pop("fixed_boundary_nodes", []))

    return Model(
        config.get("meta", {}),
        graph.nodes,
        procedure,
        initial_temperature,
        fixed_boundary_nodes,
        solver_options,
    )


def load_fridge(fridge: Mapping[str, Any]) -> Graph:
    """Builds the full graph, with every link active, from the `[fridge]` table.

    Args:
        fridge: An array of `nodes` tables, each with a `name` and a `composition` of [material name, mass in kg]
            pairs, and an optional array of `links` tables, each accepted by `load_link`.

    Returns:
        The graph of the configured nodes followed by any `internal_nodes` supplied by links, and the links, in
        configuration order.

    Raises:
        ConfigurationError: A node or link is missing a key, or a name refers to an unknown material, node, or link
            type.
    """
    links = tuple(load_link(link) for link in fridge.get("links", []))
    configured = (_load_node(node) for node in fridge["nodes"])
    nodes = (*configured, *(node for link in links for node in getattr(link, "internal_nodes", ())))

    return Graph(nodes=nodes, links=links)


def load_link(config: Mapping[str, Any]) -> Link:
    """Builds a link from its configuration table with the loader registered for its type.

    Args:
        config: The link's `type`, a key of `LINK_TYPES`, and the fields that type's loader reads:

            name = "Rod"
            type = "conduction"
            nodes = ["CFP", "STP"]
            material = "304 SS"
            area = 1.0e-4
            length = 0.1

    Returns:
        The link.

    Raises:
        ConfigurationError: The type is unknown, a key is missing, or the type's loader rejects the table.
    """
    link_type = config.get("type")

    if link_type not in LINK_TYPES:
        raise ConfigurationError(
            f"Link {config.get('name')!r} has unknown type {link_type!r}; known types are {', '.join(LINK_TYPES)}."
        )

    try:
        return LINK_TYPES[link_type](config)
    except KeyError as error:
        raise ConfigurationError(f"Link {config.get('name')!r} is missing {error}.") from error


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
