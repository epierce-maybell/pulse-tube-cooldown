"""Links that transfer heat between nodes, built from configuration.

A link's configuration is a table holding its `type`, a key of `LINK_TYPES`, and the fields that type's loader reads:

    name = "Rod"
    type = "conduction"
    nodes = ["CFP", "STP"]
    material = "304 SS"
    area = 1.0e-4
    length = 0.1
"""

from collections.abc import Callable, Mapping
from typing import Any

from pulse_tube_cooldown.core.errors import ConfigurationError
from pulse_tube_cooldown.core.types import Link, Node
from pulse_tube_cooldown.library.links.conduction import load_conduction
from pulse_tube_cooldown.library.links.cooler.link import load_cooler, load_matrix_nodes
from pulse_tube_cooldown.library.links.exchange_gas import load_exchange_gas
from pulse_tube_cooldown.library.links.gas_gap_switch import load_gas_gap_switch
from pulse_tube_cooldown.library.links.gas_pushing import load_gas_pushing
from pulse_tube_cooldown.library.links.joint_contact import load_joint_contact
from pulse_tube_cooldown.library.links.multi_layer_insulation import load_multi_layer_insulation
from pulse_tube_cooldown.library.links.radiation import load_radiation

__all__ = ["LINK_TYPES", "load_link", "load_link_nodes"]

LINK_TYPES: Mapping[str, Callable[[Mapping[str, Any]], Link]] = {
    "conduction": load_conduction,
    "cooler": load_cooler,
    "exchange_gas": load_exchange_gas,
    "gas_gap_switch": load_gas_gap_switch,
    "gas_pushing": load_gas_pushing,
    "joint_contact": load_joint_contact,
    "multi_layer_insulation": load_multi_layer_insulation,
    "radiation": load_radiation,
}


def load_link(config: Mapping[str, Any]) -> Link:
    """Builds a link from its configuration table with the loader registered for its type.

    Args:
        config: The link's `type` and the fields that type's loader reads.

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


def load_link_nodes(config: Mapping[str, Any]) -> tuple[Node, ...]:
    """Builds the nodes that a link's configuration defines, such as a cooler's regenerator matrix nodes.

    Args:
        config: A link's configuration table that `load_link` accepts.

    Returns:
        The nodes, in the order the link's type defines them, or none for a type that defines none.
    """
    if config["type"] == "cooler":
        return load_matrix_nodes(config)
    return ()
