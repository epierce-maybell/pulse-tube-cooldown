"""Conditions that end phases and procedures, built from configuration.

A condition's configuration is a table holding its `condition`, a key of `CONDITIONS`, and the arguments that
condition's factory takes after the graph:

    until = { condition = "gas_pushing_begins", node = "CFP", threshold = 70.0 }

Each factory returns a function of (time, temperature) that is positive before its event, zero at the event, and
negative after it.
"""

from collections.abc import Callable, Mapping

import numpy as np
from numpy.typing import NDArray

from pulse_tube_cooldown.core.types import Graph

__all__ = ["CONDITIONS", "gas_pushing_begins"]


def gas_pushing_begins(graph: Graph, node: str, threshold: float) -> Callable[[float, NDArray[np.float64]], float]:
    """Gas pushing begins once a node cools to a threshold temperature.

    Args:
        graph: A graph that holds the node, with the node ordering the condition will see.
        node: The name of the node.
        threshold: The node's temperature at which gas pushing begins, in K.

    Raises:
        ConfigurationError: The graph has no node with this name.
    """
    idx = graph.get_node_index(node)
    return lambda time, temperature: temperature[idx] - threshold


CONDITIONS: Mapping[str, Callable[..., Callable[[float, NDArray[np.float64]], float]]] = {
    "gas_pushing_begins": gas_pushing_begins,
}
