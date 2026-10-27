"""Thermal properties and heat balances of a graph at given node temperatures.

Temperature arrays hold one value in K per node, in node index order.
"""

from collections.abc import Sequence

import numpy as np
from numpy.typing import NDArray

from pulse_tube_cooldown.core.types import Graph, Node

__all__ = ["calc_heat_capacity", "calc_link_heat_flows", "calc_net_heat_flow"]


def calc_heat_capacity(nodes: Sequence[Node], temperature: NDArray[np.float64]) -> NDArray[np.float64]:
    """Calculates the heat capacity of each node at its temperature.

    Args:
        nodes: Nodes in index order.
        temperature: Temperature in K of each node.

    Returns:
        Heat capacity in J/K of each node.
    """
    return np.array(
        [
            sum(mass * material.specific_heat(temperature[i : i + 1])[0] for material, mass in node.composition)
            for i, node in enumerate(nodes)
        ],
        dtype=np.float64,
    )


def calc_link_heat_flows(graph: Graph, temperature: NDArray[np.float64]) -> NDArray[np.float64]:
    """Calculates the heat flow from each link into each node.

    Args:
        graph: Nodes and the links acting on them.
        temperature: Temperature in K of each node.

    Returns:
        Heat flow in W, with one row per link and one column per node, in index order. Positive values indicate heat
        flow into the node; entries are NaN where a link does not touch a node, matching a time slice of the Result's
        `heat_flow`.
    """
    flows = np.full((len(graph.links), len(graph.nodes)), np.nan, dtype=np.float64)

    for row, (link, indices) in enumerate(zip(graph.links, graph.connectivity, strict=True)):
        flows[row, list(indices)] = link.calc_heat_flow(temperature[list(indices)])

    return flows


def calc_net_heat_flow(graph: Graph, temperature: NDArray[np.float64]) -> NDArray[np.float64]:
    """Calculates the total heat flow into each node from all links.

    Args:
        graph: Nodes and the links acting on them.
        temperature: Temperature in K of each node.

    Returns:
        Heat flow in W into each node. Positive values indicate heat flow into the node. A NaN from a link propagates
        to its nodes.
    """
    flows = calc_link_heat_flows(graph, temperature)
    net = np.zeros(len(graph.nodes), dtype=np.float64)

    for row, indices in enumerate(graph.connectivity):
        net[list(indices)] += flows[row, list(indices)]

    return net
