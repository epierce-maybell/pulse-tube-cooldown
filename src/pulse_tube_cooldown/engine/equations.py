"""Thermal properties and energy balances of a graph at given node temperatures.

Each node's energy balance is C dT/dt = Q, for its heat capacity C and the heat flow Q into it. Temperature arrays hold
one value in K per node, in node index order.
"""

from collections.abc import Sequence

import numpy as np
from numpy.typing import NDArray

from pulse_tube_cooldown.core.types import EnergyBalanceTerms, Graph, Node

__all__ = [
    "calc_energy_balance_terms",
    "calc_heat_capacity",
    "calc_link_heat_flows",
    "calc_net_heat_flow",
    "calc_temperature_rate_of_change",
]


def calc_heat_capacity(nodes: Sequence[Node], temperature: NDArray[np.float64]) -> NDArray[np.float64]:
    """Calculates the heat capacity of each node at its temperature.

    Args:
        nodes: Nodes in index order.
        temperature: Temperature in K of each node.

    Returns:
        Heat capacity in J/K of each node.
    """
    return np.array([node.calc_heat_capacity(float(temperature[i])) for i, node in enumerate(nodes)], dtype=np.float64)


def calc_energy_balance_terms(graph: Graph, temperature: NDArray[np.float64]) -> EnergyBalanceTerms:
    """Assembles each node's energy balance from its own heat capacity and its links' terms.

    Args:
        graph: Nodes and the links acting on them.
        temperature: Temperature in K of each node.

    Returns:
        The net heat flow into each node, and each node's heat capacity plus that of the links' gas there, in node index
        order. A NaN from a link propagates to its nodes.
    """
    heat_flow = np.zeros(len(graph.nodes), dtype=np.float64)
    heat_capacity = calc_heat_capacity(graph.nodes, temperature)

    for link, indices in zip(graph.links, graph.connectivity, strict=True):
        terms = link.calc_energy_balance_terms(temperature[list(indices)])
        heat_flow[list(indices)] += terms.heat_flow
        heat_capacity[list(indices)] += terms.heat_capacity

    return EnergyBalanceTerms(heat_flow=heat_flow, heat_capacity=heat_capacity)


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
        flows[row, list(indices)] = link.calc_energy_balance_terms(temperature[list(indices)]).heat_flow

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
    return calc_energy_balance_terms(graph, temperature).heat_flow


def calc_temperature_rate_of_change(
    graph: Graph, fixed_boundary_nodes: tuple[int, ...], time: float, temperature: NDArray[np.float64]
) -> NDArray[np.float64]:
    """Calculates the rate of change of each node's temperature.

    Args:
        graph: Nodes and the links acting on them.
        fixed_boundary_nodes: Indices of nodes held at constant temperature.
        time: Time in s. The rate depends only on temperature; `time` matches the ODE solver's call signature.
        temperature: Temperature in K of each node.

    Returns:
        Rate of change in K/s of each node's temperature from its energy balance, and zero at fixed boundary nodes.
        Links still exchange heat with fixed boundary nodes, so their neighbors' rates include it.
    """
    terms = calc_energy_balance_terms(graph, temperature)
    rate_of_change = terms.heat_flow / terms.heat_capacity
    rate_of_change[list(fixed_boundary_nodes)] = 0.0
    return rate_of_change
