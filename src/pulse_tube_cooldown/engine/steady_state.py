"""Steady states of a graph: node temperatures at which every free node's energy balance is zero."""

from collections.abc import Mapping

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import least_squares

from pulse_tube_cooldown.core.errors import SimulationError
from pulse_tube_cooldown.core.types import Graph
from pulse_tube_cooldown.engine.equations import calc_net_heat_flow

__all__ = ["solve_steady_state"]


def solve_steady_state(
    graph: Graph,
    temperature: NDArray[np.float64],
    /,
    *,
    fixed_boundary_nodes: tuple[int, ...] = (),
    heat_input: Mapping[str, float] | None = None,
    temperature_range: tuple[float, float] = (2.0, 1000.0),
    tolerance: float = 1e-6,
) -> NDArray[np.float64]:
    """Solves for the temperatures at which the net heat flow into every free node is zero.

    Args:
        graph: Nodes and the links acting on them.
        temperature: Temperature in K of each node: the starting point for free nodes, and the held value of fixed
            nodes.
        fixed_boundary_nodes: Indices of nodes held at their values in `temperature`.
        heat_input: Heat in W applied to named nodes, such as heater power, in addition to the links' heat flows.
        temperature_range: Lowest and highest temperatures in K the search may give a free node.
        tolerance: Largest net heat flow in W into any free node at an accepted steady state.

    Returns:
        Temperature in K of each node at steady state.

    Raises:
        SimulationError: No steady state within `tolerance` is found from `temperature`.
    """
    applied = np.zeros(len(graph.nodes), dtype=np.float64)
    for name, power in (heat_input or {}).items():
        applied[graph.get_node_index(name)] += power

    free = np.setdiff1d(np.arange(len(graph.nodes)), fixed_boundary_nodes)
    state = np.array(temperature, dtype=np.float64)

    def imbalance(log_temperature: NDArray[np.float64]) -> NDArray[np.float64]:
        state[free] = np.exp(log_temperature)
        return (calc_net_heat_flow(graph, state) + applied)[free]

    bounds = np.log(temperature_range)
    start = np.clip(np.log(state[free]), *bounds)
    solution = least_squares(imbalance, start, bounds=tuple(bounds), x_scale="jac", ftol=1e-15, xtol=1e-15, gtol=1e-15)
    residual = imbalance(solution.x)

    largest = float(np.max(np.abs(residual), initial=0.0))
    if not largest <= tolerance:
        raise SimulationError(
            f"No steady state found: the net heat flow into {graph.nodes[free[np.argmax(np.abs(residual))]].name!r} is "
            f"{largest:.3g} W, above the {tolerance:.3g} W tolerance ({solution.message})."
        )
    return state
