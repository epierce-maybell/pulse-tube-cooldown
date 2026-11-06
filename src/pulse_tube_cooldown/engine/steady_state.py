"""Steady states of a graph: node temperatures at which every free node's energy balance is zero."""

from collections.abc import Mapping

import numpy as np
from numpy.typing import NDArray

from pulse_tube_cooldown.core.errors import SimulationError
from pulse_tube_cooldown.core.types import Graph
from pulse_tube_cooldown.engine.equations import calc_energy_balance_terms

__all__ = ["solve_steady_state"]


def solve_steady_state(
    graph: Graph,
    temperature: NDArray[np.float64],
    /,
    *,
    fixed_boundary_nodes: tuple[int, ...] = (),
    heat_input: Mapping[str, float] | None = None,
    temperature_floor: float = 2.1,
    tolerance: float = 1e-6,
    max_steps: int = 500,
) -> NDArray[np.float64]:
    """Solves for the temperatures at which the net heat flow into every free node is zero.

    The free nodes relax along their energy balances, C dT/dt = Q, in log temperature. Each step is a linearly implicit
    Euler step whose time step grows as the heat flows fall, so the final steps are Newton steps. The steps follow the
    dynamics only while the time step is shorter than the growth time of any unstable mode, so the steady state reached
    can be unstable. Successive solves continue from one another by passing each solution as the next `temperature`.

    Args:
        graph: Nodes and the links acting on them.
        temperature: Temperature in K of each node: the starting point for free nodes, and the held value of fixed
            nodes. Free nodes start no colder than `temperature_floor`.
        fixed_boundary_nodes: Indices of nodes held at their values in `temperature`.
        heat_input: Heat in W applied to named nodes, such as heater power, in addition to the links' heat flows.
        temperature_floor: Lowest temperature in K a free node may reach. The default is the lower limit of the helium
            properties.
        tolerance: Largest net heat flow in W into any free node at an accepted steady state.
        max_steps: Most relaxation steps taken before the search is abandoned.

    Returns:
        Temperature in K of each node at steady state.

    Raises:
        SimulationError: A step takes a free node below `temperature_floor`, or no steady state within `tolerance` is
            reached in `max_steps` steps.
    """
    applied = np.zeros(len(graph.nodes), dtype=np.float64)
    for name, power in (heat_input or {}).items():
        applied[graph.get_node_index(name)] += power

    free = np.setdiff1d(np.arange(len(graph.nodes)), fixed_boundary_nodes)
    start = np.array(temperature, dtype=np.float64)
    start[free] = np.maximum(start[free], temperature_floor)

    def with_free(log_temperature: NDArray[np.float64]) -> NDArray[np.float64]:
        state = start.copy()
        state[free] = np.exp(log_temperature)
        return state

    def evaluate(log_temperature: NDArray[np.float64]) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Net heat flow in W into each free node, and its heat capacity per unit log temperature, C T, in J."""
        state = with_free(log_temperature)
        terms = calc_energy_balance_terms(graph, state)
        return (terms.heat_flow + applied)[free], (terms.heat_capacity * state)[free]

    def heat_flow_jacobian(
        log_temperature: NDArray[np.float64], heat_flow: NDArray[np.float64], step: float = 1e-7
    ) -> NDArray[np.float64]:
        """Derivative in W of each free node's heat flow with respect to each free node's log temperature."""
        columns = []
        for index in range(len(log_temperature)):
            shifted = log_temperature.copy()
            shifted[index] += step
            columns.append((evaluate(shifted)[0] - heat_flow) / step)
        return np.column_stack(columns)

    log_temperature = np.log(start[free])
    heat_flow, log_capacity = evaluate(log_temperature)
    if np.max(np.abs(heat_flow), initial=0.0) <= tolerance:
        return with_free(log_temperature)

    # The first time step moves the fastest-relaxing node by 0.1 in log temperature; no step moves any node by more
    # than 0.5.
    relaxation_rate = np.linalg.norm(heat_flow / log_capacity)
    time_step = 0.1 / np.max(np.abs(heat_flow / log_capacity))
    for _ in range(max_steps):
        jacobian = heat_flow_jacobian(log_temperature, heat_flow)
        change = np.linalg.solve(np.diag(log_capacity / time_step) - jacobian, heat_flow)
        change *= min(1.0, 0.5 / np.max(np.abs(change)))
        log_temperature = log_temperature + change

        if np.min(log_temperature) < np.log(temperature_floor):
            node = graph.nodes[free[np.argmin(log_temperature)]].name
            raise SimulationError(f"No steady state above the {temperature_floor:.3g} K floor: {node!r} reaches it.")

        heat_flow, log_capacity = evaluate(log_temperature)
        largest = np.max(np.abs(heat_flow))
        if largest <= tolerance:
            return with_free(log_temperature)

        previous_rate, relaxation_rate = relaxation_rate, np.linalg.norm(heat_flow / log_capacity)
        time_step *= previous_rate / relaxation_rate

    raise SimulationError(
        f"No steady state found in {max_steps} steps: the net heat flow into "
        f"{graph.nodes[free[np.argmax(np.abs(heat_flow))]].name!r} is {largest:.3g} W, above the {tolerance:.3g} W "
        "tolerance."
    )
