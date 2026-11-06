"""Local stage coupling and stability at a test-stand operating point."""

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from pulse_tube_cooldown.core.errors import SimulationError
from pulse_tube_cooldown.core.types import Graph
from pulse_tube_cooldown.engine.equations import calc_energy_balance_terms, calc_net_heat_flow


@dataclass(frozen=True, slots=True, kw_only=True)
class Coupling:
    """Stage response with unheated internal nodes at equilibrium.

    Attributes:
        stage_heat_flow: Heat into each stage in W, in the supplied stage order.
        stage_derivative: Entry i, j is dQ_i/dT_j in W/K after the internal nodes
            re-equilibrate, with the other stage and boundary temperatures held fixed.
        matrix_growth_rate: Largest real eigenvalue in 1/s with both stages fixed.
        free_node_growth_rate: Largest real eigenvalue in 1/s with the stage heat
            inputs held constant at the values that balance stage_heat_flow.
        matrix_imbalance: Largest absolute internal-node heat imbalance in W.
    """

    stage_heat_flow: NDArray[np.float64]
    stage_derivative: NDArray[np.float64]
    matrix_growth_rate: float
    free_node_growth_rate: float
    matrix_imbalance: float


def analyze_coupling(
    graph: Graph,
    temperature: NDArray[np.float64],
    *,
    stages: tuple[str, str],
    fixed_boundary_nodes: tuple[int, ...],
    relative_step: float = 1e-4,
    tolerance: float = 1e-6,
) -> Coupling:
    """Linearizes heat flow and eliminates the unheated internal-node response.

    Args:
        graph: Stages, internal nodes, boundaries, and their thermal links.
        temperature: Positive node temperatures in K at an operating point.
        stages: The two stage names, in output order.
        fixed_boundary_nodes: Indices of nodes held at fixed temperature.
        relative_step: Fractional temperature step for centered differences.
        tolerance: Maximum accepted internal-node heat imbalance in W.

    Returns:
        Coupling and growth rates. Positive growth rates indicate instability.

    Raises:
        SimulationError: An internal node is out of equilibrium.
        numpy.linalg.LinAlgError: The internal equilibrium Jacobian is singular.
    """
    stage_indices = [graph.get_node_index(name) for name in stages]
    free = np.setdiff1d(np.arange(len(graph.nodes)), fixed_boundary_nodes)
    stage_positions = [list(free).index(index) for index in stage_indices]
    matrix_positions = np.setdiff1d(np.arange(len(free)), stage_positions)
    matrix_indices = free[matrix_positions]
    terms = calc_energy_balance_terms(graph, temperature)
    imbalance = np.abs(terms.heat_flow[matrix_indices])
    largest = float(np.max(imbalance, initial=0.0))
    if not largest <= tolerance:
        name = graph.nodes[matrix_indices[np.argmax(imbalance)]].name
        raise SimulationError(f"Internal node {name!r} has {largest:.6g} W imbalance, above {tolerance:.6g} W.")

    columns = []
    for index in free:
        step = relative_step * temperature[index]
        upper, lower = temperature.copy(), temperature.copy()
        upper[index] += step
        lower[index] -= step
        columns.append((calc_net_heat_flow(graph, upper)[free] - calc_net_heat_flow(graph, lower)[free]) / (2 * step))
    jacobian = np.column_stack(columns)
    matrix_block = jacobian[np.ix_(matrix_positions, matrix_positions)]
    matrix_response = -np.linalg.solve(matrix_block, jacobian[np.ix_(matrix_positions, stage_positions)])
    stage_derivative = (
        jacobian[np.ix_(stage_positions, stage_positions)]
        + jacobian[np.ix_(stage_positions, matrix_positions)] @ matrix_response
    )
    dynamics = jacobian / terms.heat_capacity[free, None]
    return Coupling(
        stage_heat_flow=terms.heat_flow[stage_indices],
        stage_derivative=stage_derivative,
        matrix_growth_rate=float(np.max(np.linalg.eigvals(dynamics[np.ix_(matrix_positions, matrix_positions)]).real)),
        free_node_growth_rate=float(np.max(np.linalg.eigvals(dynamics).real)),
        matrix_imbalance=largest,
    )
