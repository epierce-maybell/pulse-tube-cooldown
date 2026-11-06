"""The test-stand template's cold root with the matrix capacity at node boundaries from task 3.3.

Run from the repository root with
`python simulations/validation/test_stand_head/boundary_capacity_stability.py`.
It writes `boundary_capacity_stability.txt` beside this file.

The run starts from the unstable cold root in `simulations/validation/test_stand_steady_state/report.txt`, with no
heater power, and marches every free node in time until the heat flows settle. A settled state is a stable steady
state, which the steady-state solver then polishes. Its growth rate is the largest real part of the eigenvalues of the
free nodes' heat-flow Jacobian, by centred differences of 1e-4 of each temperature, with each row divided by that
node's heat capacity. It runs with the heat-exchanger voids zeroed, as in the report, and with PT1.5's.
"""

import dataclasses
from pathlib import Path

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

import pulse_tube_cooldown as ptc
from pulse_tube_cooldown.core.types import Graph
from pulse_tube_cooldown.engine.equations import calc_energy_balance_terms, calc_net_heat_flow
from pulse_tube_cooldown.engine.steady_state import solve_steady_state

# The report's cold root, in K, in the template's node order.
REPORT_ROOT = np.array(
    [295, 295, 43.21967056, 3.13403475, 222.70565885, 84.18527372, 37.72324017, 24.56445317, 13.01793574, 7.91567591]
    + [7.84192735, 7.75419686, 7.64807061, 7.51697982, 7.35054704, 7.13079306, 5.38051750]
)


def without_voids(graph: Graph) -> Graph:
    """The graph with every head's heat-exchanger voids set to zero."""

    def strip(link):
        if not hasattr(link, "heads"):
            return link
        heads = tuple(
            dataclasses.replace(
                mounted,
                head=dataclasses.replace(
                    mounted.head, model=dataclasses.replace(mounted.head.model, heat_exchanger_voids=(0.0, 0.0))
                ),
            )
            for mounted in link.heads
        )
        return dataclasses.replace(link, heads=heads)

    return dataclasses.replace(graph, links=tuple(strip(link) for link in graph.links))


def settle(graph: Graph, start: NDArray[np.float64], fixed: tuple[int, ...], duration: float) -> NDArray[np.float64]:
    """Temperatures in K after marching the free nodes from `start` for `duration` seconds."""
    free = np.setdiff1d(np.arange(len(start)), fixed)
    state = start.copy()

    def rate(_, temperature):
        state[free] = temperature
        terms = calc_energy_balance_terms(graph, state)
        return terms.heat_flow[free] / terms.heat_capacity[free]

    solution = solve_ivp(rate, (0.0, duration), start[free], method="BDF", rtol=1e-6, atol=1e-6)
    state[free] = solution.y[:, -1]
    return state


def growth_rate(graph: Graph, temperature: NDArray[np.float64], fixed: tuple[int, ...]) -> complex:
    """The eigenvalue in 1/s of the free nodes' linearized dynamics with the largest real part."""
    free = np.setdiff1d(np.arange(len(temperature)), fixed)
    columns = []
    for index in free:
        step = 1e-4 * temperature[index]
        upper, lower = temperature.copy(), temperature.copy()
        upper[index] += step
        lower[index] -= step
        columns.append((calc_net_heat_flow(graph, upper)[free] - calc_net_heat_flow(graph, lower)[free]) / (2 * step))
    capacity = calc_energy_balance_terms(graph, temperature).heat_capacity[free]
    eigenvalues = np.linalg.eigvals(np.column_stack(columns) / capacity[:, None])
    return complex(eigenvalues[np.argmax(eigenvalues.real)])


def main() -> None:
    """Writes each void setting's settled root and growth rate."""
    model = ptc.configure_model("data/test_stand/fridge.toml")
    fixed = model.fixed_boundary_nodes
    names = [node.name for node in model.nodes]
    lines = ["Test-stand cold root with mean-specific-heat boundary capacity (task 3.3)", ""]
    for label, graph in (("Voids zeroed", without_voids(model.graphs[0])), ("PT1.5 voids", model.graphs[0])):
        settled = settle(graph, REPORT_ROOT, fixed, duration=3e5)
        root = solve_steady_state(graph, settled, fixed_boundary_nodes=fixed)
        growth = growth_rate(graph, root, fixed)
        lines += [
            f"{label}: growth rate {growth.real:+.5f} {growth.imag:+.5f}j /s",
            *(f"  {name:<22} {temperature:9.4f} K" for name, temperature in zip(names, root, strict=True)),
            "",
        ]
    Path(__file__).with_suffix(".txt").write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    main()
