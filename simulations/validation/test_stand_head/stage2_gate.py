"""Stage 2's slope at T1 = 50 K and T2 = 3.4 K on the cold matrix branch.

Run from the repository root with `python simulations/validation/test_stand_head/stage2_gate.py`.
It writes `stage2_gate.txt` beside this file.

The run uses the test-stand template, `data/test_stand/test_stand.toml`, with PT1.5's heat-exchanger voids and the
mean-specific-heat matrix capacity at node boundaries, at a1 = 0.15, R1 = 1e10 Pa s/m^3, and c = 0. The stages are
held and the matrix nodes are solved. From the template's a2 = 0.05 and R2 = 3e10 Pa s/m^3, each path moves a2 to
its value and then raises R2, starting each solve from the previous state, which keeps the matrix on its cold branch.
The couplings and growth rates are those of `analyze_coupling` in
`simulations/validation/test_stand_steady_state/coupling.py`.
"""

import dataclasses
import sys
import warnings
from pathlib import Path

import numpy as np

import pulse_tube_cooldown as ptc
from pulse_tube_cooldown.core.types import Graph
from pulse_tube_cooldown.engine.steady_state import solve_steady_state

sys.path.insert(0, str(Path(__file__).parents[1] / "test_stand_steady_state"))
from coupling import analyze_coupling  # noqa: E402


def with_stage_2(graph: Graph, packing: float, resistance: float) -> Graph:
    """The graph with every head's a2 and R2 replaced."""

    def replace(link):
        if not hasattr(link, "heads"):
            return link
        heads = tuple(
            dataclasses.replace(
                mounted,
                head=dataclasses.replace(
                    mounted.head,
                    packing=(mounted.head.packing[0], packing),
                    resistance=(mounted.head.resistance[0], resistance),
                ),
            )
            for mounted in link.heads
        )
        return dataclasses.replace(link, heads=heads)

    return dataclasses.replace(graph, links=tuple(replace(link) for link in graph.links))


def main(
    paths: tuple[tuple[float, tuple[float, ...]], ...] = (
        (0.05, (3e10,)),
        (0.038, (3e10, 4.0e10, 4.3e10)),
        (0.034, (3e10, 4.0e10, 4.4e10, 4.8e10)),
    ),
    stage_temperature: tuple[float, float] = (50.0, 3.4),
    matrix_start: tuple[float, ...] = (217, 93, 45, 26, 13, 4.7, 3.9, 3.7, 3.6, 3.55, 3.5, 3.45, 3.42),
    tolerance: float = 1e-5,
) -> None:
    """Writes the stage heater powers, couplings, and growth rates at each path's R2 values.

    Args:
        paths: Each a2 with the R2 values in Pa s/m^3 to report, in the order reached.
        stage_temperature: T1 and T2 in K.
        matrix_start: Matrix temperatures in K, warm end first, from which the first solve starts.
        tolerance: Largest net heat flow in W into any matrix node at an accepted state.
    """
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    model = ptc.configure_model("data/test_stand/test_stand.toml")
    template = model.graphs[0]
    stages = ("Stage 1", "Stage 2")
    stage_indices = [template.get_node_index(name) for name in stages]
    held = (*model.fixed_boundary_nodes, *stage_indices)
    matrix = np.setdiff1d(np.arange(len(model.nodes)), held)
    start = np.full(len(model.nodes), 295.0)
    start[stage_indices] = stage_temperature
    start[matrix] = matrix_start

    def solve(a2, r2, state):
        return solve_steady_state(with_stage_2(template, a2, r2), state, fixed_boundary_nodes=held, tolerance=tolerance)

    origin = solve(0.05, 3e10, start)
    lines = [
        f"Stage 2's slope on the cold matrix branch at T1 = {stage_temperature[0]} K, T2 = {stage_temperature[1]} K.",
        "a1 = 0.15, R1 = 1e10 Pa s/m^3, c = 0. P1 and P2 are the stage heater powers in W. The couplings",
        "are in W/K with heat into the node positive. The growth rates in 1/s are with the stages held (matrix) and",
        "with the heater powers held (free).",
        "",
        f"{'a2':>6} {'R2':>8} {'P1':>7} {'P2':>6} {'dQ1/dT1':>8} {'dQ1/dT2':>8} {'dQ2/dT1':>8} {'dQ2/dT2':>8} "
        f"{'matrix':>9} {'free':>9}",
    ]
    for a2, resistances in paths:
        state = origin
        for step in np.linspace(0.05, a2, 9)[1:]:
            state = solve(step, 3e10, state)
        previous = 3e10
        for r2 in resistances:
            for step in np.linspace(previous, r2, 9)[1:]:
                state = solve(a2, step, state)
            previous = r2
            coupling = analyze_coupling(
                with_stage_2(template, a2, r2),
                state,
                stages=stages,
                fixed_boundary_nodes=model.fixed_boundary_nodes,
                tolerance=tolerance,
            )
            heater, slope = -coupling.stage_heat_flow, coupling.stage_derivative
            lines.append(
                f"{a2:6.3f} {r2:8.3g} {heater[0]:+7.2f} {heater[1]:+6.3f} {slope[0, 0]:+8.3f} {slope[0, 1]:+8.3f} "
                f"{slope[1, 0]:+8.4f} {slope[1, 1]:+8.3f} {coupling.matrix_growth_rate:+9.2e} "
                f"{coupling.free_node_growth_rate:+9.2e}"
            )
    Path(__file__).with_suffix(".txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
