"""Rechecks archived coupling cases: run with python -m simulations.validation.test_stand_steady_state.check_cases."""

import dataclasses
import json
from pathlib import Path

import numpy as np

import pulse_tube_cooldown as ptc
from pulse_tube_cooldown.core.types import Graph
from simulations.validation.test_stand_steady_state.coupling import analyze_coupling


def main() -> None:
    """Reports derivatives and stability for the archived, equilibrated cooler-only cases."""
    model = ptc.configure_model("data/test_stand/fridge.toml")
    cooler = model.graphs[0].links[0]
    mounted = cooler.heads[0]
    chain = [model.graphs[0].get_node_index(name) for name in mounted.calc_chain(cooler.warm_flange)]
    cases = json.loads(Path(__file__).with_name("coupling_cases.json").read_text(encoding="utf-8"))
    for case in cases:
        head = dataclasses.replace(mounted.head, packing=tuple(case["packing"]), resistance=tuple(case["resistance"]))
        link = dataclasses.replace(cooler, heads=(dataclasses.replace(mounted, head=head),))
        graph = Graph(nodes=model.nodes, links=(link,))
        temperature = model.initial_temperature.copy()
        temperature[chain] = np.array(case["temperature"])
        response = analyze_coupling(
            graph, temperature, stages=mounted.stages, fixed_boundary_nodes=model.fixed_boundary_nodes
        )
        print(
            json.dumps(
                {
                    "case": case["label"],
                    "stage_derivative_W_per_K": response.stage_derivative.tolist(),
                    "matrix_growth_rate_per_s": response.matrix_growth_rate,
                    "free_node_growth_rate_per_s": response.free_node_growth_rate,
                    "matrix_imbalance_W": response.matrix_imbalance,
                }
            ),
            flush=True,
        )


if __name__ == "__main__":
    main()
