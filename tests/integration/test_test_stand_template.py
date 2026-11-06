"""The test-stand template builds with its boundaries at room temperature and reaches a steady state."""

import numpy as np

import pulse_tube_cooldown as ptc
from pulse_tube_cooldown.engine.equations import calc_net_heat_flow
from pulse_tube_cooldown.engine.steady_state import solve_steady_state


def test_template_holds_its_boundaries_at_room_temperature():
    model = ptc.configure_model("data/test_stand/fridge.toml")
    names = [node.name for node in model.nodes]

    assert [names[index] for index in model.fixed_boundary_nodes] == ["Warm flange", "Vacuum can"]
    np.testing.assert_array_equal(model.initial_temperature[list(model.fixed_boundary_nodes)], 295.0)


def test_template_reaches_a_steady_state_at_a_load_map_point():
    model = ptc.configure_model("data/test_stand/fridge.toml")
    graph = model.graphs[0]
    # A cold profile, warm end first along each regenerator section: flange, can, stages, then the matrix nodes.
    start = np.array([295, 295, 43, 3.1, 223, 84, 38, 25, 13, 7.9, 7.8, 7.8, 7.6, 7.5, 7.4, 7.1, 5.4], dtype=float)
    # S0009's 2026-05-14 load map: 20 W on stage 1 and 1 W on stage 2.
    powers = {"Stage 1": 20.0, "Stage 2": 1.0}

    temperature = solve_steady_state(
        graph, start, fixed_boundary_nodes=model.fixed_boundary_nodes, heat_input=powers, tolerance=1e-5
    )

    imbalance = calc_net_heat_flow(graph, temperature)
    for name, power in powers.items():
        imbalance[graph.get_node_index(name)] += power
    free = np.setdiff1d(np.arange(len(model.nodes)), model.fixed_boundary_nodes)
    assert np.max(np.abs(imbalance[free])) <= 1e-5
    np.testing.assert_array_equal(temperature[list(model.fixed_boundary_nodes)], 295.0)
