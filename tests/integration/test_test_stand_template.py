"""The test-stand template builds with its boundaries at room temperature."""

import numpy as np

import pulse_tube_cooldown as ptc


def test_template_holds_its_boundaries_at_room_temperature():
    model = ptc.configure_model("data/test_stand/fridge.toml")
    names = [node.name for node in model.nodes]

    assert [names[index] for index in model.fixed_boundary_nodes] == ["Warm flange", "Vacuum can"]
    np.testing.assert_array_equal(model.initial_temperature[list(model.fixed_boundary_nodes)], 295.0)
