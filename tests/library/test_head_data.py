"""Requirements for the head-model loader and the library's PT1.1 and PT1.5 design files."""

import pytest

from pulse_tube_cooldown.core.errors import ConfigurationError
from pulse_tube_cooldown.library.heads import load_head_model

# Hand-computed from data/heads/pt1_1/geometry.csv and data/heads/pt1_5/geometry.csv: bores and walls
# from each tube's outer radius, length, and the 0.5 mm wall; reservoirs as tank plus hose; clamp fractions as the
# clamp's distance from the stage-2 regenerator's top over the tube's length; heat-exchanger voids as recorded, none for
# PT1.1.
_DESIGNS = {
    "PT1.1": {
        "clamp_fraction": 163.4 / 244.5,
        "dead_volume": 45e-6,
        "regenerator_length": (0.1800, 0.2443),
        "regenerator_wall_area": (1.0009e-4, 5.1082e-5),
        "pulse_tube_length": (0.1810, 0.4441),
        "pulse_tube_volume": (250.765e-6, 200.906e-6),
        "pulse_tube_wall_area": (6.6759e-5, 3.8485e-5),
        "reservoir_volume": (5292e-6, 5292e-6),
        "heat_exchanger_voids": (0.0, 0.0),
        "layers": (
            (("304 SS", 0.9764, 2.8478e-4, 0.72), ("Pb", 1.1123, 6.0119e-5, 0.28)),
            (
                ("Pb", 0.2213, 1.1961e-5, 0.16),
                ("Er3Ni", 0.3532, 2.3922e-5, 0.32),
                ("HoCu2", 0.3599, 2.3922e-5, 0.32),
                ("Gd2O2S", 0.1792, 1.4951e-5, 0.20),
            ),
        ),
    },
    "PT1.5": {
        "clamp_fraction": 151.3 / 244.5,
        "dead_volume": 45e-6,
        "regenerator_length": (0.1800, 0.2445),
        "regenerator_wall_area": (1.0006e-4, 5.1051e-5),
        "pulse_tube_length": (0.1810, 0.4440),
        "pulse_tube_volume": (250.765e-6, 200.861e-6),
        "pulse_tube_wall_area": (6.6759e-5, 3.8485e-5),
        "reservoir_volume": (5280e-6, 5280e-6),
        "heat_exchanger_voids": (39.3e-6, 10.8e-6),
        "layers": (
            (("304 SS", 0.9758, 2.8459e-4, 0.72), ("Pb", 1.1116, 6.0081e-5, 0.28)),
            (
                ("Pb", 0.2212, 1.1956e-5, 0.16),
                ("Er3Ni", 0.3530, 2.3911e-5, 0.32),
                ("HoCu2", 0.3597, 2.3911e-5, 0.32),
                ("Gd2O2S", 0.1791, 1.4945e-5, 0.20),
            ),
        ),
    },
}


@pytest.mark.parametrize("name", _DESIGNS)
def test_design_geometry_matches_hand_computed_values(name):
    expected = _DESIGNS[name]

    model = load_head_model(name)

    assert model.name == name
    assert model.clamp_fraction == pytest.approx(expected["clamp_fraction"], abs=5e-4)
    assert model.dead_volume == pytest.approx(expected["dead_volume"])
    assert model.heat_exchanger_voids == pytest.approx(expected["heat_exchanger_voids"])
    for stage in (0, 1):
        regenerator, tube = model.regenerators[stage], model.pulse_tubes[stage]
        assert regenerator.length == pytest.approx(expected["regenerator_length"][stage])
        assert regenerator.wall_area == pytest.approx(expected["regenerator_wall_area"][stage], rel=1e-3)
        assert regenerator.wall_material.name == "304 SS"
        assert tube.length == pytest.approx(expected["pulse_tube_length"][stage])
        assert tube.volume == pytest.approx(expected["pulse_tube_volume"][stage], rel=1e-3)
        assert tube.wall_area == pytest.approx(expected["pulse_tube_wall_area"][stage], rel=1e-3)
        assert tube.wall_material.name == "304 SS"
        assert tube.reservoir_volume == pytest.approx(expected["reservoir_volume"][stage], rel=1e-3)


@pytest.mark.parametrize("name", _DESIGNS)
def test_matrix_filling_matches_hand_computed_values(name):
    model = load_head_model(name)

    for regenerator, layers in zip(model.regenerators, _DESIGNS[name]["layers"], strict=True):
        assert [layer.material.name for layer in regenerator.layers] == [layer[0] for layer in layers]
        assert [layer.mass for layer in regenerator.layers] == pytest.approx([layer[1] for layer in layers], rel=1e-3)
        assert [layer.pore_volume for layer in regenerator.layers] == pytest.approx(
            [layer[2] for layer in layers], rel=1e-3
        )
        assert [layer.length_fraction for layer in regenerator.layers] == pytest.approx([layer[3] for layer in layers])


def _write_design(directory, layer):
    tube = 'length = 0.2\nvolume = 1e-4\nwall_material = "304 SS"\nwall_area = 1e-5\nreservoir_volume = 5e-3\n'
    regenerator = 'length = 0.2\nwall_material = "304 SS"\nwall_area = 1e-4\n'
    (directory / "design.toml").write_text(
        'name = "Test"\nclamp_fraction = 0.5\ndead_volume = 4e-5\nwaveform_factor = 1.4\n'
        "buffer_fraction = 0.5\nstroke_spread = 0.1\n"
        f"[stage_1_regenerator]\n{regenerator}[[stage_1_regenerator.layers]]\n{layer}"
        f"[stage_2_regenerator]\n{regenerator}"
        '[[stage_2_regenerator.layers]]\nindex = 0\nmaterial = "Pb"\nmass = 1.0\npore_volume = 1e-5\n'
        "length_fraction = 1.0\n"
        f"[stage_1_pulse_tube]\n{tube}[stage_2_pulse_tube]\n{tube}",
        encoding="utf-8",
    )


@pytest.mark.parametrize(
    ("field", "value"), [("index", "2"), ("length_fraction", "0.5"), ("node_count", "0"), ("node_count", "1.5")]
)
def test_invalid_layer_is_rejected(tmp_path, field, value):
    fields = {"index": "0", "material": '"Pb"', "mass": "1.0", "pore_volume": "1e-5", "length_fraction": "1.0"}
    fields["node_count"] = "1"
    fields[field] = value
    _write_design(tmp_path, "".join(f"{key} = {entry}\n" for key, entry in fields.items()))

    with pytest.raises(ConfigurationError, match="stage_1_regenerator"):
        load_head_model("Test", directory=tmp_path)


def test_layer_defaults_to_one_unnamed_node(tmp_path):
    _write_design(tmp_path, 'index = 0\nmaterial = "Pb"\nmass = 1.0\npore_volume = 1e-5\nlength_fraction = 1.0\n')

    layer = load_head_model("Test", directory=tmp_path).regenerators[0].layers[0]

    assert layer.node_count == 1
    assert layer.name is None


def test_design_without_voids_loads_them_as_zero(tmp_path):
    _write_design(tmp_path, 'index = 0\nmaterial = "Pb"\nmass = 1.0\npore_volume = 1e-5\nlength_fraction = 1.0\n')

    assert load_head_model("Test", directory=tmp_path).heat_exchanger_voids == (0.0, 0.0)


@pytest.mark.parametrize(("line", "expected"), [("", 0.0), ("stage_coupling = 1.5\n", 1.5)])
def test_stage_coupling_loads_from_design_and_defaults_to_zero(tmp_path, line, expected):
    _write_design(tmp_path, 'index = 0\nmaterial = "Pb"\nmass = 1.0\npore_volume = 1e-5\nlength_fraction = 1.0\n')
    design = tmp_path / "design.toml"
    design.write_text(line + design.read_text(encoding="utf-8"), encoding="utf-8")

    assert load_head_model("Test", directory=tmp_path).stage_coupling == expected
