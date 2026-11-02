"""Tests for the library's head model files against the head geometry they were derived from."""

import itertools
from pathlib import Path

import pytest

import pulse_tube_cooldown
from pulse_tube_cooldown.library.materials.load import load_material
from pulse_tube_cooldown.library.toml_files import read_named_toml

_HEADS = Path(pulse_tube_cooldown.__file__).parent / "library/heads/data"


@pytest.fixture(scope="module")
def pt1_1():
    return read_named_toml("PT1.1", directory=_HEADS, kind="heads")


def test_stage_1_pore_volumes_sum_to_geometry_record(pt1_1):
    # data/sn09/pulse_tube/geometry.csv, stage1_regenerator_gas_volume: 344.9 cm^3.
    total = sum(layer["pore_volume"] for layer in pt1_1["stage_1_column"])

    assert total == pytest.approx(344.9e-6, abs=0.05e-6)


def test_stage_2_pore_volumes_sum_to_geometry_record(pt1_1):
    # data/sn09/pulse_tube/geometry.csv, the stage-2 regenerator gas volumes on either side of the clamp: 49.2 cm^3
    # and 25.6 cm^3, each rounded to 0.1 cm^3.
    total = sum(layer["pore_volume"] for layer in pt1_1["stage_2_column"])

    assert total == pytest.approx(74.8e-6, abs=0.1e-6)


# Hand-computed from geometry.csv and the library's densities; stage 1 is DESIGN "Cooler configuration"'s example.
@pytest.mark.parametrize(
    ("column", "material", "mass", "pore_volume"),
    [
        ("stage_1_column", "304 SS", 0.9764, 284.8e-6),
        ("stage_1_column", "Pb", 1.1123, 60.12e-6),
        ("stage_2_column", "Pb", 0.2213, 11.96e-6),
        ("stage_2_column", "Er3Ni", 0.3532, 23.92e-6),
        ("stage_2_column", "HoCu2", 0.3599, 23.92e-6),
        ("stage_2_column", "Gd2O2S", 0.1792, 14.95e-6),
    ],
)
def test_layers_match_hand_computed_values(pt1_1, column, material, mass, pore_volume):
    layers = [layer for layer in pt1_1[column] if layer["material"] == material]

    assert sum(layer["mass"] for layer in layers) == pytest.approx(mass, rel=1e-3)
    assert sum(layer["pore_volume"] for layer in layers) == pytest.approx(pore_volume, rel=1e-3)


@pytest.mark.parametrize("column", ["stage_1_column", "stage_2_column"])
def test_columns_take_the_design_form(pt1_1, column):
    layers = pt1_1[column]

    assert sorted(layer["index"] for layer in layers) == list(range(len(layers)))
    assert sum(layer["length_fraction"] for layer in layers) == pytest.approx(1.0)
    for layer in layers:
        load_material(layer["material"])


def test_head_values_match_geometry_record(pt1_1):
    # data/sn09/pulse_tube/geometry.csv: clamp_fraction_x_R, and the stage-1 and stage-2 pulse-tube gas volumes, each
    # rounded to 0.1 cm^3.
    assert pt1_1["clamp_fraction"] == 0.658
    assert pt1_1["pulse_tube_volumes"] == pytest.approx([59.7e-6, 200.9e-6], abs=0.1e-6)


def test_clamp_lies_in_the_hocu2_layer(pt1_1):
    # geometry.csv, stage2_regenerator_layer_HoCu2 notes.
    layers = sorted(pt1_1["stage_2_column"], key=lambda layer: layer["index"])
    ends = itertools.accumulate(layer["length_fraction"] for layer in layers)

    clamp_layer = next(layer for layer, end in zip(layers, ends, strict=True) if pt1_1["clamp_fraction"] < end)
    assert clamp_layer["material"] == "HoCu2"
