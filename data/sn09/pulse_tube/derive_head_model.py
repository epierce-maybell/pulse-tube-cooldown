"""Writes the PT1.1 head model file from the SN09 head's geometry.csv."""

import csv
import math
from pathlib import Path

from pulse_tube_cooldown.library.materials.load import load_material

GEOMETRY = Path(__file__).with_name("geometry.csv")
OUTPUT = Path(__file__).parents[3] / "src/pulse_tube_cooldown/library/heads/data/pt1_1.toml"

# Each column's tube, and its layers from the warm end as (length fraction, material, porosity, entries). A layer
# divides into its number of entries of equal length, mass, and pore volume.
COLUMNS = {
    "stage_1_column": (
        "stage1_regenerator",
        [
            ("stage1_regenerator_layer_SS", "304 SS", "matrix_porosity_screens", 1),
            ("stage1_regenerator_layer_Pb", "Pb", "matrix_porosity_spheres", 1),
        ],
    ),
    "stage_2_column": (
        "stage2_regenerator",
        [
            ("stage2_regenerator_layer_Pb", "Pb", "matrix_porosity_spheres", 1),
            ("stage2_regenerator_layer_Er3Ni", "Er3Ni", "matrix_porosity_spheres", 1),
            ("stage2_regenerator_layer_HoCu2", "HoCu2", "matrix_porosity_spheres", 1),
            ("stage2_regenerator_layer_GOS", "Gd2O2S", "matrix_porosity_spheres", 1),
        ],
    ),
}
PULSE_TUBES = ["stage1_pulse_tube", "stage2_pulse_tube"]


def read_geometry(path):
    """Returns geometry.csv's rows keyed by quantity."""
    with open(path, newline="", encoding="utf-8") as file:
        return {row["quantity"]: row for row in csv.DictReader(file)}


def value(rows, quantity):
    """Returns a geometry row's value in its csv unit."""
    return float(rows[quantity]["value"])


def bore_volume(rows, tube):
    """Returns the volume inside a tube's wall in m^3."""
    inner_radius = (value(rows, f"{tube}_outer_radius") - value(rows, "tube_wall_thickness")) * 1e-3
    return math.pi * inner_radius**2 * value(rows, f"{tube}_length") * 1e-3


def render(rows, columns=COLUMNS):
    """Returns the head model file's text, with each column's layers as in `COLUMNS`."""
    estimates = [quantity for quantity, row in rows.items() if row["source"] == "estimate"]
    pulse_tube_volumes = ", ".join(f"{bore_volume(rows, tube):.3e}" for tube in PULSE_TUBES)

    lines = [
        "# Written by data/sn09/pulse_tube/derive_head_model.py from geometry.csv.",
        f"# Estimates: {', '.join(estimates)}.",
        'name = "PT1.1"',
        f"clamp_fraction = {rows['clamp_fraction_x_R']['value']}",
        f"pulse_tube_volumes = [{pulse_tube_volumes}]  # m^3, stage 1 and stage 2",
        'source = "CAD 910-000018; Falter 2024 Fig. 2; derived in data/sn09/pulse_tube/"',
    ]
    for column, (tube, layers) in columns.items():
        column_volume = bore_volume(rows, tube)
        entries = [
            (fraction, material, porosity, count)
            for fraction, material, porosity, count in layers
            for _ in range(count)
        ]
        for index, (fraction, material, porosity, count) in enumerate(entries):
            entry_fraction = value(rows, fraction) / count
            entry_volume = column_volume * entry_fraction
            pore_fraction = value(rows, porosity)
            lines += [
                "",
                f"[[{column}]]",
                f"index = {index}",
                f'material = "{material}"',
                f"mass = {load_material(material).density * entry_volume * (1 - pore_fraction):.4g}  # kg",
                f"pore_volume = {entry_volume * pore_fraction:.3e}  # m^3",
                f"length_fraction = {entry_fraction:.6g}",
            ]
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    OUTPUT.write_text(render(read_geometry(GEOMETRY)), encoding="utf-8")
