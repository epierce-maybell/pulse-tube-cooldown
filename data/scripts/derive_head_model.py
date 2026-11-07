"""Writes the library's PT1.1 and PT1.5 head model files from the heads' geometry records."""

import csv
import math
from pathlib import Path

from pulse_tube_cooldown.library.materials.load import load_material

ROOT = Path(__file__).resolve().parents[2]
PT1_1_GEOMETRY = "data/heads/pt1_1/geometry.csv"
PT1_5_GEOMETRY = "data/heads/pt1_5/geometry.csv"
OUTPUT = ROOT / "src/pulse_tube_cooldown/library/heads/data"

# Each regenerator section's tube, and its layers from the warm end as (length fraction, material, porosity, nodes).
REGENERATORS = {
    "stage_1_regenerator": (
        "stage1_regenerator",
        [
            ("stage1_regenerator_layer_SS", "304 SS", "matrix_porosity_screens", 1),
            ("stage1_regenerator_layer_Pb", "Pb", "matrix_porosity_spheres", 1),
        ],
    ),
    "stage_2_regenerator": (
        "stage2_regenerator",
        [
            ("stage2_regenerator_layer_Pb", "Pb", "matrix_porosity_spheres", 1),
            ("stage2_regenerator_layer_Er3Ni", "Er3Ni", "matrix_porosity_spheres", 1),
            ("stage2_regenerator_layer_HoCu2", "HoCu2", "matrix_porosity_spheres", 1),
            ("stage2_regenerator_layer_GOS", "Gd2O2S", "matrix_porosity_spheres", 8),
        ],
    ),
}
PULSE_TUBES = {"stage_1_pulse_tube": "stage1_pulse_tube", "stage_2_pulse_tube": "stage2_pulse_tube"}
HEAT_EXCHANGER_VOIDS = ("stage1_heat_exchanger_void_volume", "stage2_heat_exchanger_void_volume")
WALL_MATERIAL = "304 SS"

# Values per design from Snodgrass et al. (2024), docs/papers/pulse_tube/, with their sources.
SNODGRASS_2024 = {
    "waveform_factor": (
        math.sqrt(2),
        "Snodgrass 2024 Supp. Methods A, eq. S2: a square-wave pressure source with large omega R C carries twice a "
        "sine's acoustic power",
    ),
    "buffer_fraction": (
        0.5,
        "Snodgrass 2024 Fig. 8: large-stroke loss sets in near a stroke 2|xi|/L of 0.45 at cold impedance and 0.6 "
        "more open, read from the figure",
    ),
    "stroke_spread": (
        0.1,
        "Snodgrass 2024 Fig. 8: spread of the large-stroke onset across both needle settings, read from the figure",
    ),
}


# Each design's stage-coupling coefficient c, with its source.
STAGE_COUPLING = {
    "PT1.5": (
        3.694,
        "simulations/validation/test_stand_head/fit.txt: fitted with a and R to data/test_stand/fit_set.csv",
    ),
}


def read_geometry(path):
    """Returns a geometry record's rows keyed by quantity, each with the record's path."""
    with open(ROOT / path, newline="", encoding="utf-8") as file:
        return {row["quantity"]: {**row, "record": path} for row in csv.DictReader(file)}


def pt1_5_rows():
    """Returns PT1.5's rows, taking each quantity its record lacks from PT1.1's."""
    return read_geometry(PT1_1_GEOMETRY) | read_geometry(PT1_5_GEOMETRY)


def value(rows, quantity):
    """Returns a geometry row's value in its record's unit."""
    return float(rows[quantity]["value"])


def cite(rows, *quantities):
    """Returns the geometry rows behind a value, grouped by record, with estimates marked."""
    records = {}
    for quantity in dict.fromkeys(quantities):
        row = rows[quantity]
        mark = " (estimate)" if row["source"] == "estimate" else ""
        records.setdefault(row["record"], []).append(row["quantity"] + mark)
    return "; ".join(f"{record}: {', '.join(names)}" for record, names in records.items())


def bore_volume(rows, tube):
    """Returns the volume inside a tube's wall in m^3."""
    inner_radius = (value(rows, f"{tube}_outer_radius") - value(rows, "tube_wall_thickness")) * 1e-3
    return math.pi * inner_radius**2 * value(rows, f"{tube}_length") * 1e-3


def wall_area(rows, tube):
    """Returns the cross-sectional area of a tube's wall in m^2."""
    outer_radius = value(rows, f"{tube}_outer_radius") * 1e-3
    inner_radius = outer_radius - value(rows, "tube_wall_thickness") * 1e-3
    return math.pi * (outer_radius**2 - inner_radius**2)


def render(name, rows):
    """Returns a head model file's text from its geometry rows under PT1.1's quantity names."""
    tube_rows = ("tube_wall_thickness", "tube_wall_material")
    clamp_rows = ("prp_to_clamp", "clamp_to_cfp")
    reservoir = (value(rows, "reservoir_volume") + value(rows, "reservoir_hose_volume")) * 1e-6
    sources = {
        "clamp_fraction": cite(rows, *clamp_rows),
        "dead_volume": cite(rows, "valve_to_head_dead_volume"),
        **{key: source for key, (_, source) in SNODGRASS_2024.items()},
    }
    lines = [
        "# Written by data/scripts/derive_head_model.py from the geometry rows in [sources].",
        f'name = "{name}"',
        f"clamp_fraction = {value(rows, 'prp_to_clamp') / sum(value(rows, row) for row in clamp_rows):.3f}",
        f"dead_volume = {value(rows, 'valve_to_head_dead_volume') * 1e-6:.4g}  # m^3",
        *(f"{key} = {entry:.6g}" for key, (entry, _) in SNODGRASS_2024.items()),
    ]
    if all(row in rows for row in HEAT_EXCHANGER_VOIDS):
        sources["heat_exchanger_voids"] = cite(rows, *HEAT_EXCHANGER_VOIDS)
        voids = ", ".join(f"{value(rows, row) * 1e-6:.4e}" for row in HEAT_EXCHANGER_VOIDS)
        lines.append(f"heat_exchanger_voids = [{voids}]  # m^3, stage 1 and stage 2")
    if name in STAGE_COUPLING:
        coupling, sources["stage_coupling"] = STAGE_COUPLING[name]
        lines.append(f"stage_coupling = {coupling:.4g}")
    for section, (tube, layers) in REGENERATORS.items():
        sources[section] = cite(
            rows, f"{tube}_outer_radius", f"{tube}_length", *tube_rows, *(row for layer in layers for row in layer[::2])
        )
        lines += [
            "",
            f"[{section}]",
            f"length = {value(rows, f'{tube}_length') * 1e-3:.4g}  # m",
            f'wall_material = "{WALL_MATERIAL}"',
            f"wall_area = {wall_area(rows, tube):.4e}  # m^2",
        ]
        for index, (fraction, material, porosity, count) in enumerate(layers):
            volume = bore_volume(rows, tube) * value(rows, fraction)
            pore_fraction = value(rows, porosity)
            lines += [
                "",
                f"[[{section}.layers]]",
                f"index = {index}",
                f'material = "{material}"',
                f"mass = {load_material(material).density * volume * (1 - pore_fraction):.4g}  # kg",
                f"pore_volume = {volume * pore_fraction:.4e}  # m^3",
                f"length_fraction = {value(rows, fraction):.6g}",
            ]
            if count > 1:
                lines.append(f"node_count = {count}")
    for section, tube in PULSE_TUBES.items():
        sources[section] = cite(
            rows, f"{tube}_outer_radius", f"{tube}_length", *tube_rows, "reservoir_volume", "reservoir_hose_volume"
        )
        lines += [
            "",
            f"[{section}]",
            f"length = {value(rows, f'{tube}_length') * 1e-3:.4g}  # m",
            f"volume = {bore_volume(rows, tube):.4e}  # m^3",
            f'wall_material = "{WALL_MATERIAL}"',
            f"wall_area = {wall_area(rows, tube):.4e}  # m^2",
            f"reservoir_volume = {reservoir:.4e}  # m^3, tank and hose",
        ]
    lines += ["", "[sources]", *(f'{key} = "{source}"' for key, source in sources.items())]
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    (OUTPUT / "pt1_1.toml").write_text(render("PT1.1", read_geometry(PT1_1_GEOMETRY)), encoding="utf-8")
    (OUTPUT / "pt1_5.toml").write_text(render("PT1.5", pt1_5_rows()), encoding="utf-8")
