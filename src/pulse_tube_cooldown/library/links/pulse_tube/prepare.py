"""Writes the packaged PT1.1 and PT1.5 designs from their geometry records."""

import csv
import math
from pathlib import Path

import tomli_w

from pulse_tube_cooldown.library.materials.load import load_material


def read_geometry(path: Path) -> dict[str, dict[str, str]]:
    """Reads geometry rows by quantity."""
    with path.open(newline="", encoding="utf-8") as file:
        return {row["quantity"]: row for row in csv.DictReader(file)}


def compile_design(name: str, rows: dict[str, dict[str, str]], source: str) -> dict:
    """Calculates layer masses in kg and length fractions, rounded to the precision of the geometry records."""

    def value(quantity: str) -> float:
        return float(rows[quantity]["value"])

    above, below = value("prp_to_clamp"), value("clamp_to_cfp")
    design = {
        "name": name,
        "clamp_fraction": float(f"{above / (above + below):.3f}"),
        "source": source,
        "regenerate": "uv run python -m pulse_tube_cooldown.library.links.pulse_tube.prepare",
    }
    for regenerator in (1, 2):
        tube = f"stage{regenerator}_regenerator"
        radius = (value(f"{tube}_outer_radius") - value("tube_wall_thickness")) * 1e-3
        volume = math.pi * radius**2 * value(f"{tube}_length") * 1e-3
        layers = []
        for quantity, row in rows.items():
            if not quantity.startswith(f"{tube}_layer_"):
                continue
            fraction = float(row["value"])
            density = load_material(row["material"]).density
            mass = density * volume * fraction * (1 - value(f"matrix_porosity_{row['packing']}"))
            layers.append(
                {"material": row["material"], "mass": float(f"{mass:.4g}"), "length_fraction": float(f"{fraction:.6g}")}
            )
        design[f"regenerator_{regenerator}"] = layers
    return design


def main() -> None:
    """Writes each head's design beside its geometry. PT1.5 takes quantities missing from its geometry from PT1.1."""
    heads = Path(__file__).parent / "data" / "heads"
    pt1_1 = read_geometry(heads / "pt1_1" / "geometry.csv")
    pt1_5 = pt1_1 | read_geometry(heads / "pt1_5" / "geometry.csv")
    for head, name, rows, source in (
        ("pt1_1", "PT1.1", pt1_1, "geometry.csv"),
        ("pt1_5", "PT1.5", pt1_5, "geometry.csv; ../pt1_1/geometry.csv"),
    ):
        (heads / head / f"{head}.toml").write_text(tomli_w.dumps(compile_design(name, rows, source)), encoding="utf-8")


if __name__ == "__main__":
    main()
