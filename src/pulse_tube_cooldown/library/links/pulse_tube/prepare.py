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
    """Calculates the design in SI units, rounded to the precision of the geometry records."""

    def value(quantity: str) -> float:
        return float(rows[quantity]["value"])

    def rounded(number: float) -> float:
        return float(f"{number:.4g}")

    wall = value("tube_wall_thickness") * 1e-3

    def areas(tube: str) -> tuple[float, float]:
        outer = value(f"{tube}_outer_radius") * 1e-3
        return math.pi * (outer - wall) ** 2, math.pi * (outer**2 - (outer - wall) ** 2)

    def hydraulic_diameter(row: dict[str, str], porosity: float) -> float:
        if row["packing"] == "spheres":
            return 2.0 * porosity * value(f"sphere_diameter_{row['material']}") * 1e-3 / (3.0 * (1.0 - porosity))
        return porosity * value("screen_wire_diameter") * 1e-3 / (1.0 - porosity)

    above, below = value("prp_to_clamp"), value("clamp_to_cfp")
    design = {
        "name": name,
        "clamp_fraction": float(f"{above / (above + below):.3f}"),
        "source": source,
        "regenerate": "uv run python -m pulse_tube_cooldown.library.links.pulse_tube.prepare",
        "wall_material": rows["tube_wall_material"]["material"],
        "dead_volume": rounded(value("valve_to_head_dead_volume") * 1e-6),
        "reservoir_volume": rounded(2 * (value("reservoir_volume") + value("reservoir_hose_volume")) * 1e-6),
        "heat_exchanger_void_volume": [
            rounded(value(f"stage{stage}_heat_exchanger_void_volume") * 1e-6) for stage in (1, 2)
        ],
    }
    for kind in ("regenerator", "pulse_tube"):
        tubes = [f"stage{stage}_{kind}" for stage in (1, 2)]
        design[f"{kind}_length"] = [rounded(value(f"{tube}_length") * 1e-3) for tube in tubes]
        design[f"{kind}_bore_area"] = [rounded(areas(tube)[0]) for tube in tubes]
        design[f"{kind}_wall_area"] = [rounded(areas(tube)[1]) for tube in tubes]

    for regenerator in (1, 2):
        tube = f"stage{regenerator}_regenerator"
        volume = areas(tube)[0] * value(f"{tube}_length") * 1e-3
        layers = []
        for quantity, row in rows.items():
            if not quantity.startswith(f"{tube}_layer_"):
                continue
            fraction = float(row["value"])
            porosity = value(f"matrix_porosity_{row['packing']}")
            mass = load_material(row["material"]).density * volume * fraction * (1 - porosity)
            layers.append(
                {
                    "material": row["material"],
                    "mass": rounded(mass),
                    "length_fraction": float(f"{fraction:.6g}"),
                    "packing": row["packing"],
                    "porosity": porosity,
                    "hydraulic_diameter": rounded(hydraulic_diameter(row, porosity)),
                }
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
