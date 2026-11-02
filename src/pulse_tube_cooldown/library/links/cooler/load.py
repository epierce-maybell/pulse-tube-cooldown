"""Builds a pulse-tube unit from its model and parameter files."""

import math
import os
import tomllib
from collections import Counter
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from pulse_tube_cooldown.core.errors import ConfigurationError
from pulse_tube_cooldown.core.types import Node
from pulse_tube_cooldown.library.links.cooler.head import PulseTubeHead, RegeneratorLayer
from pulse_tube_cooldown.library.materials.load import load_material
from pulse_tube_cooldown.library.toml_files import read_named_toml

__all__ = ["load_pulse_tube_head"]


def load_pulse_tube_head(
    *,
    name: str,
    stages: tuple[str, str],
    model: str,
    parameters: str | os.PathLike[str],
    directory: str | os.PathLike[str] = Path(__file__).parents[2] / "heads/data",
) -> PulseTubeHead:
    """Loads one unit, with each model layer supplying its matrix Node and gas volume.

    Args:
        name: Unit's name, used in generated matrix node names.
        stages: Names of the stage-1 and stage-2 nodes.
        model: Model's name in the library, such as "PT1.1".
        parameters: Per-unit TOML file with packing_constants, conductances, exchanger_effectiveness, and frequency.
        directory: Directory containing the head model TOML files.

    Returns:
        The unit, including its generated matrix nodes.

    Raises:
        ConfigurationError: A model is unknown, its column indices or fractions are invalid,
            or its clamp position is outside the stage-2 column.
    """
    data = read_named_toml(model, directory=directory, kind="heads")
    unit = tomllib.loads(Path(parameters).read_text(encoding="utf-8"))
    columns = tuple(
        _load_column(data[column], name, stage, column)
        for stage, column in enumerate(("stage_1_column", "stage_2_column"), start=1)
    )
    clamp_fraction = float(data["clamp_fraction"])
    if not 0 <= clamp_fraction <= 1:
        raise ConfigurationError(f"{model}: clamp_fraction must lie between 0 and 1.")
    first_volume, second_volume = data["pulse_tube_volumes"]
    first_packing, second_packing = unit["packing_constants"]
    first_conductance, second_conductance = unit["conductances"]
    return PulseTubeHead(
        name=name,
        stages=stages,
        columns=(columns[0], columns[1]),
        pulse_tube_volumes=(float(first_volume), float(second_volume)),
        clamp_fraction=clamp_fraction,
        packing_constants=(float(first_packing), float(second_packing)),
        conductances=(float(first_conductance), float(second_conductance)),
        exchanger_effectiveness=float(unit["exchanger_effectiveness"]),
        frequency=float(unit["frequency"]),
    )


def _load_column(
    entries: Sequence[Mapping[str, Any]], head_name: str, stage: int, label: str
) -> tuple[RegeneratorLayer, ...]:
    """Parses a column's indexed layers into matrix Nodes and their gas volumes.

    A layer's node is named from the head, stage, and material, numbered from the warm end where the material fills
    several layers of the column, unless the layer gives its own name.
    """
    layers = sorted(entries, key=lambda layer: layer["index"])
    if not layers or [layer["index"] for layer in layers] != list(range(len(layers))):
        raise ConfigurationError(f"{head_name} {label}: indices must run from 0 to {len(layers) - 1}, each once.")
    fractions = [float(layer["length_fraction"]) for layer in layers]
    if any(fraction <= 0 for fraction in fractions) or not math.isclose(sum(fractions), 1.0):
        raise ConfigurationError(f"{head_name} {label}: positive length fractions must sum to one.")
    materials = Counter(layer["material"] for layer in layers)
    seen: Counter[str] = Counter()
    names = []
    for layer in layers:
        material = layer["material"]
        seen[material] += 1
        number = f" {seen[material]}" if materials[material] > 1 else ""
        names.append(layer.get("name", f"{head_name} stage {stage} {material}{number}"))
    return tuple(
        RegeneratorLayer(
            node=Node(
                name=name,
                composition=((load_material(layer["material"]), float(layer["mass"])),),
            ),
            pore_volume=float(layer["pore_volume"]),
            length_fraction=float(layer["length_fraction"]),
        )
        for layer, name in zip(layers, names, strict=True)
    )
