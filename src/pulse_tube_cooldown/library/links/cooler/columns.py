"""Expands physical regenerator layers into equal-volume matrix nodes."""

import math
from collections import Counter
from collections.abc import Mapping, Sequence
from typing import Any

from pulse_tube_cooldown.core.errors import ConfigurationError
from pulse_tube_cooldown.core.types import Node
from pulse_tube_cooldown.library.links.cooler.head import RegeneratorLayer
from pulse_tube_cooldown.library.materials.load import load_material


def load_column(entries: Sequence[Mapping[str, Any]], *, head_name: str, stage: int) -> tuple[RegeneratorLayer, ...]:
    """Parses indexed physical layers into their matrix nodes and gas volumes.

    Each layer's positive integer node_count defaults to one. Its nodes share its
    mass, pore volume, and length equally. Generated names count each material's
    nodes from the warm end; a split named layer numbers its own nodes.

    Args:
        entries: Physical layers with indices, materials, masses, pore volumes, and length fractions.
        head_name: Head's name in generated node names and configuration errors.
        stage: Column's stage number, 1 or 2.

    Raises:
        ConfigurationError: Indices, length fractions, or node counts are invalid.
    """
    label = f"{head_name} stage_{stage}_column"
    layers = sorted(entries, key=lambda layer: layer["index"])
    if not layers or [layer["index"] for layer in layers] != list(range(len(layers))):
        raise ConfigurationError(f"{label}: indices must run from 0 to {len(layers) - 1}, each once.")
    fractions = [float(layer["length_fraction"]) for layer in layers]
    if any(fraction <= 0 for fraction in fractions) or not math.isclose(sum(fractions), 1.0):
        raise ConfigurationError(f"{label}: positive length fractions must sum to one.")
    counts = [layer.get("node_count", 1) for layer in layers]
    if any(not isinstance(count, int) or isinstance(count, bool) or count < 1 for count in counts):
        raise ConfigurationError(f"{label}: node_count must be a positive integer.")
    materials: Counter[str] = Counter()
    for layer, count in zip(layers, counts, strict=True):
        materials[layer["material"]] += count
    seen: Counter[str] = Counter()
    nodes = []
    for layer, count in zip(layers, counts, strict=True):
        material_name = layer["material"]
        material = load_material(material_name)
        for index in range(count):
            seen[material_name] += 1
            if "name" in layer:
                name = layer["name"] + (f" {index + 1}" if count > 1 else "")
            else:
                number = f" {seen[material_name]}" if materials[material_name] > 1 else ""
                name = f"{head_name} stage {stage} {material_name}{number}"
            nodes.append(
                RegeneratorLayer(
                    node=Node(name=name, composition=((material, float(layer["mass"]) / count),)),
                    pore_volume=float(layer["pore_volume"]) / count,
                    length_fraction=float(layer["length_fraction"]) / count,
                )
            )
    return tuple(nodes)
