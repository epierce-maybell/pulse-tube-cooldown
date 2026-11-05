"""Loads head models from the library's TOML files."""

import math
import os
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from pulse_tube_cooldown.core.errors import ConfigurationError
from pulse_tube_cooldown.library.heads.model import HeadModel, MatrixLayer, PulseTube, RegeneratorSection
from pulse_tube_cooldown.library.materials.load import load_material
from pulse_tube_cooldown.library.toml_files import read_named_toml

__all__ = ["load_head_model"]


def load_head_model(name: str, /, *, directory: str | os.PathLike[str] = Path(__file__).parent / "data") -> HeadModel:
    """Loads the head model with the given name.

    Args:
        name: The model's `name` field, such as "PT1.1".
        directory: The directory whose `*.toml` files are searched. Defaults to the library's head models.

    Raises:
        ConfigurationError: No file or more than one file defines the name, or a section's layer indices, length
            fractions, or node counts are invalid.
    """
    data = read_named_toml(name, directory=directory, kind="heads")
    regenerators = tuple(
        _load_regenerator(data[section], label=f"{name} {section}")
        for section in ("stage_1_regenerator", "stage_2_regenerator")
    )
    pulse_tubes = tuple(_load_pulse_tube(data[tube]) for tube in ("stage_1_pulse_tube", "stage_2_pulse_tube"))
    return HeadModel(
        name=name,
        regenerators=(regenerators[0], regenerators[1]),
        pulse_tubes=(pulse_tubes[0], pulse_tubes[1]),
        clamp_fraction=float(data["clamp_fraction"]),
        dead_volume=float(data["dead_volume"]),
        waveform_factor=float(data["waveform_factor"]),
        buffer_fraction=float(data["buffer_fraction"]),
        stroke_spread=float(data["stroke_spread"]),
    )


def _load_regenerator(data: Mapping[str, Any], *, label: str) -> RegeneratorSection:
    """Parses a regenerator section, with its layers ordered by index from the warm end."""
    layers = sorted(data["layers"], key=lambda layer: layer["index"])
    if [layer["index"] for layer in layers] != list(range(len(layers))):
        raise ConfigurationError(f"{label}: layer indices must run from 0 to {len(layers) - 1}, each once.")
    fractions = [float(layer["length_fraction"]) for layer in layers]
    if any(fraction <= 0 for fraction in fractions) or not math.isclose(sum(fractions), 1.0):
        raise ConfigurationError(f"{label}: positive length fractions must sum to one.")
    counts = [layer.get("node_count", 1) for layer in layers]
    if any(not isinstance(count, int) or isinstance(count, bool) or count < 1 for count in counts):
        raise ConfigurationError(f"{label}: node_count must be a positive integer.")
    return RegeneratorSection(
        length=float(data["length"]),
        wall_material=load_material(data["wall_material"]),
        wall_area=float(data["wall_area"]),
        layers=tuple(
            MatrixLayer(
                material=load_material(layer["material"]),
                mass=float(layer["mass"]),
                pore_volume=float(layer["pore_volume"]),
                length_fraction=fraction,
                node_count=count,
                name=layer.get("name"),
            )
            for layer, fraction, count in zip(layers, fractions, counts, strict=True)
        ),
    )


def _load_pulse_tube(data: Mapping[str, Any]) -> PulseTube:
    """Parses a pulse tube and its reservoir."""
    return PulseTube(
        length=float(data["length"]),
        volume=float(data["volume"]),
        wall_material=load_material(data["wall_material"]),
        wall_area=float(data["wall_area"]),
        reservoir_volume=float(data["reservoir_volume"]),
    )
