"""Pulse-tube head designs: regenerator fillings and the geometry of the head's tubes and gas volumes."""

from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from pulse_tube_cooldown.core.types import Material
from pulse_tube_cooldown.library.materials.load import load_material
from pulse_tube_cooldown.library.toml_files import read_named_toml

__all__ = ["HeadDesignParameters", "RegeneratorLayer", "load_head_design"]


@dataclass(frozen=True, slots=True, kw_only=True)
class RegeneratorLayer:
    """One material fill of a regenerator.

    Attributes:
        material: The fill material.
        mass: Mass of the fill.
        length_fraction: Fraction of the regenerator's length.
        packing: Form of the fill, "spheres" or "screens".
        porosity: Void fraction of the fill.
        hydraulic_diameter: Hydraulic diameter of the fill's pores.
    """

    material: Material
    mass: float
    length_fraction: float
    packing: str
    porosity: float
    hydraulic_diameter: float


@dataclass(frozen=True, slots=True, kw_only=True)
class HeadDesignParameters:
    """A head design's regenerator fillings, clamp position, and tube and gas-volume geometry.

    Pairs are in stage order: regenerator or pulse tube 1, then 2.

    Attributes:
        name: Design name, such as "PT1.5".
        clamp_fraction: Clamp position as a fraction of regenerator 2's length from its warm end.
        regenerator_1: Layers from warm to cold.
        regenerator_2: Layers from warm to cold.
        regenerator_length: Length of each regenerator.
        regenerator_bore_area: Bore cross-section of each regenerator.
        regenerator_wall_area: Wall cross-section of each regenerator.
        pulse_tube_length: Length of each pulse tube, from its stage to the warm flange.
        pulse_tube_bore_area: Bore cross-section of each pulse tube.
        pulse_tube_wall_area: Wall cross-section of each pulse tube.
        heat_exchanger_void_volume: Gas volume of each stage's heat exchanger.
        dead_volume: Gas volume between the rotary valve and regenerator 1.
        reservoir_volume: Gas volume of the head's reservoirs and their hoses.
        wall_material: Material of the regenerator and pulse-tube walls.
    """

    name: str
    clamp_fraction: float
    regenerator_1: tuple[RegeneratorLayer, ...]
    regenerator_2: tuple[RegeneratorLayer, ...]
    regenerator_length: tuple[float, float]
    regenerator_bore_area: tuple[float, float]
    regenerator_wall_area: tuple[float, float]
    pulse_tube_length: tuple[float, float]
    pulse_tube_bore_area: tuple[float, float]
    pulse_tube_wall_area: tuple[float, float]
    heat_exchanger_void_volume: tuple[float, float]
    dead_volume: float
    reservoir_volume: float
    wall_material: Material


def load_head_design(name: str, /) -> HeadDesignParameters:
    """Loads a packaged head design by name, such as "PT1.5".

    Raises:
        ConfigurationError: If no file or more than one file defines the name.
    """
    table = read_named_toml(name, directory=Path(__file__).parent / "data" / "heads", kind="head designs")
    return HeadDesignParameters(
        name=table["name"],
        clamp_fraction=table["clamp_fraction"],
        regenerator_1=regenerator_layers(table["regenerator_1"]),
        regenerator_2=regenerator_layers(table["regenerator_2"]),
        regenerator_length=tuple(table["regenerator_length"]),
        regenerator_bore_area=tuple(table["regenerator_bore_area"]),
        regenerator_wall_area=tuple(table["regenerator_wall_area"]),
        pulse_tube_length=tuple(table["pulse_tube_length"]),
        pulse_tube_bore_area=tuple(table["pulse_tube_bore_area"]),
        pulse_tube_wall_area=tuple(table["pulse_tube_wall_area"]),
        heat_exchanger_void_volume=tuple(table["heat_exchanger_void_volume"]),
        dead_volume=table["dead_volume"],
        reservoir_volume=table["reservoir_volume"],
        wall_material=load_material(table["wall_material"]),
    )


def regenerator_layers(tables: list[Mapping[str, Any]]) -> tuple[RegeneratorLayer, ...]:
    """Builds regenerator layers from their tables."""
    return tuple(
        RegeneratorLayer(
            material=load_material(table["material"]),
            mass=table["mass"],
            length_fraction=table["length_fraction"],
            packing=table["packing"],
            porosity=table["porosity"],
            hydraulic_diameter=table["hydraulic_diameter"],
        )
        for table in tables
    )
