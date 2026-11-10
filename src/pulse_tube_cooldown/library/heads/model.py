"""Shared head design, from `docs/adrs/archive/03-coupled-capacity-maps.md`."""

from dataclasses import dataclass

from pulse_tube_cooldown.core.types import Material

__all__ = ["HeadModel", "MatrixLayer", "PulseTube", "RegeneratorSection"]


@dataclass(frozen=True, slots=True, kw_only=True)
class MatrixLayer:
    """One layer of matrix filling along a regenerator section.

    Attributes:
        material: Matrix material.
        mass: Matrix mass in kg.
        pore_volume: Gas volume in m^3 between the matrix particles.
        length_fraction: Fraction of the section's length the layer fills.
        node_count: Number of equal matrix nodes the layer divides into.
        name: Name for the layer's nodes, or None to name them from the head, stage, and material.
    """

    material: Material
    mass: float
    pore_volume: float
    length_fraction: float
    node_count: int = 1
    name: str | None = None


@dataclass(frozen=True, slots=True, kw_only=True)
class RegeneratorSection:
    """One regenerator section: warm flange to stage 1, or stage 1 to stage 2.

    Attributes:
        length: Section length in m.
        wall_material: Material of the tube wall.
        wall_area: Cross-sectional area of the tube wall in m^2.
        layers: Matrix layers ordered from the section's warm end.
    """

    length: float
    wall_material: Material
    wall_area: float
    layers: tuple[MatrixLayer, ...]


@dataclass(frozen=True, slots=True, kw_only=True)
class PulseTube:
    """One pulse tube, from the warm flange to its stage, and its reservoir.

    Attributes:
        length: Tube length in m.
        volume: Gas volume inside the tube in m^3.
        wall_material: Material of the tube wall.
        wall_area: Cross-sectional area of the tube wall in m^2.
        reservoir_volume: Gas volume in m^3 beyond the tube's needle valve, at the warm flange.
    """

    length: float
    volume: float
    wall_material: Material
    wall_area: float
    reservoir_volume: float


@dataclass(frozen=True, slots=True, kw_only=True)
class HeadModel:
    """The design of one head model.

    Attributes:
        name: Model name, such as "PT1.1".
        regenerators: Stage-1 and stage-2 regenerator sections.
        pulse_tubes: Stage-1 and stage-2 pulse tubes.
        clamp_fraction: Clamp position x_R along the stage-2 section, as a fraction of its length from stage 1.
        dead_volume: Gas volume in m^3 between the rotary valve and the head, at the warm flange.
        waveform_factor: Pressure amplitude over half the swing, kappa.
        buffer_fraction: Fraction of a pulse tube's gas that buffers its stroke, lambda.
        stroke_spread: Spread of the stroke as a fraction of a pulse tube's gas, sigma.
        heat_exchanger_voids: Gas volumes in m^3 of the stage-1 and stage-2 cold heat exchangers, outside the tubes.
        stage_coupling: Coefficient c of stage 1's heat flow c E_2 x (1 - x), for stage 2's acoustic power E_2 and
            gross-cooling fraction x.
    """

    name: str
    regenerators: tuple[RegeneratorSection, RegeneratorSection]
    pulse_tubes: tuple[PulseTube, PulseTube]
    clamp_fraction: float
    dead_volume: float
    waveform_factor: float
    buffer_fraction: float
    stroke_spread: float
    heat_exchanger_voids: tuple[float, float] = (0.0, 0.0)
    stage_coupling: float = 0.0
