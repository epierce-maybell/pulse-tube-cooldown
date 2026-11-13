"""Loading an energy-flow pulse-tube head from its link table, its design, and the packaged parameters."""

import itertools
import tomllib
from collections.abc import Mapping
from pathlib import Path
from typing import Any, NamedTuple

import numpy as np

from pulse_tube_cooldown.core.errors import ConfigurationError
from pulse_tube_cooldown.core.types import Node
from pulse_tube_cooldown.library.links.pulse_tube.design import RegeneratorLayer, load_head_design
from pulse_tube_cooldown.library.links.pulse_tube.energy_flow.head import PulseTubeHead
from pulse_tube_cooldown.library.links.pulse_tube.energy_flow.manifold import Manifold

__all__ = ["load_energy_flow_pulse_tube_head"]


class Section(NamedTuple):
    """A regenerator section: its fill, length, mass, and warm-end position along its regenerator."""

    layer: RegeneratorLayer
    length: float
    mass: float
    start: float


def load_energy_flow_pulse_tube_head(config: Mapping[str, Any]) -> PulseTubeHead:
    """Loads a head from its link table, its packaged design, and the packaged energy-flow parameters.

    Args:
        config: An `energy_flow_pulse_tube_head` link table:

            name = "Head"
            type = "energy_flow_pulse_tube_head"
            design = "PT1.5"
            stages = ["PRP", "CFP"]
            clamp = "RGP"
            contact_resistance = 2.0
            warm_flange_temperature = 295.0
            valve_frequency = 1.1

            [manifold]
            compressor = "F-70"
            head_count = 1
            compressor_count = 2
            fill_pressure = 1.7e6
            supply_volume = 12.8e-3
            return_volume = 15.8e-3

            `clamp` and `contact_resistance` are given together, or neither for an unclamped head.

    Raises:
        ConfigurationError: Only one of `clamp` and `contact_resistance` is given, or the resistance is not positive.
    """
    clamp = config.get("clamp")
    resistance = config.get("contact_resistance")
    if (clamp is None) != (resistance is None) or (resistance is not None and resistance <= 0):
        raise ConfigurationError("A clamped head needs both a clamp node and a positive contact_resistance.")

    parameters = tomllib.loads((Path(__file__).parent / "data" / "parameters.toml").read_text(encoding="utf-8"))
    design = load_head_design(config["design"])
    coefficients = parameters["designs"][design.name]
    manifold = config["manifold"]
    compressor = parameters["compressors"][manifold["compressor"]]
    cells = parameters["cells_per_layer"]

    clamp_position = design.clamp_fraction * design.regenerator_length[1]
    regenerator_1 = regenerator_sections(design.regenerator_1, design.regenerator_length[0], cells)
    regenerator_2 = regenerator_sections(design.regenerator_2, design.regenerator_length[1], cells, clamp_position)
    upper = [section for section in regenerator_2 if section.start < clamp_position]
    lower = [section for section in regenerator_2 if section.start >= clamp_position]
    sections = regenerator_1 + regenerator_2
    counts = (len(regenerator_1), len(regenerator_2))

    names = (f"{config['name']} regenerator node {number}" for number in itertools.count(1))

    def nodes(sections: list[Section]) -> tuple[Node, ...]:
        return tuple(Node(next(names), ((section.layer.material, section.mass),)) for section in sections)

    return PulseTubeHead(
        name=config["name"],
        stages=tuple(config["stages"]),
        regenerator_plate=clamp,
        warm_flange_temperature=float(config["warm_flange_temperature"]),
        clamp_contact_conductance=0.0 if resistance is None else 1.0 / resistance,
        effective_impedance=coefficients["effective_impedance"],
        stage_1_split=coefficients["stage_1_split"],
        valve_frequency=float(config["valve_frequency"]),
        nusselt_coeff=parameters["nusselt"]["coefficient"],
        nusselt_peclet_exp=parameters["nusselt"]["peclet_exponent"],
        nusselt_porosity_exp=parameters["nusselt"]["porosity_exponent"],
        axial_conductivity_coeff=parameters["axial_conductivity"]["coefficient"],
        axial_conductivity_peclet_exp=parameters["axial_conductivity"]["peclet_exponent"],
        axial_conductivity_porosity_exp=parameters["axial_conductivity"]["porosity_exponent"],
        manifold=Manifold(
            head_count=manifold["head_count"],
            compressor_count=manifold["compressor_count"],
            displacement_rate=compressor["displacement_rate"],
            volumetric_efficiency=coefficients["volumetric_efficiency"],
            relief_amplitude=compressor["relief_amplitude"],
            fill_pressure=manifold["fill_pressure"],
            supply_volume=manifold["supply_volume"],
            return_volume=manifold["return_volume"],
        ),
        regenerator_1=nodes(regenerator_1),
        regenerator_2_upper=nodes(upper),
        regenerator_2_lower=nodes(lower),
        section_length=np.array([section.length for section in sections]),
        section_flow_area=np.repeat(design.regenerator_bore_area, counts),
        section_porosity=np.array([section.layer.porosity for section in sections]),
        section_hydraulic_diameter=np.array([section.layer.hydraulic_diameter for section in sections]),
        section_wall_area=np.repeat(design.regenerator_wall_area, counts),
        pulse_tube_wall_area=np.array(design.pulse_tube_wall_area),
        pulse_tube_length=np.array(design.pulse_tube_length),
        pulse_tube_bore_area=np.array(design.pulse_tube_bore_area),
        pulse_tube_segments=parameters["pulse_tube_segments"],
        heat_exchanger_void_volume=np.array(design.heat_exchanger_void_volume),
        dead_volume=design.dead_volume,
        reservoir_volume=design.reservoir_volume,
        wall_material=design.wall_material,
    )


def regenerator_sections(
    layers: tuple[RegeneratorLayer, ...], length: float, cells_per_layer: int, cut: float | None = None
) -> list[Section]:
    """Return a regenerator's sections from warm to cold.

    Each layer splits into equal cells, and the cell holding `cut` splits there.
    """
    sections = []
    start = 0.0
    for layer in layers:
        layer_length = layer.length_fraction * length
        bounds = start + layer_length * np.arange(cells_per_layer + 1) / cells_per_layer
        if cut is not None and bounds[0] < cut < bounds[-1]:
            bounds = np.unique(np.append(bounds, cut))
        for warm, cold in itertools.pairwise(bounds):
            sections.append(Section(layer, cold - warm, layer.mass * (cold - warm) / layer_length, warm))
        start += layer_length
    return sections
