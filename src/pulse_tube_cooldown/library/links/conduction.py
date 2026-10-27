"""Conduction through a solid of uniform cross-section."""

import math
from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import NDArray

from pulse_tube_cooldown.core.errors import ConfigurationError
from pulse_tube_cooldown.core.types import Material
from pulse_tube_cooldown.library.materials.load import load_material

__all__ = ["Conduction", "load_conduction"]


@dataclass(frozen=True, slots=True, kw_only=True)
class Conduction:
    """Conduction between two nodes through a solid of uniform cross-section.

    Attributes:
        name: Unique link identifier within a graph.
        nodes: Names of the two nodes the solid joins.
        material: The solid's material.
        area: Cross-sectional area in m^2.
        length: Length in m between the nodes.

    Raises:
        ConfigurationError: There are other than two nodes, or the area or length is not positive and finite.
    """

    name: str
    nodes: tuple[str, str]
    material: Material
    area: float
    length: float

    def __post_init__(self) -> None:
        """Validates the nodes and geometry."""
        if len(self.nodes) != 2:
            raise ConfigurationError(f"Conduction link {self.name!r} joins exactly two nodes; got {self.nodes}.")

        for key, value in (("area", self.area), ("length", self.length)):
            if not 0.0 < value < math.inf:
                raise ConfigurationError(f"Conduction link {self.name!r} needs a positive, finite {key}; got {value}.")

    def calc_heat_flow(self, temperature: NDArray[np.float64]) -> NDArray[np.float64]:
        """Calculates the heat flow into each of the two nodes.

        Args:
            temperature: Temperatures in K of the two nodes, in the order of `nodes`.

        Returns:
            Heat flows in W, in the order of `nodes`. Positive values indicate heat flow into the node.
        """
        first, second = self.material.conductivity_integral(temperature)
        flow = self.area / self.length * (second - first)
        return np.array([flow, -flow], dtype=np.float64)


def load_conduction(config: Mapping[str, Any]) -> Conduction:
    """Builds a conduction link from its configuration table.

    Args:
        config: The link's `name`, `type`, `nodes`, `material` name, `area` in m^2, and `length` in m.

    Raises:
        KeyError: A key is missing.
        ConfigurationError: The material is unknown, or the link is invalid.
    """
    return Conduction(
        name=config["name"],
        nodes=tuple(config["nodes"]),
        material=load_material(config["material"]),
        area=config["area"],
        length=config["length"],
    )
