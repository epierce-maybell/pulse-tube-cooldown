"""Conduction through a solid of uniform cross-section."""

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import NDArray

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
    """

    name: str
    nodes: tuple[str, str]
    material: Material
    area: float
    length: float

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
        ValueError: `nodes` does not hold exactly two entries.
        ConfigurationError: The material is unknown.
    """
    first, second = config["nodes"]
    return Conduction(
        name=config["name"],
        nodes=(first, second),
        material=load_material(config["material"]),
        area=config["area"],
        length=config["length"],
    )
