"""Gray-body radiation between an inner surface and the outer surface enclosing it."""

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import NDArray
from scipy.constants import Stefan_Boltzmann

from pulse_tube_cooldown.core.errors import ConfigurationError
from pulse_tube_cooldown.core.types import Surface
from pulse_tube_cooldown.library.surfaces.load import load_surface

__all__ = ["Radiation", "load_radiation"]


@dataclass(frozen=True, slots=True, kw_only=True)
class Radiation:
    """Gray-body radiation between two nodes, the inner one enclosed by the outer one.

    Parallel plates are the case of equal areas.

    Attributes:
        name: Unique link identifier within a graph.
        nodes: Names of the inner and outer nodes, in that order.
        area: Radiating areas in m^2 of the inner and outer surfaces.
        surface: Finishes of the inner and outer surfaces.

    Raises:
        ConfigurationError: The inner area is larger than the outer area.
    """

    name: str
    nodes: tuple[str, str]
    area: tuple[float, float]
    surface: tuple[Surface, Surface]

    def __post_init__(self) -> None:
        """Checks that the outer surface can enclose the inner one."""
        inner, outer = self.area
        if inner > outer:
            raise ConfigurationError(
                f"Radiation link {self.name!r} has inner area {inner} larger than outer area {outer}."
            )

    def calc_heat_flow(self, temperature: NDArray[np.float64]) -> NDArray[np.float64]:
        """Calculates the heat flow into each of the two nodes.

        Args:
            temperature: Temperatures in K of the inner and outer nodes.

        Returns:
            Heat flows in W, in the order of `nodes`. Positive values indicate heat flow into the node.
        """
        t_inner, t_outer = temperature
        a_inner, a_outer = self.area
        e_inner = self.surface[0].emissivity(temperature[:1])[0]
        e_outer = self.surface[1].emissivity(temperature[1:])[0]

        flow = (
            Stefan_Boltzmann
            * a_inner
            * (t_outer**4 - t_inner**4)
            / (1.0 / e_inner + a_inner / a_outer * (1.0 / e_outer - 1.0))
        )
        return np.array([flow, -flow], dtype=np.float64)


def load_radiation(config: Mapping[str, Any]) -> Radiation:
    """Builds a radiation link from its configuration table.

    Args:
        config: The link's `name`, `type`, `nodes` as `[inner, outer]`, `area` in m^2 and `surface` names, each as
            `[inner, outer]`.

    Raises:
        KeyError: A key is missing.
        ValueError: `nodes`, `area`, or `surface` does not hold exactly two entries.
        ConfigurationError: A surface is unknown, or the link is invalid.
    """
    inner_node, outer_node = config["nodes"]
    inner_area, outer_area = config["area"]
    inner_surface, outer_surface = config["surface"]
    return Radiation(
        name=config["name"],
        nodes=(inner_node, outer_node),
        area=(inner_area, outer_area),
        surface=(load_surface(inner_surface), load_surface(outer_surface)),
    )
