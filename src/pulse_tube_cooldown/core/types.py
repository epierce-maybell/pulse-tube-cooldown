"""Component types for thermal networks.

All physical quantities are in SI units, for example:

- temperature: K
- time: s
- mass: kg
- heat flow: W
- specific heat: J/(kg K)
- thermal conductivity: W/(m K)
- length: m; area: m^2; volume: m^3
- pressure: Pa
- molar mass: kg/mol
"""

from collections.abc import Callable
from dataclasses import dataclass
from typing import Protocol

import numpy as np
from numpy.typing import NDArray

__all__ = ["Gas", "Link", "Material", "Node"]


@dataclass(frozen=True, slots=True, kw_only=True)
class Material:
    """A substance with temperature-dependent thermal properties.

    Attributes:
        name: Identifier used by configuration files, such as "OFHC Cu RRR100".
        specific_heat: Maps temperatures in K to specific heats in J/(kg K), elementwise.
        thermal_conductivity: Maps temperatures in K to thermal conductivities in W/(m K), elementwise.
    """

    name: str
    specific_heat: Callable[[NDArray[np.float64]], NDArray[np.float64]]
    thermal_conductivity: Callable[[NDArray[np.float64]], NDArray[np.float64]]


@dataclass(frozen=True, slots=True, kw_only=True)
class Gas:
    """A gas used as a heat-transfer medium.

    Attributes:
        name: Identifier used by configuration files, such as "He-4".
        molar_mass: Molar mass in kg/mol.
        heat_capacity_ratio: Ratio of specific heats at constant pressure and constant volume.
    """

    name: str
    molar_mass: float
    heat_capacity_ratio: float


@dataclass(frozen=True, slots=True)
class Node:
    """A thermal mass with a single temperature.

    Attributes:
        name: Unique node identifier within a graph, such as "PRP".
        composition: Pairs of a material and its mass in kg.
    """

    name: str
    composition: tuple[tuple[Material, float], ...]


class Link(Protocol):
    """Anything that transfers heat to or from its nodes.

    A link connects two nodes, couples several nodes through a shared medium, or acts as a source or sink at its
    nodes. It operates only on its own nodes, in the order given by `nodes`.
    """

    @property
    def name(self) -> str:
        """Unique link identifier within a graph."""
        ...

    @property
    def nodes(self) -> tuple[str, ...]:
        """Names of the nodes this link touches, in link-local order."""
        ...

    def calc_heat_flow(self, temperature: NDArray[np.float64]) -> NDArray[np.float64]:
        """Calculates the heat flow into each of this link's nodes.

        Args:
            temperature: Temperatures in K of this link's nodes, in the order of `nodes`.

        Returns:
            Heat flows in W, in the order of `nodes`. Positive values indicate heat flow into the node.
        """
        ...
