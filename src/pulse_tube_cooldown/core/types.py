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
from dataclasses import dataclass, field
from typing import Protocol

import numpy as np
from numpy.typing import NDArray

from pulse_tube_cooldown.core.errors import ConfigurationError

__all__ = ["Gas", "Graph", "Link", "Material", "Node"]


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


@dataclass(frozen=True, slots=True)
class Graph:
    """Nodes and the links acting on them.

    Node order defines the graph's temperature indices. Links refer to nodes by name, and `connectivity` translates
    those names into indices.

    Attributes:
        nodes: Nodes with unique names, in index order.
        links: Links with unique names, in index order. Each refers only to nodes in `nodes`.
        connectivity: For each link, the graph indices of its nodes in link-local order, derived from `nodes` and
            `links`. `temperature[connectivity[i]]` is the input to `links[i].calc_heat_flow`, and its output adds
            to the nodes at those same indices.

    Raises:
        ConfigurationError: A node or link name is duplicated, or a link refers to an unknown node.
    """

    nodes: tuple[Node, ...]
    links: tuple[Link, ...]
    connectivity: tuple[tuple[int, ...], ...] = field(init=False)

    def __post_init__(self) -> None:
        """Validates names and derives connectivity."""
        node_names = [node.name for node in self.nodes]

        for kind, names in (("node", node_names), ("link", [link.name for link in self.links])):
            if len(set(names)) < len(names):
                raise ConfigurationError(f"{kind.capitalize()} names are not unique: {names}.")

        for link in self.links:
            unknown = [name for name in link.nodes if name not in node_names]
            if unknown:
                raise ConfigurationError(f"Link {link.name!r} references unknown nodes {unknown}.")

        connectivity = tuple(tuple(node_names.index(name) for name in link.nodes) for link in self.links)
        object.__setattr__(self, "connectivity", connectivity)

    def get_node_index(self, name: str) -> int:
        """Returns the index of the node named `name`.

        Raises:
            ConfigurationError: No node has this name.
        """
        names = [node.name for node in self.nodes]

        if name not in names:
            raise ConfigurationError(f"No node named {name!r}; nodes are {names}.")

        return names.index(name)

    def get_link_index(self, name: str) -> int:
        """Returns the index of the link named `name`.

        Raises:
            ConfigurationError: No link has this name.
        """
        names = [link.name for link in self.links]

        if name not in names:
            raise ConfigurationError(f"No link named {name!r}; links are {names}.")

        return names.index(name)
