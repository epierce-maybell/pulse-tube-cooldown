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

from collections import defaultdict
from collections.abc import Callable, Mapping
from dataclasses import dataclass, field
from typing import Any, Protocol

import numpy as np
import xarray as xr
from numpy.typing import NDArray

from pulse_tube_cooldown.core.errors import ConfigurationError

__all__ = ["Gas", "Graph", "Link", "Material", "Model", "Node", "Phase", "Procedure", "Result"]


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


@dataclass(frozen=True, slots=True)
class Phase:
    """A stage of a run during which one set of links is active.

    Attributes:
        name: Unique phase identifier within a procedure, such as "precool".
        links: The links active during this phase. A link may reappear in later phases with different parameters.
        until: Maps time in s and node temperatures in K to a value that is positive while the phase continues, zero
            where it ends, and negative afterwards. None for the final phase, which continues until the procedure
            stops.
    """

    name: str
    links: tuple[Link, ...]
    until: Callable[[float, NDArray[np.float64]], float] | None


@dataclass(frozen=True, slots=True)
class Procedure:
    """The sequence of phases making up a run.

    Attributes:
        phases: Phases in the order they run. Every phase except the last has an `until` condition. A link name
            refers to the same nodes, in the same order, in every phase.
        stop_condition: Maps time in s and node temperatures in K to a value that is positive while the run continues
            and zero where it ends, during any phase. None when only the solver's time limit ends the run.

    Raises:
        ConfigurationError: There are no phases, phase names repeat, a phase other than the last lacks `until`, the
            last phase has `until`, or a link name refers to different nodes or node orders in different phases.
    """

    phases: tuple[Phase, ...]
    stop_condition: Callable[[float, NDArray[np.float64]], float] | None

    def __post_init__(self) -> None:
        """Validates phase names, end conditions, and link identities."""
        if not self.phases:
            raise ConfigurationError("A procedure needs at least one phase.")

        names = [phase.name for phase in self.phases]
        repeated = sorted({name for name in names if names.count(name) > 1})
        if repeated:
            raise ConfigurationError(f"Phase names {repeated} repeat; phases are {names}.")

        *leading, last = self.phases
        unending = [phase.name for phase in leading if phase.until is None]
        if unending:
            raise ConfigurationError(f"Phases {unending} need an `until` condition because later phases follow them.")
        if last.until is not None:
            raise ConfigurationError(f"Final phase {last.name!r} has an `until` condition; no phase follows it.")

        nodes_by_link: defaultdict[str, set[tuple[str, ...]]] = defaultdict(set)
        for phase in self.phases:
            for link in phase.links:
                nodes_by_link[link.name].add(link.nodes)

        inconsistent = {name: sorted(nodes) for name, nodes in nodes_by_link.items() if len(nodes) > 1}
        if inconsistent:
            raise ConfigurationError(
                "Each link needs the same node names, in the same order, in every phase. These links differ in names "
                f"or order between phases: {inconsistent}."
            )


@dataclass(frozen=True, slots=True)
class Model:
    """Everything needed to simulate one run.

    Attributes:
        nodes: Nodes with unique names, in index order.
        procedure: The phases to run and the condition that stops the run.
        initial_temperature: Temperature in K of each node at the start of the run.
        fixed_boundary_nodes: Indices of nodes held at their initial temperatures.
        solver_options: Options passed to the ODE solver, such as the integration method.
        graphs: For each phase in `procedure`, the graph of `nodes` and that phase's links, derived from the other
            fields.

    Raises:
        ConfigurationError: A phase's graph is invalid, such as a link referring to an unknown node.
    """

    nodes: tuple[Node, ...]
    procedure: Procedure
    initial_temperature: NDArray[np.float64]
    fixed_boundary_nodes: tuple[int, ...]
    solver_options: Mapping[str, Any]
    graphs: tuple[Graph, ...] = field(init=False)

    def __post_init__(self) -> None:
        """Builds one graph per phase."""
        graphs = []
        for phase in self.procedure.phases:
            try:
                graphs.append(Graph(nodes=self.nodes, links=phase.links))
            except ConfigurationError as error:
                raise ConfigurationError(f"Phase {phase.name!r}: {error}") from error

        object.__setattr__(self, "graphs", tuple(graphs))


@dataclass(frozen=True, slots=True)
class Result:
    """The outcome of simulating a model.

    Attributes:
        status: How the run ended.
        model: The simulated model.
        data: The simulated history. Dimensions are `time`, `node`, and `link`, with coordinates of the model's node
            names and every link name used in any phase. Variables are `phase(time)` (phase names),
            `temperature(time, node)` in K, `heat_flow(time, link, node)` in W into the node (NaN where the link does
            not touch the node or is inactive), and `connected(link, node)`. A phase switch repeats its time, once for
            each phase. Attributes record status, units, and run provenance.
    """

    status: str
    model: Model
    data: xr.Dataset
