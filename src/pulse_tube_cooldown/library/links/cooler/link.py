"""One gas manifold's heads and drive, from `docs/adrs/archive/03-coupled-capacity-maps.md`.

A cooler's configuration names its warm flange, a drive file, and its heads, each with a name, two stage nodes, a head
model, and a parameter file:

    name = "Cooler"
    type = "cooler"
    warm_flange = "Warm flange"
    drive = "data/sn09/pulse_tube/drive.toml"

    [[heads]]
    name = "PT"
    stages = ["PRP", "CFP"]
    model = "PT1.1"
    parameters = "data/sn09/pulse_tube/parameters.toml"

Relative paths start from the repository root. The drive file holds the manifold's `warm_volume` in m^3, its
`fill_pressure` in Pa and `fill_temperature` in K at rest, and a `compressors` array of tables, each with a
`displacement_rate` in m^3/s and a `relief_swing` in Pa. A parameter file holds the head's stage-1 and stage-2
`packing` coefficients and needle-valve `resistance` in Pa s/m^3, and its valve `frequency` in Hz.

The configuration also defines each head's regenerator matrix nodes, which `load_matrix_nodes` builds. Each holds only
its share of the matrix and is named "<head> stage <i> <material>", numbered from 1 at the warm end where the material
fills more than one node of the section, or "<head> <name>" from the layer's `name`, numbered where the layer has more
than one node. The cooler's energy balance terms at each matrix node include the heat capacity of its pore gas at the
manifold's mean pressure.
"""

from collections import Counter
from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import Any

import numpy as np
from numpy.typing import NDArray

from pulse_tube_cooldown.core.types import EnergyBalanceTerms, Node
from pulse_tube_cooldown.library.heads import HeadModel, load_head_model
from pulse_tube_cooldown.library.links.cooler.drive import (
    Compressor,
    calc_charge,
    calc_mean_pressure,
    calc_pressure_swing,
)
from pulse_tube_cooldown.library.links.cooler.head import (
    Head,
    calc_gas_mass,
    calc_heat_flow,
    calc_mass_swing,
    calc_pore_heat_capacity,
)
from pulse_tube_cooldown.library.toml_files import read_toml

__all__ = ["Cooler", "MountedHead", "load_cooler", "load_matrix_nodes"]


@dataclass(frozen=True, slots=True, kw_only=True)
class MountedHead:
    """One head on a manifold, with the nodes along its chain.

    Attributes:
        head: The head's design and unit values.
        stages: Names of the stage-1 and stage-2 nodes.
        matrix_nodes: Names of the matrix nodes of the stage-1 and stage-2 regenerator sections, each from its warm end.
    """

    head: Head
    stages: tuple[str, str]
    matrix_nodes: tuple[tuple[str, ...], tuple[str, ...]]

    def calc_chain(self, warm_flange: str) -> tuple[str, ...]:
        """Names of the nodes along the head's chain, from the warm flange to stage 2."""
        first, second = self.matrix_nodes
        return (warm_flange, *first, self.stages[0], *second, self.stages[1])


@dataclass(frozen=True, slots=True, kw_only=True)
class Cooler:
    """One gas manifold's heads and compressors, as a link over its warm flange and each head's chain.

    Attributes:
        name: Unique link name in the graph.
        warm_flange: Name of the node at the heads' warm end, whose temperature the warm gas takes.
        heads: Heads on the manifold.
        compressors: Compressors on the manifold.
        charge: Helium mass in kg on the manifold.
        warm_volume: Gas volume in m^3 of the manifold's lines and buffers.
        chain_indices: For each head, the indices into `nodes` of its chain, derived from `warm_flange` and `heads`.
    """

    name: str
    warm_flange: str
    heads: tuple[MountedHead, ...]
    compressors: tuple[Compressor, ...]
    charge: float
    warm_volume: float
    chain_indices: tuple[tuple[int, ...], ...] = field(init=False)

    def __post_init__(self) -> None:
        """Derives each head's chain indices."""
        index = {name: position for position, name in enumerate(self.nodes)}
        chain_indices = tuple(tuple(index[name] for name in head.calc_chain(self.warm_flange)) for head in self.heads)
        object.__setattr__(self, "chain_indices", chain_indices)

    @property
    def nodes(self) -> tuple[str, ...]:
        """The warm flange, then each head's chain, each node once."""
        return tuple(dict.fromkeys(name for head in self.heads for name in head.calc_chain(self.warm_flange)))

    def calc_energy_balance_terms(self, temperature: NDArray[np.float64]) -> EnergyBalanceTerms:
        """Calculates the heat flow into each of `nodes` and the heat capacity of the gas at its temperature.

        Both follow from the manifold's mean pressure and pressure swing at these temperatures.

        Args:
            temperature: Temperatures in K of `nodes`.

        Returns:
            The heat flow in W into each of `nodes`, and the heat capacity in J/K of the gas at each node's temperature:
            pore gas at matrix nodes, heat-exchanger voids at stages, and none at the warm flange.
        """
        chains = [list(chain) for chain in self.chain_indices]
        warm = float(temperature[0])

        def head_gas_mass(pressure: float) -> float:
            return sum(
                calc_gas_mass(mounted.head.model, temperature[chain], pressure)
                for mounted, chain in zip(self.heads, chains, strict=True)
            )

        mean_pressure = calc_mean_pressure(
            head_gas_mass, self.charge, warm_volume=self.warm_volume, warm_temperature=warm
        )

        def head_mass_flow(swing: float) -> float:
            return sum(
                mounted.head.frequency
                * float(calc_mass_swing(mounted.head, temperature[chain], mean_pressure, swing).sum())
                for mounted, chain in zip(self.heads, chains, strict=True)
            )

        swing = calc_pressure_swing(head_mass_flow, mean_pressure, warm, compressors=self.compressors)

        heat_flow = np.zeros(len(temperature), dtype=np.float64)
        heat_capacity = np.zeros(len(temperature), dtype=np.float64)
        for mounted, chain in zip(self.heads, chains, strict=True):
            np.add.at(heat_flow, chain, calc_heat_flow(mounted.head, temperature[chain], mean_pressure, swing))
            np.add.at(
                heat_capacity, chain, calc_pore_heat_capacity(mounted.head.model, temperature[chain], mean_pressure)
            )
        return EnergyBalanceTerms(heat_flow=heat_flow, heat_capacity=heat_capacity)


def load_cooler(config: Mapping[str, Any]) -> Cooler:
    """Builds a cooler from its configuration table, reading its drive and head parameter files.

    Args:
        config: The link's `name`, `type`, `warm_flange`, `drive` file, and `heads` tables.

    Raises:
        KeyError: A key is missing from the table or its files.
        ValueError: `stages`, `packing`, or `resistance` does not hold exactly two entries.
        ConfigurationError: A head model is unknown or invalid.
    """
    drive = read_toml(config["drive"])
    heads = tuple(_load_head(head) for head in config["heads"])
    warm_volume = float(drive["warm_volume"])
    return Cooler(
        name=config["name"],
        warm_flange=config["warm_flange"],
        heads=heads,
        compressors=tuple(
            Compressor(
                displacement_rate=float(compressor["displacement_rate"]),
                relief_swing=float(compressor["relief_swing"]),
            )
            for compressor in drive["compressors"]
        ),
        charge=calc_charge(
            [mounted.head.model for mounted in heads],
            warm_volume=warm_volume,
            fill_pressure=float(drive["fill_pressure"]),
            fill_temperature=float(drive["fill_temperature"]),
        ),
        warm_volume=warm_volume,
    )


def load_matrix_nodes(config: Mapping[str, Any]) -> tuple[Node, ...]:
    """Builds the regenerator matrix nodes of each head on a cooler, in head and chain order.

    Args:
        config: The cooler's configuration table, as `load_cooler` reads it.

    Raises:
        KeyError: A head is missing its `name` or `model`.
        ConfigurationError: A head model is unknown or invalid.
    """
    return tuple(
        node
        for head in config["heads"]
        for section in _matrix_nodes(load_head_model(head["model"]), head["name"])
        for node in section
    )


def _load_head(config: Mapping[str, Any]) -> MountedHead:
    """Builds one head from its table and parameter file."""
    model = load_head_model(config["model"])
    parameters = read_toml(config["parameters"])
    first_stage, second_stage = config["stages"]
    first_packing, second_packing = parameters["packing"]
    first_resistance, second_resistance = parameters["resistance"]
    first_section, second_section = _matrix_nodes(model, config["name"])
    return MountedHead(
        head=Head(
            model=model,
            packing=(float(first_packing), float(second_packing)),
            resistance=(float(first_resistance), float(second_resistance)),
            frequency=float(parameters["frequency"]),
        ),
        stages=(first_stage, second_stage),
        matrix_nodes=(
            tuple(node.name for node in first_section),
            tuple(node.name for node in second_section),
        ),
    )


def _matrix_nodes(model: HeadModel, head_name: str) -> tuple[tuple[Node, ...], tuple[Node, ...]]:
    """Each regenerator section's matrix nodes from its warm end, named as the module describes."""
    sections = []
    for stage, section in enumerate(model.regenerators, start=1):
        node_counts = Counter[str]()
        for layer in section.layers:
            if layer.name is None:
                node_counts[layer.material.name] += layer.node_count
        numbered = Counter[str]()
        nodes = []
        for layer in section.layers:
            for number in range(1, layer.node_count + 1):
                if layer.name is not None:
                    name = f"{head_name} {layer.name}" + (f" {number}" if layer.node_count > 1 else "")
                else:
                    numbered[layer.material.name] += 1
                    name = f"{head_name} stage {stage} {layer.material.name}"
                    if node_counts[layer.material.name] > 1:
                        name += f" {numbered[layer.material.name]}"
                nodes.append(Node(name, ((layer.material, layer.mass / layer.node_count),)))
        sections.append(tuple(nodes))
    return sections[0], sections[1]
