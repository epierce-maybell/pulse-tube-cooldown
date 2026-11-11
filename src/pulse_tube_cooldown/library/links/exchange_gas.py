"""Exchange gas in the vacuum can, conducting heat across the gaps between facing surfaces."""

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import NDArray

from pulse_tube_cooldown.core.errors import ConfigurationError
from pulse_tube_cooldown.core.types import EnergyBalanceTerms, Gas
from pulse_tube_cooldown.library.gases.load import load_gas
from pulse_tube_cooldown.library.links.gas_gap_switch import calc_gap_conductance

__all__ = ["ExchangeGas", "FacingSurfaces", "load_exchange_gas"]


@dataclass(frozen=True, slots=True, kw_only=True)
class FacingSurfaces:
    """Two facing surfaces with exchange gas in the gap between them.

    Attributes:
        nodes: Names of the two nodes the surfaces belong to.
        area: Facing surface area in m^2.
        gap: Distance between the surfaces in m.
    """

    nodes: tuple[str, str]
    area: float
    gap: float


@dataclass(frozen=True, slots=True, kw_only=True)
class ExchangeGas:
    """Gas at one pressure throughout the vacuum can, set by a cubic in one node's temperature.

    Attributes:
        name: Unique link identifier within a graph.
        nodes: Names of every node in a facing pair, and the pressure node.
        gas: Gas with conductivity, accommodation coefficient, heat-capacity ratio, and molar mass.
        pressure_node: Node whose temperature sets the pressure.
        pressure_coefficients: Cubic coefficients in Pa / K^3, Pa / K^2, Pa / K, and Pa, highest power first.
        surfaces: Facing surface pairs the gas conducts across.
    """

    name: str
    nodes: tuple[str, ...]
    gas: Gas
    pressure_node: str
    pressure_coefficients: tuple[float, float, float, float]
    surfaces: tuple[FacingSurfaces, ...]

    def calc_pressure(self, temperature: NDArray[np.float64]) -> float:
        """Calculates gas pressure in Pa, clamped at zero, from node temperatures in K in `nodes` order."""
        pressure_temperature = temperature[self.nodes.index(self.pressure_node)]
        return max(float(np.polyval(self.pressure_coefficients, pressure_temperature)), 0.0)

    def calc_energy_balance_terms(self, temperature: NDArray[np.float64]) -> EnergyBalanceTerms:
        """Calculates heat flows in W into each node at node temperatures in K, both in `nodes` order."""
        flows = np.zeros(len(self.nodes), dtype=np.float64)
        pressure = self.calc_pressure(temperature)
        if pressure == 0.0:
            return EnergyBalanceTerms(heat_flow=flows)
        for surfaces in self.surfaces:
            first, second = (self.nodes.index(node) for node in surfaces.nodes)
            mean = float((temperature[first] + temperature[second]) / 2)
            conductance = calc_gap_conductance(self.gas, mean, pressure, area=surfaces.area, gap=surfaces.gap)
            flow = conductance * (temperature[second] - temperature[first])
            flows[first] += flow
            flows[second] -= flow
        return EnergyBalanceTerms(heat_flow=flows)


def load_exchange_gas(config: Mapping[str, Any]) -> ExchangeGas:
    """Builds an exchange-gas link from its TOML link table.

    Args:
        config: The link's name, nodes, gas name, pressure node, four pressure coefficients, and surfaces. Each
            surface pair gives its two nodes, area, and gap. Physical values are SI.

    Returns:
        The link with its gas loaded.

    Raises:
        KeyError: A required key is missing.
        ValueError: A surface pair or the coefficients have the wrong number of entries.
        ConfigurationError: The gas is unknown, or the pressure node or a surface node is not among `nodes`.
    """
    nodes = tuple(config["nodes"])
    a, b, c, d = config["pressure_coefficients"]
    surfaces = []
    for entry in config["surfaces"]:
        first, second = entry["nodes"]
        surfaces.append(FacingSurfaces(nodes=(first, second), area=entry["area"], gap=entry["gap"]))
    referenced = {config["pressure_node"], *(node for pair in surfaces for node in pair.nodes)}
    if missing := sorted(referenced - set(nodes)):
        raise ConfigurationError(f"Link {config['name']!r} refers to {missing!r}, which are not among {nodes!r}.")
    return ExchangeGas(
        name=config["name"],
        nodes=nodes,
        gas=load_gas(config["gas"]),
        pressure_node=config["pressure_node"],
        pressure_coefficients=(a, b, c, d),
        surfaces=tuple(surfaces),
    )
