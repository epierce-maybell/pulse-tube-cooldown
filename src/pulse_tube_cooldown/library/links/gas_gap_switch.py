"""Sealed gas-gap heat switches with charcoal adsorption and series brackets."""

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import NDArray
from scipy.constants import R

from pulse_tube_cooldown.core.errors import ConfigurationError
from pulse_tube_cooldown.core.types import EnergyBalanceTerms, Gas, Material
from pulse_tube_cooldown.library.gases.load import load_gas
from pulse_tube_cooldown.library.materials.load import load_material

__all__ = ["GasGapSwitch", "SwitchBracket", "calc_gap_conductance", "load_gas_gap_switch"]


def calc_gap_conductance(gas: Gas, temperature: float, pressure: float, *, area: float, gap: float) -> float:
    """Calculates gap conductance in W/K with gas properties at the mean temperature.

    Args:
        gas: Gas with conductivity, accommodation coefficient, heat-capacity ratio, and molar mass.
        temperature: Mean surface temperature in K.
        pressure: Gas pressure in Pa.
        area: Facing surface area in m^2.
        gap: Distance between the surfaces in m.

    Returns:
        Conductance combining viscous and molecular resistances in series, in W/K.
    """
    at_temperature = np.array([temperature], dtype=np.float64)
    conductivity = float(gas.thermal_conductivity(at_temperature)[0])
    accommodation = float(gas.accommodation_coeff(at_temperature)[0])
    ratio = gas.heat_capacity_ratio
    molecular = accommodation * (ratio + 1) / (ratio - 1) * np.sqrt(R / (8 * np.pi * gas.molar_mass * temperature))
    molecular *= pressure
    return float(area * molecular / (1 + gap * molecular / conductivity))


@dataclass(frozen=True, slots=True, kw_only=True)
class SwitchBracket:
    """A solid bracket in series with the switch's gas gap.

    Attributes:
        material: Bracket material, with thermal conductivity in W/(m K).
        area: Cross-sectional area along the heat path in m^2.
        length: Heat-path length in m.
    """

    material: Material
    area: float
    length: float


@dataclass(frozen=True, slots=True, kw_only=True)
class GasGapSwitch:
    """A single or double switch with one sealed gas inventory and charcoal temperature.

    Attributes:
        name: Unique link identifier within a graph.
        nodes: Two or three node names in heat-path order. Each adjacent pair bounds one gap.
        gas: Gas properties and Langmuir adsorption constants on charcoal.
        fill_pressure: Absolute gas pressure when filled, in Pa.
        fill_temperature: Gas temperature when filled, in K.
        gas_volume: Total sealed gas volume in m^3, divided equally between the gaps.
        gap: Width of each gas gap in m.
        area: Facing surface area of each gap in m^2.
        charcoal_mass: Total adsorbent mass in kg.
        brackets: Solid brackets in series with each gas gap.
        charcoal_node: Node whose temperature sets adsorption. Required for a double switch.
            None uses the mean endpoint temperature of a single switch.

    Raises:
        ConfigurationError: The switch has other than two or three nodes, a double has no charcoal node,
            or the charcoal node is outside the switch.
    """

    name: str
    nodes: tuple[str, ...]
    gas: Gas
    fill_pressure: float
    fill_temperature: float
    gas_volume: float
    gap: float
    area: float
    charcoal_mass: float
    brackets: tuple[SwitchBracket, ...] = ()
    charcoal_node: str | None = None

    def __post_init__(self) -> None:
        """Resolves the allowed switch topology and charcoal mounting."""
        if len(self.nodes) not in (2, 3):
            raise ConfigurationError(f"Switch {self.name!r} needs two or three nodes; got {self.nodes}.")
        if len(self.nodes) == 3 and self.charcoal_node is None:
            raise ConfigurationError(f"Double switch {self.name!r} needs a charcoal_node.")
        if self.charcoal_node is not None and self.charcoal_node not in self.nodes:
            raise ConfigurationError(
                f"Switch {self.name!r} charcoal node {self.charcoal_node!r} must be one of {self.nodes}."
            )

    def calc_gap_conductance(self, temperature: float, pressure: float) -> float:
        """Calculates gas-gap conductance before the series brackets.

        Args:
            temperature: Mean surface temperature in K.
            pressure: Gas pressure in Pa.

        Returns:
            Gas-gap conductance in W/K.
        """
        return calc_gap_conductance(self.gas, temperature, pressure, area=self.area, gap=self.gap)

    def calc_pressure(self, temperature: NDArray[np.float64]) -> float:
        """Calculates the shared equilibrium pressure in Pa from the ordered node temperatures in K."""
        gap_temperature = (temperature[:-1] + temperature[1:]) / 2
        charcoal_temperature = (
            float(np.mean(temperature))
            if self.charcoal_node is None
            else float(temperature[self.nodes.index(self.charcoal_node)])
        )
        amount = self.fill_pressure * self.gas_volume / (R * self.fill_temperature)
        capacity = self.charcoal_mass * self.gas.adsorption_capacity
        free_per_pressure = self.gas_volume / (R * len(gap_temperature)) * np.sum(1 / gap_temperature)
        inverse_affinity = (
            np.exp(-self.gas.adsorption_energy / (R * charcoal_temperature)) / self.gas.adsorption_constant
        )
        linear = free_per_pressure * inverse_affinity + capacity - amount
        discriminant = np.hypot(linear, 2 * np.sqrt(free_per_pressure * amount * inverse_affinity))
        if linear > 0:
            pressure = 2 * amount * inverse_affinity / (linear + discriminant)
        else:
            pressure = (discriminant - linear) / (2 * free_per_pressure)
        saturation = self.gas.saturation_pressure(np.array([np.min(temperature)], dtype=np.float64))[0]
        return float(min(pressure, saturation))

    def calc_energy_balance_terms(self, temperature: NDArray[np.float64]) -> EnergyBalanceTerms:
        """Sums the adjacent-gap heat flows in W into each node at its temperature in K."""
        pressure = self.calc_pressure(temperature)
        heat_flow = np.zeros(len(self.nodes), dtype=np.float64)
        for index in range(len(self.nodes) - 1):
            first, second = temperature[index : index + 2]
            mean = float((first + second) / 2)
            conductance = self.calc_gap_conductance(mean, pressure)
            at_mean = np.array([mean], dtype=np.float64)
            bracket_resistance = sum(
                bracket.length / (bracket.area * bracket.material.thermal_conductivity(at_mean)[0])
                for bracket in self.brackets
            )
            conductance /= 1 + conductance * bracket_resistance
            flow = conductance * (second - first)
            heat_flow[index] += flow
            heat_flow[index + 1] -= flow
        return EnergyBalanceTerms(heat_flow=heat_flow)


def load_gas_gap_switch(config: Mapping[str, Any]) -> GasGapSwitch:
    """Builds a gas-gap switch from its TOML link table.

    Args:
        config: The link's name, two or three ordered nodes, gas name, fill pressure and temperature,
            total gas volume, gap width and area, total charcoal mass, charcoal node, and optional brackets.
            Each bracket gives its material name, area, and length. Charcoal node may be omitted for a
            two-node switch to use its mean temperature. Physical values are SI.

    Returns:
        The switch with its gas and bracket materials loaded.

    Raises:
        KeyError: A required key is missing.
        ConfigurationError: The topology or charcoal mounting is invalid, or a gas or material name is unknown.
    """
    return GasGapSwitch(
        name=config["name"],
        nodes=tuple(config["nodes"]),
        gas=load_gas(config["gas"]),
        fill_pressure=config["fill_pressure"],
        fill_temperature=config["fill_temperature"],
        gas_volume=config["gas_volume"],
        gap=config["gap"],
        area=config["area"],
        charcoal_mass=config["charcoal_mass"],
        charcoal_node=config.get("charcoal_node"),
        brackets=tuple(
            SwitchBracket(material=load_material(bracket["material"]), area=bracket["area"], length=bracket["length"])
            for bracket in config.get("brackets", [])
        ),
    )
