"""Sealed gas-gap heat switches with charcoal adsorption and series brackets."""

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import NDArray
from scipy.constants import R

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
    """A two-node switch whose sealed gas adsorbs on charcoal as it cools.

    Attributes:
        name: Unique link identifier within a graph.
        nodes: Names of the two nodes joined by the switch.
        gas: Gas properties and Langmuir adsorption constants on charcoal.
        fill_pressure: Absolute gas pressure when filled, in Pa.
        fill_temperature: Gas temperature when filled, in K.
        gas_volume: Sealed gas volume in m^3.
        gap: Gas-gap width in m.
        area: Facing surface area in m^2.
        charcoal_mass: Adsorbent mass in kg.
        brackets: Solid brackets in series with the gas gap.
    """

    name: str
    nodes: tuple[str, str]
    gas: Gas
    fill_pressure: float
    fill_temperature: float
    gas_volume: float
    gap: float
    area: float
    charcoal_mass: float
    brackets: tuple[SwitchBracket, ...] = ()

    def calc_pressure(self, temperature: NDArray[np.float64]) -> float:
        """Calculates equilibrium gas pressure in Pa at the two node temperatures in K."""
        mean = float(np.mean(temperature))
        amount = self.fill_pressure * self.gas_volume / (R * self.fill_temperature)
        capacity = self.charcoal_mass * self.gas.adsorption_capacity
        free_per_pressure = self.gas_volume / (R * mean)
        inverse_affinity = np.exp(-self.gas.adsorption_energy / (R * mean)) / self.gas.adsorption_constant
        linear = free_per_pressure * inverse_affinity + capacity - amount
        discriminant = np.hypot(linear, 2 * np.sqrt(free_per_pressure * amount * inverse_affinity))
        if linear > 0:
            pressure = 2 * amount * inverse_affinity / (linear + discriminant)
        else:
            pressure = (discriminant - linear) / (2 * free_per_pressure)
        saturation = self.gas.saturation_pressure(np.array([np.min(temperature)], dtype=np.float64))[0]
        return float(min(pressure, saturation))

    def calc_energy_balance_terms(self, temperature: NDArray[np.float64]) -> EnergyBalanceTerms:
        """Calculates heat flows in W into the two nodes at their temperatures in K."""
        mean = float(np.mean(temperature))
        conductance = calc_gap_conductance(
            self.gas, mean, self.calc_pressure(temperature), area=self.area, gap=self.gap
        )
        at_mean = np.array([mean], dtype=np.float64)
        bracket_resistance = sum(
            bracket.length / (bracket.area * bracket.material.thermal_conductivity(at_mean)[0])
            for bracket in self.brackets
        )
        conductance /= 1 + conductance * bracket_resistance
        flow = conductance * (temperature[1] - temperature[0])
        return EnergyBalanceTerms(heat_flow=np.array([flow, -flow], dtype=np.float64))


def load_gas_gap_switch(config: Mapping[str, Any]) -> GasGapSwitch:
    """Builds a gas-gap switch from its TOML link table.

    Args:
        config: The link's name, two nodes, gas name, fill pressure and temperature,
            gas volume, gap, area, charcoal mass, and optional brackets. Each bracket
            gives its material name, area, and length. Physical values are SI.

    Returns:
        The switch with its gas and bracket materials loaded.

    Raises:
        KeyError: A required key is missing.
        ValueError: `nodes` does not hold exactly two entries.
        ConfigurationError: A gas or material name is unknown.
    """
    first, second = config["nodes"]
    return GasGapSwitch(
        name=config["name"],
        nodes=(first, second),
        gas=load_gas(config["gas"]),
        fill_pressure=config["fill_pressure"],
        fill_temperature=config["fill_temperature"],
        gas_volume=config["gas_volume"],
        gap=config["gap"],
        area=config["area"],
        charcoal_mass=config["charcoal_mass"],
        brackets=tuple(
            SwitchBracket(material=load_material(bracket["material"]), area=bracket["area"], length=bracket["length"])
            for bracket in config.get("brackets", [])
        ),
    )
