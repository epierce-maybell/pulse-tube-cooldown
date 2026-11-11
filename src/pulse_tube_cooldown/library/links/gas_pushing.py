"""Cycle-averaged heat transfer from docs/adrs/02-gas-pushing-routine.md."""

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import NDArray
from scipy.constants import R

from pulse_tube_cooldown.core.errors import ConfigurationError
from pulse_tube_cooldown.core.types import EnergyBalanceTerms, Link

__all__ = ["GasPushing", "GasVolume", "load_gas_pushing"]


@dataclass(frozen=True, slots=True, kw_only=True)
class GasVolume:
    """A gas volume in m^3 sharing heat equally between its one or two thermal nodes."""

    nodes: tuple[str, ...]
    volume: float


@dataclass(frozen=True, slots=True, kw_only=True)
class GasPushing(Link):
    """Heat transfer through an ordered gas path with wave and buffer-pressure limits.

    Attributes:
        name: Link name.
        nodes: Thermal nodes, warm manifold first.
        volumes: Gas volumes ordered from the manifold inward.
        warm_volume: Open manifold volume, m^3.
        buffer_volume: Buffer volume at the manifold temperature, m^3.
        pump_speed: Pumping speed at the manifold temperature, m^3/s.
        flow_conductance: Buffer-to-path molar conductance, mol/(s Pa).
        fill_pressure: Absolute still-side fill pressure, Pa.
        wave_strength: Buffer-minus-still pressure that opens the wave path, Pa.
        refill_pressure: Drop below fill pressure that triggers refill, Pa.
        buffer_pressure_cap: Absolute buffer pressure that opens the wave path, Pa.
        depressurization_time: Time the wave path stays open, s.
    """

    name: str
    nodes: tuple[str, ...]
    volumes: tuple[GasVolume, ...]
    warm_volume: float
    buffer_volume: float
    pump_speed: float
    flow_conductance: float
    fill_pressure: float
    wave_strength: float
    refill_pressure: float
    buffer_pressure_cap: float
    depressurization_time: float

    def calc_energy_balance_terms(self, temperature: NDArray[np.float64]) -> EnergyBalanceTerms:
        """Returns heat flow into each node in W for node temperatures in K, both in `nodes` order."""
        volume_node_indices = [tuple(self.nodes.index(node) for node in volume.nodes) for volume in self.volumes]
        gas_temperature = np.array([_log_mean(temperature[list(indices)]) for indices in volume_node_indices])
        volume_molar_capacity = np.array([volume.volume for volume in self.volumes]) / (R * gas_temperature)
        downstream_molar_capacity = np.cumsum(volume_molar_capacity[::-1])[::-1]
        warm_temperature = float(temperature[0])
        path_molar_capacity = self.warm_volume / (R * warm_temperature) + downstream_molar_capacity[0]
        buffer_molar_capacity = self.buffer_volume / (R * warm_temperature)
        cycle_pressure_rate = self.calc_pressure_rate(path_molar_capacity, buffer_molar_capacity, warm_temperature)

        upstream_temperature = np.concatenate(([warm_temperature], gas_temperature[:-1]))
        molar_heat_capacity = 2.5 * R
        boundary_heat_flow = (
            molar_heat_capacity
            * cycle_pressure_rate
            * downstream_molar_capacity
            * (upstream_temperature - gas_temperature)
        )
        volume_heat_flow = boundary_heat_flow - np.concatenate((boundary_heat_flow[1:], [0.0]))
        node_heat_flow = np.zeros(len(self.nodes), dtype=np.float64)
        node_heat_flow[0] = -boundary_heat_flow[0]
        for indices, heat_flow in zip(volume_node_indices, volume_heat_flow, strict=True):
            np.add.at(node_heat_flow, list(indices), heat_flow / len(indices))
        return EnergyBalanceTerms(heat_flow=node_heat_flow)

    def calc_pressure_rate(
        self, path_molar_capacity: float, buffer_molar_capacity: float, warm_temperature: float
    ) -> float:
        """Returns cycle frequency times still-side pressure excursion, in Pa/s.

        Args:
            path_molar_capacity: Still-side gas inventory per unit pressure, mol/Pa.
            buffer_molar_capacity: Buffer gas inventory per unit pressure, mol/Pa.
            warm_temperature: Manifold and buffer temperature, K.
        """
        if self.pump_speed == 0.0 or self.flow_conductance == 0.0 or self.depressurization_time == 0.0:
            return 0.0
        equalization_time = (
            path_molar_capacity
            * buffer_molar_capacity
            / (self.flow_conductance * (path_molar_capacity + buffer_molar_capacity))
        )
        excursion_fraction = (
            buffer_molar_capacity
            / (path_molar_capacity + buffer_molar_capacity)
            * -np.expm1(-self.depressurization_time / equalization_time)
        )
        pressure_excursion = self.calc_pressure_excursion(excursion_fraction)
        if pressure_excursion == 0.0:
            return 0.0
        maximum_path_pressure = self.fill_pressure - max(self.refill_pressure - pressure_excursion, 0.0) / 2.0
        pumping_time = (
            -path_molar_capacity
            * R
            * warm_temperature
            / self.pump_speed
            * np.log1p(-pressure_excursion / maximum_path_pressure)
        )
        return pressure_excursion / (pumping_time + self.depressurization_time)

    def calc_pressure_excursion(self, excursion_fraction: float) -> float:
        """Returns the still-side pressure excursion allowed by both opening limits, in Pa.

        Args:
            excursion_fraction: Still-side pressure rise per unit buffer-minus-still
                pressure difference during the valve-open interval.
        """
        cap_excursion = (
            excursion_fraction
            * (self.buffer_pressure_cap - self.fill_pressure + self.refill_pressure / 2.0)
            / (1.0 - excursion_fraction / 2.0)
        )
        if cap_excursion > self.refill_pressure:
            cap_excursion = (
                excursion_fraction * (self.buffer_pressure_cap - self.fill_pressure) / (1.0 - excursion_fraction)
            )
        return max(0.0, min(excursion_fraction * self.wave_strength, cap_excursion))


def _log_mean(temperature: NDArray[np.float64]) -> float:
    """Returns the log-mean temperature of one or two thermal nodes, in K."""
    first, last = float(temperature[0]), float(temperature[-1])
    if first == last:
        return first
    difference = last - first
    return difference / np.log1p(difference / first)


def load_gas_pushing(config: Mapping[str, Any]) -> GasPushing:
    """Builds a gas-pushing link from its configuration table.

    Args:
        config: Fields of `GasPushing`, with each `volumes` entry holding `nodes` and `volume`.

    Raises:
        KeyError: A required field is missing.
        ConfigurationError: A volume touches other than one or two of the link's nodes, or no volumes are given.
    """
    nodes = tuple(config["nodes"])
    volumes = tuple(GasVolume(nodes=tuple(volume["nodes"]), volume=volume["volume"]) for volume in config["volumes"])
    if not volumes:
        raise ConfigurationError(f"Link {config['name']!r} needs at least one gas volume.")
    for volume in volumes:
        if len(volume.nodes) not in (1, 2) or not set(volume.nodes).issubset(nodes):
            raise ConfigurationError(
                f"Link {config['name']!r} volume {volume.nodes!r} must touch one or two nodes from {nodes!r}."
            )
    return GasPushing(
        name=config["name"],
        nodes=nodes,
        volumes=volumes,
        warm_volume=config["warm_volume"],
        buffer_volume=config["buffer_volume"],
        pump_speed=config["pump_speed"],
        flow_conductance=config["flow_conductance"],
        fill_pressure=config["fill_pressure"],
        wave_strength=config["wave_strength"],
        refill_pressure=config["refill_pressure"],
        buffer_pressure_cap=config["buffer_pressure_cap"],
        depressurization_time=config["depressurization_time"],
    )
