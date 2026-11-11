"""Models cooling power from the pulse tube head."""

from collections.abc import Callable
from dataclasses import dataclass
from typing import Literal

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import brentq

from pulse_tube_cooldown.core.types import EnergyBalanceTerms, Node

__all__ = ["PulseTubeHead"]


@dataclass(frozen=True, slots=True, kw_only=True)
class HeadParameters:
    """Coefficients fitted for a specific head and drive condition."""

    stage_1_warm_amplitude: float  # [W]
    stage_1_warm_exponent: float  # unitless
    stage_2_warm_amplitude: float  # [W]
    stage_2_warm_exponent: float  # unitless
    regenerator_1_exchange_coefficient: float  # [W/(kg K)]
    regenerator_2_exchange_coefficient: float  # [W/(kg K)]
    energy_coupling_coefficient: float  # unitless


@dataclass(frozen=True, slots=True, kw_only=True)
class RegeneratorSegment:
    """One piece of the regenerator matrix, modeled as a Node."""

    node: Node
    regenerator: Literal[1, 2]
    span: tuple[float, float]
    specific_energy: Callable[[float], float]
    inverse_energy: Callable[[float], float]
    target_specific_energy: Callable[[float, float, float], float]


@dataclass(frozen=True, slots=True, kw_only=True)
class LoadMap:
    """A cold heater-power map with its continuation, stand correction, and warm blending.

    Attributes:
        minimum_temperature: Measured temperature minima, in stage order.
        maximum_temperature: Measured temperature maxima, in stage order.
        tail_lengths: Exponential continuation lengths, in stage order.
        blend_width: Width of the warm blend in coordinates normalized by the measured temperature spans.
        heater_power: Maps normalized temperature coordinates of shape (n, 2) to heater powers of the same shape.
        stand_heating: Maps two actual stage temperatures to two stand parasitic heating powers, in stage order.
    """

    minimum_temperature: NDArray[np.float64]
    maximum_temperature: NDArray[np.float64]
    tail_lengths: NDArray[np.float64]
    blend_width: float
    heater_power: Callable[[NDArray[np.float64]], NDArray[np.float64]]
    stand_heating: Callable[[NDArray[np.float64]], NDArray[np.float64]]

    def extend_temperature(
        self,
        temperature: NDArray[np.float64],
    ) -> NDArray[np.float64]:
        """Extends the two stage temperatures with bounded exponential tails outside the measured ranges."""
        above = np.maximum(temperature - self.maximum_temperature, 0.0)
        below = np.maximum(self.minimum_temperature - temperature, 0.0)
        tail = self.tail_lengths

        return (
            np.clip(temperature, self.minimum_temperature, self.maximum_temperature)
            + tail * (-np.expm1(-above / tail))
            - tail * (-np.expm1(-below / tail))
        )

    def base_cooling_power(
        self,
        stage_1_temperature: float,
        stage_2_temperature: float,
    ) -> tuple[float, float]:
        """Calculates base cooling powers from the continued cold map and inferred stand parasitics."""
        temperature = np.array([stage_1_temperature, stage_2_temperature], dtype=np.float64)
        continued = self.extend_temperature(temperature)
        span = self.maximum_temperature - self.minimum_temperature
        coordinates = (continued - self.minimum_temperature) / span

        heater_power = self.heater_power(coordinates[None, :])[0]
        cooling_power = heater_power + self.stand_heating(temperature)
        return float(cooling_power[0]), float(cooling_power[1])

    def warm_blend(
        self,
        stage_1_temperature: float,
        stage_2_temperature: float,
    ) -> float:
        """Calculates the warm-correction blending weight from the distances outside the cold-map bounds."""
        temperature = np.array([stage_1_temperature, stage_2_temperature], dtype=np.float64)
        span = self.maximum_temperature - self.minimum_temperature
        outside = (
            np.maximum(temperature - self.maximum_temperature, 0.0)
            + np.maximum(self.minimum_temperature - temperature, 0.0)
        ) / span

        return float(-np.expm1(-np.sum(outside**2) / self.blend_width**2))


@dataclass(frozen=True, slots=True, kw_only=True)
class ClampResponse:
    """Published clamp-load response and the matrix piece carrying clamp thermal memory."""

    stage_temperature_shifts: tuple[Callable[[float], float], Callable[[float], float]]
    temperature_increment: Callable[[float], float]
    profile_fraction: float
    zero_load_slope: float
    maximum_load: float
    memory_segment: int

    def stage_temperature_shift(self, clamp_load: float) -> tuple[float, float]:
        """Returns the stage temperature shifts caused by clamp loading."""
        stage_1_shift, stage_2_shift = self.stage_temperature_shifts
        return stage_1_shift(clamp_load), stage_2_shift(clamp_load)

    def equilibrium_clamp_temperature(
        self,
        stage_1_temperature: float,
        stage_2_temperature: float,
        *,
        clamp_load: float,
    ) -> float:
        """Returns the head-side clamp temperature at internal equilibrium."""
        return (
            self.profile_fraction * stage_1_temperature
            + (1.0 - self.profile_fraction) * stage_2_temperature
            + self.temperature_increment(clamp_load)
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class PulseTubeHead:
    """Cooling and matrix storage of one pulse-tube head.

    Attributes:
        name: Unique link identifier within the graph.
        stages: Stage-1 and stage-2 node names, in that order.
        segments: The seven matrix pieces, ordered from warm to cold through regenerator 1, then regenerator 2.
        clamp: External clamp-node name, or None for the unloaded stand.
        contact_resistance: Resistance between the external clamp node and the head-side clamp, or None on the stand.
        parameters: Coefficients fitted for the head and drive condition.
        load_map: The cold-map response, stand correction, and warm blending.
        clamp_response: The published clamp-load response and its thermal-memory segment.
        warm_flange_temperature: Prescribed temperature at the warm end of regenerator 1.
        energy_scale: Fixed calibration normalization for excess regenerator-1 energy.
    """

    name: str
    stages: tuple[str, str]
    segments: tuple[RegeneratorSegment, ...]
    clamp: str | None
    contact_resistance: float | None
    parameters: HeadParameters
    load_map: LoadMap
    clamp_response: ClampResponse
    warm_flange_temperature: float
    energy_scale: float

    @property
    def internal_nodes(self) -> tuple[Node, ...]:
        """Matrix nodes in segment order."""
        return tuple(segment.node for segment in self.segments)

    @property
    def nodes(self) -> tuple[str, ...]:
        """Node names in stage order, followed by the optional clamp and the matrix pieces."""
        external = self.stages if self.clamp is None else (*self.stages, self.clamp)
        return (*external, *(segment.node.name for segment in self.segments))

    def targets(
        self,
        stage_1_temperature: float,
        stage_2_temperature: float,
        *,
        clamp_load: float,
    ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Returns target specific energies and caloric temperatures in segment order."""
        specific_energy = np.array(
            [
                segment.target_specific_energy(stage_1_temperature, stage_2_temperature, clamp_load)
                for segment in self.segments
            ],
            dtype=np.float64,
        )
        temperature = np.array(
            [
                segment.inverse_energy(float(energy))
                for segment, energy in zip(self.segments, specific_energy, strict=True)
            ],
            dtype=np.float64,
        )
        return specific_energy, temperature

    def target_temperature(
        self,
        stage_1_temperature: float,
        stage_2_temperature: float,
        *,
        clamp_load: float,
    ) -> NDArray[np.float64]:
        """Calculates the seven target caloric temperatures from the material energies in the target profiles."""
        _, temperature = self.targets(stage_1_temperature, stage_2_temperature, clamp_load=clamp_load)
        return temperature

    def calc_clamp_load(
        self,
        stage_1_temperature: float,
        stage_2_temperature: float,
        *,
        clamp_temperature: float,
        segment_temperature: NDArray[np.float64],
    ) -> float:
        """Solves the mounted clamp's contact equation using the matrix's current caloric temperatures."""
        contact_resistance = self.contact_resistance
        if contact_resistance is None:
            return 0.0

        response = self.clamp_response
        unloaded = response.equilibrium_clamp_temperature(stage_1_temperature, stage_2_temperature, clamp_load=0.0)
        lower = -unloaded / response.zero_load_slope
        upper = response.maximum_load

        def residual(clamp_load: float) -> float:
            target = self.target_temperature(stage_1_temperature, stage_2_temperature, clamp_load=clamp_load)
            equilibrium = response.equilibrium_clamp_temperature(
                stage_1_temperature, stage_2_temperature, clamp_load=clamp_load
            )
            memory = response.memory_segment
            head_temperature = equilibrium * segment_temperature[memory] / target[memory]
            return float(head_temperature + contact_resistance * clamp_load - clamp_temperature)

        while residual(upper) < 0.0:
            upper *= 2.0
        return float(brentq(residual, lower, upper))

    def unloaded_cooling_power(
        self,
        stage_1_temperature: float,
        stage_2_temperature: float,
    ) -> NDArray[np.float64]:
        """Returns the internally equilibrated stage cooling powers at zero clamp load."""
        temperature = np.array([stage_1_temperature, stage_2_temperature], dtype=np.float64)
        base = np.array(self.load_map.base_cooling_power(stage_1_temperature, stage_2_temperature), dtype=np.float64)
        amplitude = np.array(
            [self.parameters.stage_1_warm_amplitude, self.parameters.stage_2_warm_amplitude], dtype=np.float64
        )
        exponent = np.array(
            [self.parameters.stage_1_warm_exponent, self.parameters.stage_2_warm_exponent], dtype=np.float64
        )
        correction = amplitude * (temperature / self.warm_flange_temperature) ** exponent
        return base + self.load_map.warm_blend(stage_1_temperature, stage_2_temperature) * correction

    def calc_energy_balance_terms(self, temperature: NDArray[np.float64]) -> EnergyBalanceTerms:
        """Calculates stage and clamp heat flows and matrix energy rates.

        Args:
            temperature: Temperatures in K in the order of `nodes`.

        Returns:
            Heat flows in W in the order of `nodes`, positive into each node. The matrix nodes carry their material
            heat capacities.
        """
        stage_1_temperature, stage_2_temperature = float(temperature[0]), float(temperature[1])
        if self.clamp is None:
            segment_temperature = temperature[2:]
            clamp_load = 0.0
        else:
            segment_temperature = temperature[3:]
            clamp_load = self.calc_clamp_load(
                stage_1_temperature,
                stage_2_temperature,
                clamp_temperature=float(temperature[2]),
                segment_temperature=segment_temperature,
            )

        target_energy, target_temperature = self.targets(
            stage_1_temperature, stage_2_temperature, clamp_load=clamp_load
        )
        mass = np.array(
            [sum(mass for _, mass in segment.node.composition) for segment in self.segments], dtype=np.float64
        )
        regenerator_1 = np.array([segment.regenerator == 1 for segment in self.segments])
        exchange = np.where(
            regenerator_1,
            self.parameters.regenerator_1_exchange_coefficient,
            self.parameters.regenerator_2_exchange_coefficient,
        )
        energy_rate = mass * exchange * (target_temperature - segment_temperature)
        actual_energy = np.array(
            [
                segment.specific_energy(float(value))
                for segment, value in zip(self.segments, segment_temperature, strict=True)
            ],
            dtype=np.float64,
        )
        excess_energy = np.sum((mass * (actual_energy - target_energy))[regenerator_1])
        multiplier = np.exp(-self.parameters.energy_coupling_coefficient * excess_energy / self.energy_scale)

        shift_1, shift_2 = self.clamp_response.stage_temperature_shift(clamp_load)
        cooling_power = self.unloaded_cooling_power(stage_1_temperature - shift_1, stage_2_temperature - shift_2)
        cooling_power[0] *= multiplier

        midpoint = np.array([sum(segment.span) / 2.0 for segment in self.segments], dtype=np.float64)
        routing = np.array(
            [np.where(regenerator_1, midpoint, 1.0 - midpoint), np.where(regenerator_1, 0.0, midpoint)],
            dtype=np.float64,
        )
        cooling_power += routing @ energy_rate
        external_heat_flow = -cooling_power if self.clamp is None else np.append(-cooling_power, -clamp_load)
        return EnergyBalanceTerms(heat_flow=np.concatenate((external_heat_flow, energy_rate)))
