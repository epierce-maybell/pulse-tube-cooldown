"""A pulse-tube head whose cooling follows from cycle-averaged energy flows."""

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from pulse_tube_cooldown.core.types import EnergyBalanceTerms, Material, Node
from pulse_tube_cooldown.library.gases import helium
from pulse_tube_cooldown.library.links.pulse_tube.energy_flow.manifold import Manifold

__all__ = ["PulseTubeHead"]


@dataclass(frozen=True, slots=True, kw_only=True)
class PulseTubeHead:
    """Models a pulse tube head as a reduced dynamical system.

    Attributes:
        name: Unique link identifier within a graph.
        stages: Names of the stage-1 and stage-2 nodes.
        regenerator_plate: Name of the regenerator-plate node, or None when nothing is clamped.
        warm_flange_temperature: Temperature of the warm flange.
        clamp_contact_conductance: Thermal conductance between the regenerator plate and the clamp, zero if unclamped.
        effective_impedance: Acoustic impedance that sets the acoustic power entering regenerator 1, p_a^2 / (2 Z).
        stage_1_split: Fraction of the acoustic power reaching stage 1 that enters regenerator 2.
        valve_frequency: Frequency of the rotary valve's pressure cycle.
        nusselt_coeff: Coefficient of the Peclet-number term in the matrix's Nusselt number.
        nusselt_peclet_exp: Exponent of the Peclet number in the matrix's Nusselt number.
        nusselt_porosity_exp: Exponent of the porosity factor in the matrix's Nusselt number.
        axial_conductivity_coeff: Coefficient of the Peclet-number term in the gas's axial conductivity ratio.
        axial_conductivity_peclet_exp: Exponent of the Peclet number in the gas's axial conductivity ratio.
        axial_conductivity_porosity_exp: Exponent of the porosity factor in the gas's axial conductivity ratio.
        manifold: The compressors and gas supply driving the head.
        regenerator_1: Nodes of the regenerator-1 sections, from the warm flange to stage 1.
        regenerator_2_upper: Nodes of the regenerator-2 sections, from stage 1 to the clamp.
        regenerator_2_lower: Nodes of the regenerator-2 sections, from the clamp to stage 2.
        section_length: Length of each section along its regenerator, from warm to cold.
        section_flow_area: Cross-section of the regenerator bore at each section.
        section_porosity: Void fraction of each section's matrix.
        section_hydraulic_diameter: Hydraulic diameter of each section's matrix pores.
        section_wall_area: Cross-section of the regenerator wall at each section.
        pulse_tube_wall_area: Wall cross-section of pulse tubes 1 and 2.
        pulse_tube_length: Length of pulse tubes 1 and 2, from their stages to the warm flange.
        pulse_tube_bore_area: Bore cross-section of pulse tubes 1 and 2.
        pulse_tube_segments: Number of equal segments of each pulse tube's linear temperature profile.
        heat_exchanger_void_volume: Gas volume of the stage-1 and stage-2 heat exchangers.
        dead_volume: Gas volume between the rotary valve and regenerator 1, at the warm-flange temperature.
        reservoir_volume: Gas volume of the head's reservoirs, at the warm-flange temperature.
        wall_material: Material of the regenerator and pulse-tube walls.
    """

    name: str
    stages: tuple[str, str]
    regenerator_plate: str | None
    warm_flange_temperature: float
    clamp_contact_conductance: float
    effective_impedance: float
    stage_1_split: float
    valve_frequency: float

    nusselt_coeff: float
    nusselt_peclet_exp: float
    nusselt_porosity_exp: float
    axial_conductivity_coeff: float
    axial_conductivity_peclet_exp: float
    axial_conductivity_porosity_exp: float

    manifold: Manifold

    regenerator_1: tuple[Node, ...]
    regenerator_2_upper: tuple[Node, ...]
    regenerator_2_lower: tuple[Node, ...]
    section_length: NDArray[np.float64]
    section_flow_area: NDArray[np.float64]
    section_porosity: NDArray[np.float64]
    section_hydraulic_diameter: NDArray[np.float64]
    section_wall_area: NDArray[np.float64]
    pulse_tube_wall_area: NDArray[np.float64]
    pulse_tube_length: NDArray[np.float64]
    pulse_tube_bore_area: NDArray[np.float64]
    pulse_tube_segments: int
    heat_exchanger_void_volume: NDArray[np.float64]
    dead_volume: float
    reservoir_volume: float
    wall_material: Material

    @property
    def nodes(self) -> tuple[str, ...]:
        """Return nodes attached to this link."""
        plate = () if self.regenerator_plate is None else (self.regenerator_plate,)
        return (*self.stages, *plate, *(node.name for node in self.internal_nodes))

    @property
    def internal_nodes(self) -> tuple[Node, ...]:
        """Return internal nodes that this link adds to the graph, from warm to cold."""
        return (*self.regenerator_1, *self.regenerator_2_upper, *self.regenerator_2_lower)

    def calc_energy_balance_terms(self, temperature: NDArray[np.float64]) -> EnergyBalanceTerms:
        """Return heat flows at every node connected to this link."""
        clamped = self.regenerator_plate is not None
        stage_1_temperature, stage_2_temperature = temperature[:2]
        regenerator_plate_temperature = temperature[2] if clamped else 0.0
        section_temperatures = temperature[3:] if clamped else temperature[2:]
        upper_start = len(self.regenerator_1)
        lower_start = upper_start + len(self.regenerator_2_upper)

        mean_pressure, pressure_amplitude = self.drive(stage_1_temperature, section_temperatures, stage_2_temperature)

        gas = helium.properties(section_temperatures, mean_pressure)
        matrix_heat_capacity = np.array(
            [node.calc_heat_capacity(t) for node, t in zip(self.internal_nodes, section_temperatures, strict=True)]
        )
        gas_heat_capacity = (
            gas.density
            * gas.isobaric_specific_heat
            * self.section_porosity
            * self.section_flow_area
            * self.section_length
        )
        matrix_fraction = matrix_heat_capacity / (matrix_heat_capacity + gas_heat_capacity)
        transmission = matrix_fraction * gas.isobaric_expansion_coefficient * section_temperatures

        section_powers, pt1_power, pt2_power = self.acoustic_powers(
            stage_1_temperature, section_temperatures, stage_2_temperature, transmission, pressure_amplitude
        )
        half_resistance, half_offset = self.half_section_transport(
            section_temperatures,
            section_powers,
            pressure_amplitude,
            gas=gas,
            matrix_fraction=matrix_fraction,
            transmission=transmission,
        )

        upper_last, lower_first = lower_start - 1, lower_start
        head_clamp_temperature = (
            (section_temperatures[upper_last] + half_offset[upper_last]) / half_resistance[upper_last]
            + (section_temperatures[lower_first] - half_offset[lower_first]) / half_resistance[lower_first]
            + self.clamp_contact_conductance * regenerator_plate_temperature
        ) / (1.0 / half_resistance[upper_last] + 1.0 / half_resistance[lower_first] + self.clamp_contact_conductance)

        chains = (
            (self.warm_flange_temperature, slice(0, upper_start), stage_1_temperature),
            (stage_1_temperature, slice(upper_start, lower_start), head_clamp_temperature),
            (head_clamp_temperature, slice(lower_start, None), stage_2_temperature),
        )
        r1_energy_flows, r2_upper_energy_flows, r2_lower_energy_flows = (
            self.chain_energy_flows(warm, section_temperatures[part], cold, half_resistance[part], half_offset[part])
            for warm, part, cold in chains
        )
        pt1_energy_flow, pt2_energy_flow = self.pulse_tube_energy_flows(
            stage_1_temperature, stage_2_temperature, pt1_power, pt2_power
        )

        stage_1_cooling = pt1_energy_flow + r2_upper_energy_flows[0] - r1_energy_flows[-1]
        clamp_cooling = r2_lower_energy_flows[0] - r2_upper_energy_flows[-1]
        stage_2_cooling = pt2_energy_flow - r2_lower_energy_flows[-1]
        section_heat_flows = np.concatenate(
            [flows[:-1] - flows[1:] for flows in (r1_energy_flows, r2_upper_energy_flows, r2_lower_energy_flows)]
        )

        return EnergyBalanceTerms(
            heat_flow=np.array(
                [-stage_1_cooling, -stage_2_cooling, *([-clamp_cooling] if clamped else []), *section_heat_flows]
            )
        )

    def drive(
        self, stage_1_temperature: float, section_temperatures: NDArray[np.float64], stage_2_temperature: float
    ) -> tuple[float, float]:
        """Return the mean pressure and pressure amplitude that the manifold sustains with this head's gas."""
        warm = self.warm_flange_temperature
        midpoints = (np.arange(self.pulse_tube_segments) + 0.5) / self.pulse_tube_segments
        pulse_tube_temperatures = np.concatenate(
            [stage + (warm - stage) * midpoints for stage in (stage_1_temperature, stage_2_temperature)]
        )
        gas = helium.properties(
            np.concatenate(
                [section_temperatures, [stage_1_temperature, stage_2_temperature, warm], pulse_tube_temperatures]
            ),
            self.manifold.fill_pressure,
        )
        volume = np.concatenate(
            [
                self.section_porosity * self.section_flow_area * self.section_length,
                self.heat_exchanger_void_volume,
                [self.dead_volume],
                np.repeat(
                    self.pulse_tube_bore_area * self.pulse_tube_length / self.pulse_tube_segments,
                    self.pulse_tube_segments,
                ),
            ]
        )
        warm_index = len(self.section_length) + 2
        pulse_tubes = slice(warm_index + 1, None)
        density_slope = gas.density * gas.isothermal_compressibility
        adiabatic_slope = density_slope.copy()
        adiabatic_slope[pulse_tubes] /= gas.heat_capacity_ratio[pulse_tubes]

        return self.manifold.pressures(
            warm_temperature=warm,
            head_volume=volume.sum() + self.reservoir_volume,
            head_mass=volume @ gas.density + self.reservoir_volume * gas.density[warm_index],
            head_mass_slope=volume @ density_slope + self.reservoir_volume * density_slope[warm_index],
            head_admittance=np.hypot(
                gas.density[warm_index] / self.effective_impedance,
                2.0 * np.pi * self.valve_frequency * (volume @ adiabatic_slope),
            ),
        )

    def nusselt_number(self, peclet: NDArray[np.float64], porosity: NDArray[np.float64]) -> NDArray[np.float64]:
        """Return the Nusselt number for heat exchange between the gas and the matrix."""
        return (1.0 + self.nusselt_coeff * peclet**self.nusselt_peclet_exp) * porosity**self.nusselt_porosity_exp

    def axial_conductivity_ratio(
        self, peclet: NDArray[np.float64], porosity: NDArray[np.float64]
    ) -> NDArray[np.float64]:
        """Return the ratio of the gas's axial conductivity in the matrix to its molecular conductivity."""
        return 1.0 + (
            self.axial_conductivity_coeff
            * peclet**self.axial_conductivity_peclet_exp
            * porosity**self.axial_conductivity_porosity_exp
        )

    def acoustic_powers(
        self,
        stage_1_temperature: float,
        section_temperatures: NDArray[np.float64],
        stage_2_temperature: float,
        transmission: NDArray[np.float64],
        pressure_amplitude: float,
    ) -> tuple[NDArray[np.float64], float, float]:
        """Return the acoustic power at each section, and entering pulse tubes 1 and 2."""
        upper_start = len(self.regenerator_1)

        r1_powers = self.chain_acoustic_powers(
            pressure_amplitude**2 / (2.0 * self.effective_impedance),
            self.warm_flange_temperature,
            section_temperatures[:upper_start],
            stage_1_temperature,
            transmission[:upper_start],
        )
        stage_1_power = r1_powers[-1]
        r2_powers = self.chain_acoustic_powers(
            self.stage_1_split * stage_1_power,
            stage_1_temperature,
            section_temperatures[upper_start:],
            stage_2_temperature,
            transmission[upper_start:],
        )

        return (
            np.concatenate([r1_powers[1:-1], r2_powers[1:-1]]),
            (1.0 - self.stage_1_split) * stage_1_power,
            r2_powers[-1],
        )

    def chain_acoustic_powers(
        self,
        warm_end_power: float,
        warm_end_temperature: float,
        section_temperatures: NDArray[np.float64],
        cold_end_temperature: float,
        transmission: NDArray[np.float64],
    ) -> NDArray[np.float64]:
        """Return the acoustic power at each point of a chain: its warm end, each section, and its cold end."""
        points = np.concatenate([[warm_end_temperature], section_temperatures, [cold_end_temperature]])
        interval_transmission = 0.5 * (
            np.append(transmission, transmission[-1]) + np.insert(transmission, 0, transmission[0])
        )
        log_steps = interval_transmission * np.diff(np.log(points))

        return warm_end_power * np.exp(np.concatenate([[0.0], np.cumsum(log_steps)]))

    def half_section_transport(
        self,
        section_temperatures: NDArray[np.float64],
        section_powers: NDArray[np.float64],
        pressure_amplitude: float,
        *,
        gas: helium.HeliumProperties,
        matrix_fraction: NDArray[np.float64],
        transmission: NDArray[np.float64],
    ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Return each section's half resistance (length / 2) / K and half offset (length / 2) * Phi / K."""
        mean_square_mass_flow = 2.0 * (gas.density * section_powers / pressure_amplitude) ** 2
        void_area = self.section_porosity * self.section_flow_area
        diameter = self.section_hydraulic_diameter
        peclet = (
            np.sqrt(mean_square_mass_flow)
            / void_area
            * diameter
            * gas.isobaric_specific_heat
            / gas.thermal_conductivity
        )
        exchange_conductance = (
            self.nusselt_number(peclet, self.section_porosity)
            * gas.thermal_conductivity
            * 4.0
            * void_area
            / diameter**2
        )

        axial_conductance = (
            matrix_fraction**2 * gas.isobaric_specific_heat**2 * mean_square_mass_flow / exchange_conductance
            + self.axial_conductivity_ratio(peclet, self.section_porosity) * gas.thermal_conductivity * void_area
            + self.wall_material.thermal_conductivity(section_temperatures) * self.section_wall_area
        )

        half_resistance = 0.5 * self.section_length / axial_conductance
        return half_resistance, half_resistance * (1.0 - transmission) * section_powers

    def pulse_tube_energy_flows(
        self, stage_1_temperature: float, stage_2_temperature: float, pt1_power: float, pt2_power: float
    ) -> tuple[float, float]:
        """Return the energy flows up pulse tubes 1 and 2, away from their stages."""
        warm = self.wall_material.conductivity_integral(np.array([self.warm_flange_temperature]))[0]
        cold = self.wall_material.conductivity_integral(np.array([stage_1_temperature, stage_2_temperature]))
        conduction = self.pulse_tube_wall_area / self.pulse_tube_length * (warm - cold)
        return pt1_power - conduction[0], pt2_power - conduction[1]

    def chain_energy_flows(
        self,
        warm_end_temperature: float,
        section_temperatures: NDArray[np.float64],
        cold_end_temperature: float,
        half_resistance: NDArray[np.float64],
        half_offset: NDArray[np.float64],
    ) -> NDArray[np.float64]:
        """Return the energy flows, positive toward the cold end, between neighboring points of a chain."""
        points = np.concatenate([[warm_end_temperature], section_temperatures, [cold_end_temperature]])
        interval_resistance = np.append(half_resistance, 0.0) + np.insert(half_resistance, 0, 0.0)
        interval_offset = np.append(half_offset, 0.0) + np.insert(half_offset, 0, 0.0)

        return (points[:-1] - points[1:] + interval_offset) / interval_resistance
