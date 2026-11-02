"""Network equations of one pulse-tube head, from `docs/adrs/coupled-capacity-maps.md`.

A head is a chain of the warm flange, the stage-1 regenerator column's nodes from its warm end, stage 1, the stage-2
column's nodes, and stage 2. Temperatures and returned heat flows follow that order. A boundary joins two neighbors in
the chain, so a head with n_1 and n_2 column nodes has n_1 + n_2 + 2 boundaries, ordered from the warm flange. Pore
gas takes its node's temperature, and each pulse tube's gas a linear profile from the warm flange to its stage.

Each boundary carries the regenerator loss toward the stage end as a conductance,

    G = 2 * a * (m * c_p)^2 / (f * (C_warm + C_cold)),

for the column's packing constant a, the mass flow m across the boundary, helium's isobaric specific heat c_p at the
mean pressure and the boundary's mean temperature, the valve frequency f, and the heat capacities C of the matrix on
either side, zero at a plate. Each stage-2 boundary also carries the real-gas enthalpy flow
(1 - beta * T) * rho(T_2) / rho(T) * H_2 toward stage 2, for helium's isobaric expansion coefficient beta and density
rho at the mean pressure. T_2 is the stage-2 temperature, and the boundary temperature T is the plate's temperature
at a plate, otherwise the mean of its two nodes. Each stage sinks its enthalpy flow H, which lands on the warm flange.
"""

from dataclasses import dataclass
from itertools import accumulate

import numpy as np
from numpy.typing import NDArray
from scipy.special import roots_legendre

from pulse_tube_cooldown.core.types import Node
from pulse_tube_cooldown.library.gases import helium

__all__ = [
    "PulseTubeHead",
    "RegeneratorLayer",
    "calc_capacity_ratio",
    "calc_enthalpy_flow",
    "calc_gas_mass",
    "calc_heat_flow",
    "calc_mass_flow",
    "calc_utilization",
]


@dataclass(frozen=True, slots=True, kw_only=True)
class RegeneratorLayer:
    """One matrix node and its gas volume and span along a column.

    Attributes:
        node: Node whose composition holds the matrix material and mass in kg.
        pore_volume: Gas volume in m^3 at the matrix temperature.
        length_fraction: Fraction of the column's length occupied by this layer.
    """

    node: Node
    pore_volume: float
    length_fraction: float


@dataclass(frozen=True, slots=True, kw_only=True)
class PulseTubeHead:
    """One unit's model geometry, matrix nodes, and per-unit parameters.

    Attributes:
        name: Unit's name in the configuration.
        stages: Stage-1 and stage-2 node names.
        columns: Stage-1 and stage-2 layers, each ordered from its warm end.
        pulse_tube_volumes: Stage-1 and stage-2 pulse-tube gas volumes in m^3.
        clamp_fraction: Clamp position along the stage-2 column, measured from stage 1.
        packing_constants: Column packing constants a_1 and a_2.
        conductances: Stage valve conductances C_1 and C_2 in W/Pa^2.
        frequency: Valve frequency in Hz.
    """

    name: str
    stages: tuple[str, str]
    columns: tuple[tuple[RegeneratorLayer, ...], tuple[RegeneratorLayer, ...]]
    pulse_tube_volumes: tuple[float, float]
    clamp_fraction: float
    packing_constants: tuple[float, float]
    conductances: tuple[float, float]
    frequency: float

    @property
    def matrix_nodes(self) -> tuple[Node, ...]:
        """Matrix nodes in column order, stage 1 followed by stage 2."""
        return tuple(layer.node for column in self.columns for layer in column)

    @property
    def clamp_node(self) -> str:
        """Name of the stage-2 layer containing the clamp; a boundary belongs to the colder layer."""
        ends = accumulate(layer.length_fraction for layer in self.columns[1])
        for layer, end in zip(self.columns[1], ends, strict=True):
            if self.clamp_fraction < end:
                return layer.node.name
        return self.columns[1][-1].node.name

    @property
    def pore_volume(self) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Pore gas volumes in m^3 at each node of the two columns."""
        first, second = (np.array([layer.pore_volume for layer in column]) for column in self.columns)
        return first, second

    def chain_nodes(self, warm_flange: str) -> tuple[str, ...]:
        """Node names from the warm flange through both columns and stages."""
        return (
            warm_flange,
            *(layer.node.name for layer in self.columns[0]),
            self.stages[0],
            *(layer.node.name for layer in self.columns[1]),
            self.stages[1],
        )

    def calc_heat_capacity(self, temperature: NDArray[np.float64]) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Matrix heat capacities in J/K, from temperatures in K along the head's chain."""
        stage_1 = len(self.columns[0]) + 1
        column_temperatures = (temperature[1:stage_1], temperature[stage_1 + 1 : -1])
        first, second = (
            np.array(
                [layer.node.calc_heat_capacity(float(t)) for layer, t in zip(column, temperatures, strict=True)],
                dtype=np.float64,
            )
            for column, temperatures in zip(self.columns, column_temperatures, strict=True)
        )
        return first, second


def calc_mass_flow(
    temperature: NDArray[np.float64],
    pore_volume: tuple[NDArray[np.float64], NDArray[np.float64]],
    pulse_tube_volume: tuple[float, float],
    supply_pressure: float,
    return_pressure: float,
    frequency: float,
) -> NDArray[np.float64]:
    """Calculates the mass flow across each boundary of a head.

    Args:
        temperature: Temperatures in K along the head's chain.
        pore_volume: Pore gas volume in m^3 of each node in the stage-1 and stage-2 columns.
        pulse_tube_volume: Gas volumes in m^3 of the stage-1 and stage-2 pulse tubes.
        supply_pressure: Supply pressure in Pa.
        return_pressure: Return pressure in Pa.
        frequency: Valve frequency in Hz.

    Returns:
        Mass flow in kg/s across each boundary, from the warm flange.
    """
    mass_swing = _gas_mass_along_chain(temperature, pore_volume, pulse_tube_volume, supply_pressure)
    mass_swing -= _gas_mass_along_chain(temperature, pore_volume, pulse_tube_volume, return_pressure)
    return frequency * np.cumsum(mass_swing[::-1])[::-1][1:]


def calc_gas_mass(
    temperature: NDArray[np.float64],
    pore_volume: tuple[NDArray[np.float64], NDArray[np.float64]],
    pulse_tube_volume: tuple[float, float],
    pressure: float,
) -> float:
    """Calculates the helium mass in a head's pores and pulse tubes at one pressure.

    Args:
        temperature: Temperatures in K along the head's chain.
        pore_volume: Pore gas volume in m^3 of each node in the stage-1 and stage-2 columns.
        pulse_tube_volume: Gas volumes in m^3 of the stage-1 and stage-2 pulse tubes.
        pressure: Pressure in Pa.

    Returns:
        Helium mass in kg.
    """
    return float(_gas_mass_along_chain(temperature, pore_volume, pulse_tube_volume, pressure).sum())


def calc_enthalpy_flow(conductance: float, pressure_swing: float, waveform_factor: float) -> float:
    """Calculates the cycle-averaged enthalpy flow away from a stage, H = C * (kappa * dP / 2)^2 / 2.

    Args:
        conductance: The stage valve's effective conductance C in W/Pa^2.
        pressure_swing: Supply minus return pressure, dP, in Pa.
        waveform_factor: The drive's waveform factor kappa.

    Returns:
        Enthalpy flow in W.
    """
    return 0.5 * conductance * (0.5 * waveform_factor * pressure_swing) ** 2


def calc_heat_flow(
    temperature: NDArray[np.float64],
    heat_capacity: tuple[NDArray[np.float64], NDArray[np.float64]],
    mass_flow: NDArray[np.float64],
    enthalpy_flow: tuple[float, float],
    packing_constant: tuple[float, float],
    mean_pressure: float,
    frequency: float,
) -> NDArray[np.float64]:
    """Calculates the heat flow into each element of a head.

    Args:
        temperature: Temperatures in K along the head's chain.
        heat_capacity: Matrix heat capacity in J/K of each node in the stage-1 and stage-2 columns, at its temperature.
        mass_flow: Mass flow in kg/s across each boundary, from `calc_mass_flow`.
        enthalpy_flow: Enthalpy flows H_1 and H_2 in W away from stages 1 and 2.
        packing_constant: Packing constants a_1 and a_2 of the stage-1 and stage-2 columns.
        mean_pressure: Mean of supply and return pressures in Pa.
        frequency: Valve frequency in Hz.

    Returns:
        Heat flow in W into each element of the chain, summing to zero.
    """
    stage_1 = _stage_1(heat_capacity[0])
    capacity = _along_chain(*heat_capacity)
    warm_side, cold_side = temperature[:-1], temperature[1:]
    mean_temperature = (warm_side + cold_side) / 2
    packing = np.repeat(packing_constant, [stage_1, len(mass_flow) - stage_1])
    specific_heat = helium.isobaric_specific_heat(mean_temperature, mean_pressure)
    conductance = 2 * packing * (mass_flow * specific_heat) ** 2 / (frequency * (capacity[:-1] + capacity[1:]))
    flow_toward_stage = conductance * (warm_side - cold_side)

    boundary_temperature = mean_temperature[stage_1:].copy()
    boundary_temperature[[0, -1]] = temperature[[stage_1, -1]]
    expansion = helium.isobaric_expansion_coefficient(boundary_temperature, mean_pressure)
    boundary_density = helium.density(boundary_temperature, mean_pressure)
    flow_toward_stage[stage_1:] += (
        (1 - expansion * boundary_temperature) * boundary_density[-1] / boundary_density * enthalpy_flow[1]
    )

    heat_flow = -np.diff(flow_toward_stage, prepend=0.0, append=0.0)
    heat_flow[stage_1] -= enthalpy_flow[0]
    heat_flow[-1] -= enthalpy_flow[1]
    heat_flow[0] += sum(enthalpy_flow)
    return heat_flow


def calc_utilization(
    temperature: NDArray[np.float64],
    heat_capacity: tuple[NDArray[np.float64], NDArray[np.float64]],
    mass_flow: NDArray[np.float64],
    mean_pressure: float,
    frequency: float,
) -> tuple[float, float]:
    """Calculates each column's utilization, m * c_p / (f * C), with m and c_p at its stage.

    Args:
        temperature: Temperatures in K along the head's chain.
        heat_capacity: Matrix heat capacity in J/K of each node in the stage-1 and stage-2 columns, at its temperature.
        mass_flow: Mass flow in kg/s across each boundary, from `calc_mass_flow`.
        mean_pressure: Mean of supply and return pressures in Pa.
        frequency: Valve frequency in Hz.

    Returns:
        Utilization of the stage-1 and stage-2 columns.
    """
    stage_1 = _stage_1(heat_capacity[0])
    specific_heat = helium.isobaric_specific_heat(temperature[[stage_1, -1]], mean_pressure)
    capacity = np.array([column.sum() for column in heat_capacity])
    first, second = mass_flow[[stage_1 - 1, -1]] * specific_heat / (frequency * capacity)
    return float(first), float(second)


def calc_capacity_ratio(
    temperature: NDArray[np.float64],
    heat_capacity: tuple[NDArray[np.float64], NDArray[np.float64]],
    pore_volume: tuple[NDArray[np.float64], NDArray[np.float64]],
    mean_pressure: float,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Calculates each column node's pore-gas to matrix heat-capacity ratio, V * rho * c_p / C.

    Args:
        temperature: Temperatures in K along the head's chain.
        heat_capacity: Matrix heat capacity in J/K of each node in the stage-1 and stage-2 columns, at its temperature.
        pore_volume: Pore gas volume in m^3 of each node in the stage-1 and stage-2 columns.
        mean_pressure: Mean of supply and return pressures in Pa.

    Returns:
        Ratio at each node of the stage-1 and stage-2 columns.
    """
    stage_1 = _stage_1(heat_capacity[0])
    gas_capacity = (
        _along_chain(*pore_volume)
        * helium.density(temperature, mean_pressure)
        * helium.isobaric_specific_heat(temperature, mean_pressure)
    )
    return gas_capacity[1:stage_1] / heat_capacity[0], gas_capacity[stage_1 + 1 : -1] / heat_capacity[1]


def _stage_1(stage_1_column: NDArray[np.float64]) -> int:
    """Index of stage 1 along the chain, from per-node values of the stage-1 column."""
    return len(stage_1_column) + 1


def _along_chain(stage_1_column: NDArray[np.float64], stage_2_column: NDArray[np.float64]) -> NDArray[np.float64]:
    """Places per-node values of the two columns along the chain, with zeros at the plates."""
    return np.concatenate(([0.0], stage_1_column, [0.0], stage_2_column, [0.0]))


def _gas_mass_along_chain(
    temperature: NDArray[np.float64],
    pore_volume: tuple[NDArray[np.float64], NDArray[np.float64]],
    pulse_tube_volume: tuple[float, float],
    pressure: float,
) -> NDArray[np.float64]:
    """Helium mass in kg at each element of the chain, with each pulse tube's gas placed at its stage."""
    mass = _along_chain(*pore_volume) * helium.density(temperature, pressure)
    for stage, volume in zip((_stage_1(pore_volume[0]), -1), pulse_tube_volume, strict=True):
        mass[stage] += volume * _mean_tube_density(temperature[0], temperature[stage], pressure)
    return mass


def _mean_tube_density(warm: float, cold: float, pressure: float) -> float:
    """Helium density in kg/m^3 averaged over a linear temperature profile from `warm` to `cold`.

    The average is integral(density dT) / integral(dT), with both integrals taken by quadrature in log temperature.
    """
    points, weights = roots_legendre(48)
    temperature = np.exp(np.log(cold) + (points + 1) / 2 * np.log(warm / cold))
    return float(weights @ (helium.density(temperature, pressure) * temperature) / (weights @ temperature))
