"""One pulse-tube head's gas inventory, mass swing, and heat flows, from `docs/adrs/03-coupled-capacity-maps.md`.

A head's temperatures run along its chain: the warm flange, the stage-1 section's nodes from its warm end, stage 1, the
stage-2 section's nodes, and stage 2. Each node's pore gas is at the node's temperature, each stage's heat-exchanger
void at the stage's, each pulse tube's gas along a linear profile from the warm flange to its stage, and the dead volume
and reservoirs at the warm flange's. Values along the chain place the dead volume at the warm flange and each pulse
tube, its reservoir, and its heat-exchanger void at its stage.

Over a cycle the pressure swings by dP about the mean pressure P. Pore gas swings at its node's temperature and void
gas at its stage's. The other gas swings at constant entropy from its temperature at P, each reservoir by
dP / sqrt(1 + (omega R B)^2) for the angular valve frequency omega, its valve resistance R, and its pressure response
B = V_r / (gamma P), with V_r its volume and gamma helium's heat capacity ratio at the warm flange's temperature and P.

Each pulse tube draws acoustic power E = e P_a^2 at its cold end, with pressure amplitude P_a = kappa dP / 2 and
e = omega^2 R B^2 / (2 [1 + (omega R B)^2]).

Gas that crosses a whole pulse tube in one cycle carries its enthalpy to the stage. Pulse tube i loses

    L_i = f [h(T_0) - h(T_i)] s(dm_i - lambda m_i)

to stage i, for helium's specific enthalpy h at P and the warm flange's and stage's temperatures, the mass dm_i entering
the tube's reservoir through its warm end over a cycle, the tube's gas mass m_i at P, and the buffer fraction lambda.
The smooth positive part s(x) = x Phi(x / w) + w phi(x / w), for the standard normal distribution Phi and its density
phi, is the mean excess over the buffer of a stroke spread normally with standard deviation w = sigma m_i.

Within a regenerator section, temperature is linear in position between the warm end, the node centres, and the cold
end. At each node boundary, with x increasing towards the warm end, the energy flow towards the cold end is

    F = (1 - beta T / (1 + epsilon)) E + a (m c_p)^2 / (f C'_m) dT/dx + A dLambda/dx,

for helium's isobaric expansion coefficient beta and specific heat c_p at the boundary's temperature T, the section's
packing coefficient a, the mass flow m of every gas volume on the boundary's cold side, the valve frequency f, and the
conductivity integral Lambda of the section's wall of area A. C'_m is the matrix heat capacity of the half-nodes on
either side of the boundary over the distance between the points on either side, at temperatures T_a and T_b. Each
half-node takes its material's mean specific heat between them, (H(T_a) - H(T_b)) / (T_a - T_b) for the specific-heat
integral H. epsilon is the half-nodes' pore gas heat capacity, at their node temperatures, over their matrix's.
Acoustic power grows towards the warm end as d(ln E) = beta T / (1 + epsilon_k) d(ln T) within node k, from E_2 at
stage 2 and from E_1 plus the stage-2 section's warm-end power at stage 1.

A node's heat flow is the energy flow at its warm boundary less that at its cold boundary. Stage i also takes the
section energy flows on either side, its tube wall's conduction from the warm flange, L_i, and -E_i. Stage 1 also takes
c E_2 x (1 - x), for the design's stage coupling c and stage 2's gross-cooling fraction x = beta T / (1 + epsilon) at
the stage-2 boundary. The warm flange's heat flow is zero.
"""

import math
from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
from numpy.polynomial.legendre import leggauss
from numpy.typing import NDArray
from scipy.special import ndtr

from pulse_tube_cooldown.core.types import Material
from pulse_tube_cooldown.library.gases import helium
from pulse_tube_cooldown.library.heads import HeadModel, PulseTube, RegeneratorSection

__all__ = [
    "Head",
    "calc_acoustic_power",
    "calc_breakthrough_loss",
    "calc_gas_mass",
    "calc_heat_flow",
    "calc_mass_swing",
    "calc_pore_heat_capacity",
]


@dataclass(frozen=True, slots=True, kw_only=True)
class Head:
    """One pulse-tube head: its design and the values particular to the unit.

    Attributes:
        model: The head's design.
        packing: Packing coefficients a of the stage-1 and stage-2 regenerator sections.
        resistance: Stage-1 and stage-2 needle-valve resistances R in Pa s/m^3.
        frequency: Valve frequency f in Hz.
    """

    model: HeadModel
    packing: tuple[float, float]
    resistance: tuple[float, float]
    frequency: float


def calc_gas_mass(
    model: HeadModel, temperature: NDArray[np.float64], pressure: float, *, quadrature_points: int = 16
) -> float:
    """Calculates the helium mass in a head at one pressure.

    Args:
        model: The head's design.
        temperature: Temperatures in K along the head's chain.
        pressure: Pressure in Pa throughout the head.
        quadrature_points: Gauss-Legendre points averaging each pulse tube's density over its profile.

    Returns:
        Helium mass in kg.
    """
    warm = temperature[0]
    isothermal = _isothermal_volume_along_chain(model) @ helium.density(temperature, pressure)
    tubes = sum(
        _tube_gas_mass(tube, warm, temperature[stage], pressure, quadrature_points)
        for tube, stage in zip(model.pulse_tubes, _stages(model), strict=True)
    )
    warm_volume = model.dead_volume + sum(tube.reservoir_volume for tube in model.pulse_tubes)
    return float(isothermal + tubes + warm_volume * helium.density(warm, pressure))


def calc_pore_heat_capacity(model: HeadModel, temperature: NDArray[np.float64], pressure: float) -> NDArray[np.float64]:
    """Calculates the heat capacity of the gas at each element's temperature along a head's chain.

    That gas is the pore gas at each node and the heat-exchanger void at each stage.

    Args:
        model: The head's design.
        temperature: Temperatures in K along the head's chain.
        pressure: Pressure in Pa throughout the head.

    Returns:
        Heat capacity in J/K at each element of the chain, zero at the warm flange.
    """
    return (
        _isothermal_volume_along_chain(model)
        * helium.density(temperature, pressure)
        * helium.isobaric_specific_heat(temperature, pressure)
    )


def calc_mass_swing(
    head: Head,
    temperature: NDArray[np.float64],
    mean_pressure: float,
    pressure_swing: float,
    *,
    quadrature_points: int = 16,
) -> NDArray[np.float64]:
    """Calculates the mass entering each element of a head's chain over one cycle.

    Args:
        head: The head.
        temperature: Temperatures in K along the head's chain.
        mean_pressure: Mean of supply and return pressures in Pa.
        pressure_swing: Supply minus return pressure in Pa.
        quadrature_points: Gauss-Legendre points averaging each pulse tube's swing over its profile.

    Returns:
        Mass in kg per cycle at each element of the chain. The sum is the head's mass swing.
    """

    def adiabatic_swing(t: NDArray[np.float64]) -> NDArray[np.float64]:
        return _adiabatic_density_swing(t, mean_pressure, pressure_swing)

    model = head.model
    warm = temperature[0]
    mass = _isothermal_volume_along_chain(model) * (
        helium.density(temperature, mean_pressure + pressure_swing / 2)
        - helium.density(temperature, mean_pressure - pressure_swing / 2)
    )
    mass[0] += model.dead_volume * float(adiabatic_swing(warm))
    reservoirs = _reservoir_mass_swing(head, warm, mean_pressure, pressure_swing)
    for tube, stage, reservoir in zip(model.pulse_tubes, _stages(model), reservoirs, strict=True):
        mass[stage] += tube.volume * _mean_over_tube(adiabatic_swing, warm, temperature[stage], quadrature_points)
        mass[stage] += reservoir
    return mass


def calc_acoustic_power(
    head: Head, mean_pressure: float, pressure_swing: float, warm_temperature: float
) -> tuple[float, float]:
    """Calculates the acoustic power each pulse tube draws at its cold end.

    Args:
        head: The head.
        mean_pressure: Mean of supply and return pressures in Pa.
        pressure_swing: Supply minus return pressure in Pa.
        warm_temperature: Warm-flange temperature in K.

    Returns:
        Stage-1 and stage-2 acoustic power in W.
    """
    amplitude = head.model.waveform_factor * pressure_swing / 2
    omega = 2 * math.pi * head.frequency
    first, second = (
        amplitude**2 * omega**2 * resistance * response**2 / (2 * (1 + (omega * resistance * response) ** 2))
        for resistance, response in zip(
            head.resistance, _pressure_responses(head.model, warm_temperature, mean_pressure), strict=True
        )
    )
    return first, second


def calc_breakthrough_loss(
    head: Head, temperature: NDArray[np.float64], mean_pressure: float, pressure_swing: float
) -> tuple[float, float]:
    """Calculates the heat each pulse tube's breakthrough carries to its stage.

    Args:
        head: The head.
        temperature: Temperatures in K along the head's chain.
        mean_pressure: Mean of supply and return pressures in Pa.
        pressure_swing: Supply minus return pressure in Pa.

    Returns:
        Stage-1 and stage-2 breakthrough loss in W.
    """
    model = head.model
    warm = temperature[0]
    stages = _stages(model)
    warm_enthalpy, *stage_enthalpy = helium.enthalpy(temperature[[0, *stages]], mean_pressure)
    first, second = (
        head.frequency
        * (warm_enthalpy - enthalpy)
        * _smooth_positive_part(stroke - model.buffer_fraction * tube_mass, model.stroke_spread * tube_mass)
        for tube, stage, enthalpy, stroke in zip(
            model.pulse_tubes,
            stages,
            stage_enthalpy,
            _reservoir_mass_swing(head, warm, mean_pressure, pressure_swing),
            strict=True,
        )
        for tube_mass in [_tube_gas_mass(tube, warm, temperature[stage], mean_pressure)]
    )
    return first, second


def calc_heat_flow(
    head: Head, temperature: NDArray[np.float64], mean_pressure: float, pressure_swing: float
) -> NDArray[np.float64]:
    """Calculates the heat flow into each element of a head's chain.

    Args:
        head: The head.
        temperature: Temperatures in K along the head's chain.
        mean_pressure: Mean of supply and return pressures in Pa.
        pressure_swing: Supply minus return pressure in Pa.

    Returns:
        Heat flow in W into each element of the chain.
    """
    model = head.model
    stage_1, stage_2 = _stages(model)
    warm = temperature[0]
    mass_swing = calc_mass_swing(head, temperature, mean_pressure, pressure_swing)
    cold_side_flow = head.frequency * np.cumsum(mass_swing[::-1])[::-1]
    tube_power = calc_acoustic_power(head, mean_pressure, pressure_swing, warm)
    loss = calc_breakthrough_loss(head, temperature, mean_pressure, pressure_swing)
    second, second_warm_power, stage_2_cooling_fraction = _energy_flow(
        model.regenerators[1],
        temperature[stage_1 : stage_2 + 1],
        cold_side_flow[stage_1 + 1 : stage_2 + 1],
        cold_power=tube_power[1],
        packing=head.packing[1],
        frequency=head.frequency,
        pressure=mean_pressure,
    )
    first, _, _ = _energy_flow(
        model.regenerators[0],
        temperature[: stage_1 + 1],
        cold_side_flow[1 : stage_1 + 1],
        cold_power=tube_power[0] + second_warm_power,
        packing=head.packing[0],
        frequency=head.frequency,
        pressure=mean_pressure,
    )
    tube_conduction = []
    for tube, stage in zip(model.pulse_tubes, (stage_1, stage_2), strict=True):
        integral = tube.wall_material.conductivity_integral(np.array([warm, temperature[stage]]))
        tube_conduction.append(tube.wall_area / tube.length * float(integral[0] - integral[1]))
    coupling = model.stage_coupling * tube_power[1] * stage_2_cooling_fraction * (1 - stage_2_cooling_fraction)

    heat = np.zeros_like(temperature, dtype=float)
    heat[1:stage_1] = first[:-1] - first[1:]
    heat[stage_1 + 1 : stage_2] = second[:-1] - second[1:]
    heat[stage_1] = first[-1] + tube_conduction[0] + loss[0] - tube_power[0] - second[0] + coupling
    heat[stage_2] = second[-1] + tube_conduction[1] + loss[1] - tube_power[1]
    return heat


def _adiabatic_density_swing(
    temperature: NDArray[np.float64] | float, mean_pressure: float, pressure_swing: float
) -> NDArray[np.float64]:
    """Density swing in kg/m^3 at constant entropy from `temperature` at the mean pressure."""
    return helium.isentropic_density(
        temperature, mean_pressure, mean_pressure + pressure_swing / 2
    ) - helium.isentropic_density(temperature, mean_pressure, mean_pressure - pressure_swing / 2)


def _reservoir_mass_swing(
    head: Head, warm_temperature: float, mean_pressure: float, pressure_swing: float
) -> list[float]:
    """Mass in kg entering each pulse tube's reservoir over one cycle."""
    omega = 2 * math.pi * head.frequency
    return [
        tube.reservoir_volume
        * float(
            _adiabatic_density_swing(
                warm_temperature, mean_pressure, pressure_swing / math.sqrt(1 + (omega * resistance * response) ** 2)
            )
        )
        for tube, resistance, response in zip(
            head.model.pulse_tubes,
            head.resistance,
            _pressure_responses(head.model, warm_temperature, mean_pressure),
            strict=True,
        )
    ]


def _smooth_positive_part(value: float, width: float) -> float:
    """Mean of max(X, 0) for X normally distributed with mean `value` and standard deviation `width`."""
    scaled = value / width
    return float(value * ndtr(scaled) + width * math.exp(-(scaled**2) / 2) / math.sqrt(2 * math.pi))


def _energy_flow(
    section: RegeneratorSection,
    temperature: NDArray[np.float64],
    mass_flow: NDArray[np.float64],
    *,
    cold_power: float,
    packing: float,
    frequency: float,
    pressure: float,
) -> tuple[NDArray[np.float64], float, float]:
    """Energy flow at a section's node boundaries, acoustic power at its warm end, and gross cooling at its cold end.

    Args:
        section: The regenerator section.
        temperature: Temperatures in K of the warm end, the nodes from the warm end, and the cold end.
        mass_flow: Mass flow in kg/s through each node boundary from the warm end.
        cold_power: Acoustic power in W at the cold end.
        packing: Packing coefficient a.
        frequency: Valve frequency in Hz.
        pressure: Mean pressure in Pa.

    Returns:
        Energy flow in W towards the cold end at each node boundary from the warm end, acoustic power in W at the warm
        end, and the gross-cooling fraction beta T / (1 + epsilon) at the cold end.
    """
    node_layers = [layer for layer in section.layers for _ in range(layer.node_count)]
    length = section.length * np.array([layer.length_fraction / layer.node_count for layer in node_layers])
    edge = np.concatenate(([0.0], np.cumsum(length)))
    point = np.concatenate(([0.0], edge[:-1] + length / 2, [edge[-1]]))
    node_temperature = temperature[1:-1]

    matrix_capacity = np.array(
        [
            layer.mass / layer.node_count * layer.material.specific_heat(np.array([t])).item()
            for layer, t in zip(node_layers, node_temperature, strict=True)
        ]
    )
    gas_capacity = (
        np.array([layer.pore_volume / layer.node_count for layer in node_layers])
        * helium.density(node_temperature, pressure)
        * helium.isobaric_specific_heat(node_temperature, pressure)
    )
    node_ratio = gas_capacity / matrix_capacity
    half_gas = np.concatenate(([0.0], gas_capacity / 2, [0.0]))
    boundary_matrix = np.zeros(len(node_layers) + 1)
    for node, layer in enumerate(node_layers):
        boundary_matrix[node : node + 2] += (
            layer.mass
            / layer.node_count
            / 2
            * _mean_specific_heat(layer.material, temperature[node : node + 2], temperature[node + 1 : node + 3])
        )
    boundary_ratio = (half_gas[:-1] + half_gas[1:]) / boundary_matrix

    span = np.diff(point)
    boundary_temperature = temperature[:-1] + np.diff(temperature) * (edge - point[:-1]) / span
    gradient = -np.diff(temperature) / span
    boundary_beta_t = boundary_temperature * helium.isobaric_expansion_coefficient(boundary_temperature, pressure)
    node_beta_t = node_temperature * helium.isobaric_expansion_coefficient(node_temperature, pressure)

    warm_half = (boundary_beta_t[:-1] + node_beta_t) / 2 * np.log(boundary_temperature[:-1] / node_temperature)
    cold_half = (node_beta_t + boundary_beta_t[1:]) / 2 * np.log(node_temperature / boundary_temperature[1:])
    log_power_step = (warm_half + cold_half) / (1 + node_ratio)
    power = cold_power * np.exp(np.concatenate((np.cumsum(log_power_step[::-1])[::-1], [0.0])))

    capacity_rate = mass_flow * helium.isobaric_specific_heat(boundary_temperature, pressure)
    cooling_fraction = boundary_beta_t / (1 + boundary_ratio)
    pressure_flow = (1 - cooling_fraction) * power
    regeneration = packing * capacity_rate**2 / (frequency * boundary_matrix / span) * gradient
    conduction = -section.wall_area * np.diff(section.wall_material.conductivity_integral(temperature)) / span
    return pressure_flow + regeneration + conduction, float(power[0]), float(cooling_fraction[-1])


def _mean_specific_heat(
    material: Material, warm: NDArray[np.float64], cold: NDArray[np.float64]
) -> NDArray[np.float64]:
    """Mean specific heat in J/(kg K) between `warm` and `cold`, (H(warm) - H(cold)) / (warm - cold).

    Where the two temperatures agree to 1 part in 10^6, it is the specific heat at their midpoint.
    """
    difference = warm - cold
    close = np.abs(difference) <= 1e-6 * warm
    secant = (material.specific_heat_integral(warm) - material.specific_heat_integral(cold)) / np.where(
        close, 1.0, difference
    )
    return np.where(close, material.specific_heat((warm + cold) / 2), secant)


def _pressure_responses(model: HeadModel, warm_temperature: float, mean_pressure: float) -> list[float]:
    """Pressure response B in m^3/Pa of each pulse tube's reservoir."""
    heat_capacity_ratio = float(helium.heat_capacity_ratio(warm_temperature, mean_pressure))
    return [tube.reservoir_volume / (heat_capacity_ratio * mean_pressure) for tube in model.pulse_tubes]


def _stages(model: HeadModel) -> tuple[int, int]:
    """Indices of stage 1 and stage 2 along the chain."""
    return 1 + _node_count(model, 0), 2 + _node_count(model, 0) + _node_count(model, 1)


def _node_count(model: HeadModel, section: int) -> int:
    return sum(layer.node_count for layer in model.regenerators[section].layers)


def _isothermal_volume_along_chain(model: HeadModel) -> NDArray[np.float64]:
    """Volume in m^3 of the gas at each element's temperature: pores at nodes, voids at stages, none at the flange."""
    first, second = (
        [layer.pore_volume / layer.node_count for layer in section.layers for _ in range(layer.node_count)]
        for section in model.regenerators
    )
    stage_1_void, stage_2_void = model.heat_exchanger_voids
    return np.array([0.0, *first, stage_1_void, *second, stage_2_void])


def _tube_gas_mass(tube: PulseTube, warm: float, cold: float, pressure: float, points: int = 16) -> float:
    """Gas mass in kg in a pulse tube at one pressure, along a linear temperature profile from `warm` to `cold`."""
    return tube.volume * _mean_over_tube(lambda t: helium.density(t, pressure), warm, cold, points)


def _mean_over_tube(
    function: Callable[[NDArray[np.float64]], NDArray[np.float64]], warm: float, cold: float, points: int
) -> float:
    """Mean of `function` over a linear temperature profile from `warm` to `cold`.

    The mean is integral(function dT) / integral(dT), with both integrals taken by quadrature in log temperature.
    """
    nodes, weights = leggauss(points)
    temperature = np.exp(np.log(cold) + (nodes + 1) / 2 * np.log(warm / cold))
    return float(weights @ (function(temperature) * temperature) / (weights @ temperature))
