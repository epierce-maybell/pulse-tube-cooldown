"""One pulse-tube head's gas inventory and mass swing, from `docs/adrs/03-coupled-capacity-maps.md`.

A head's temperatures run along its chain: the warm flange, the stage-1 section's nodes from its warm end, stage 1, the
stage-2 section's nodes, and stage 2. Each node's pore gas is at the node's temperature, each pulse tube's gas along a
linear profile from the warm flange to its stage, and the dead volume and reservoirs at the warm flange's. Values along
the chain place the dead volume at the warm flange and each pulse tube and its reservoir at its stage.

Over a cycle the pressure swings by dP about the mean pressure P. Pore gas swings at its node's temperature. The other
gas swings at constant entropy from its temperature at P, each reservoir by dP / sqrt(1 + (omega R B)^2) for the angular
valve frequency omega, its valve resistance R, and its pressure response B = V_r / (gamma P), with V_r its volume and
gamma helium's heat capacity ratio at the warm flange's temperature and P.
"""

import math
from collections.abc import Callable

import numpy as np
from numpy.polynomial.legendre import leggauss
from numpy.typing import NDArray

from pulse_tube_cooldown.library.gases import helium
from pulse_tube_cooldown.library.heads import HeadModel

__all__ = ["calc_gas_mass", "calc_mass_swing"]


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
    pores = _pore_volume_along_chain(model) @ helium.density(temperature, pressure)
    tubes = sum(
        tube.volume
        * _mean_over_tube(lambda t: helium.density(t, pressure), warm, temperature[stage], quadrature_points)
        for tube, stage in zip(model.pulse_tubes, _stages(model), strict=True)
    )
    warm_volume = model.dead_volume + sum(tube.reservoir_volume for tube in model.pulse_tubes)
    return float(pores + tubes + warm_volume * helium.density(warm, pressure))


def calc_mass_swing(
    model: HeadModel,
    temperature: NDArray[np.float64],
    mean_pressure: float,
    pressure_swing: float,
    *,
    resistance: tuple[float, float],
    frequency: float,
    quadrature_points: int = 16,
) -> NDArray[np.float64]:
    """Calculates the mass entering each element of a head's chain over one cycle.

    Args:
        model: The head's design.
        temperature: Temperatures in K along the head's chain.
        mean_pressure: Mean of supply and return pressures in Pa.
        pressure_swing: Supply minus return pressure in Pa.
        resistance: Stage-1 and stage-2 needle-valve resistances in Pa s/m^3.
        frequency: Valve frequency in Hz.
        quadrature_points: Gauss-Legendre points averaging each pulse tube's swing over its profile.

    Returns:
        Mass in kg per cycle at each element of the chain. The sum is the head's mass swing.
    """

    def adiabatic_swing(t: NDArray[np.float64] | float, swing: float = pressure_swing) -> NDArray[np.float64]:
        return helium.isentropic_density(t, mean_pressure, mean_pressure + swing / 2) - helium.isentropic_density(
            t, mean_pressure, mean_pressure - swing / 2
        )

    warm = temperature[0]
    mass = _pore_volume_along_chain(model) * (
        helium.density(temperature, mean_pressure + pressure_swing / 2)
        - helium.density(temperature, mean_pressure - pressure_swing / 2)
    )
    mass[0] += model.dead_volume * float(adiabatic_swing(warm))
    heat_capacity_ratio = float(helium.heat_capacity_ratio(warm, mean_pressure))
    for tube, stage, valve_resistance in zip(model.pulse_tubes, _stages(model), resistance, strict=True):
        response = (
            2 * math.pi * frequency * valve_resistance * tube.reservoir_volume / (heat_capacity_ratio * mean_pressure)
        )
        mass[stage] += tube.volume * _mean_over_tube(adiabatic_swing, warm, temperature[stage], quadrature_points)
        mass[stage] += tube.reservoir_volume * float(adiabatic_swing(warm, pressure_swing / math.sqrt(1 + response**2)))
    return mass


def _stages(model: HeadModel) -> tuple[int, int]:
    """Indices of stage 1 and stage 2 along the chain."""
    return 1 + _node_count(model, 0), 2 + _node_count(model, 0) + _node_count(model, 1)


def _node_count(model: HeadModel, section: int) -> int:
    return sum(layer.node_count for layer in model.regenerators[section].layers)


def _pore_volume_along_chain(model: HeadModel) -> NDArray[np.float64]:
    """Pore gas volume in m^3 at each element of the chain, zero at the warm flange and stages."""
    first, second = (
        [layer.pore_volume / layer.node_count for layer in section.layers for _ in range(layer.node_count)]
        for section in model.regenerators
    )
    return np.array([0.0, *first, 0.0, *second, 0.0])


def _mean_over_tube(
    function: Callable[[NDArray[np.float64]], NDArray[np.float64]], warm: float, cold: float, points: int
) -> float:
    """Mean of `function` over a linear temperature profile from `warm` to `cold`.

    The mean is integral(function dT) / integral(dT), with both integrals taken by quadrature in log temperature.
    """
    nodes, weights = leggauss(points)
    temperature = np.exp(np.log(cold) + (nodes + 1) / 2 * np.log(warm / cold))
    return float(weights @ (function(temperature) * temperature) / (weights @ temperature))
