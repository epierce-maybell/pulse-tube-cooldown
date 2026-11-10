"""Drive of one gas manifold, from `docs/adrs/archive/03-coupled-capacity-maps.md`.

The heads, compressors, and warm lines on a manifold share its supply and return pressures. The mean pressure P is the
one at which the manifold's gas volumes hold its helium charge. The pressure swing dP between supply P + dP / 2 and
return P - dP / 2 follows from the supply equation,

    sum of head mass flows = rho(T_0, P - dP / 2) * sum of displacement rates,

for helium's density rho at the warm-flange temperature T_0. The relief valve of the compressor with the lowest relief
swing dP_r holds the swing at dP_r. With dP_s the swing that solves the supply equation, the swing is the smooth minimum

    dP = dP_s - w * ln(1 + exp((dP_s - dP_r) / w)),

for a relief width w.
"""

from collections.abc import Callable, Sequence
from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq

from pulse_tube_cooldown.library.gases import helium
from pulse_tube_cooldown.library.heads import HeadModel

__all__ = ["Compressor", "calc_charge", "calc_mean_pressure", "calc_pressure_swing"]


@dataclass(frozen=True, slots=True, kw_only=True)
class Compressor:
    """One compressor on a manifold.

    Attributes:
        displacement_rate: Volume per second drawn in at the low-pressure inlet, in m^3/s.
        relief_swing: Supply minus return pressure in Pa at which the relief valve opens.
    """

    displacement_rate: float
    relief_swing: float


def calc_charge(
    models: Sequence[HeadModel], *, warm_volume: float, fill_pressure: float, fill_temperature: float
) -> float:
    """Calculates a manifold's helium charge from a uniform fill at rest.

    Args:
        models: Design of each head on the manifold.
        warm_volume: Gas volume in m^3 of the manifold's lines and buffers.
        fill_pressure: Absolute pressure in Pa throughout the manifold at rest.
        fill_temperature: Temperature in K of every gas volume at rest.

    Returns:
        Helium mass in kg.
    """
    volume = warm_volume + sum(_gas_volume(model) for model in models)
    return volume * float(helium.density(fill_temperature, fill_pressure))


def calc_mean_pressure(
    head_gas_mass: Callable[[float], float],
    charge: float,
    *,
    warm_volume: float,
    warm_temperature: float,
    pressure_bracket: tuple[float, float] = (1.0e4, 2.5e6),
) -> float:
    """Calculates the mean pressure at which a manifold holds its helium charge.

    Args:
        head_gas_mass: Helium mass in kg in every head on the manifold at a pressure in Pa.
        charge: Helium mass in kg on the manifold.
        warm_volume: Gas volume in m^3 of the manifold's lines and buffers.
        warm_temperature: Temperature in K of the lines and buffers.
        pressure_bracket: Pressures in Pa that bound the mean pressure.

    Returns:
        Mean pressure in Pa.
    """

    def excess_mass(pressure: float) -> float:
        return head_gas_mass(pressure) + warm_volume * float(helium.density(warm_temperature, pressure)) - charge

    return brentq(excess_mass, *pressure_bracket)


def calc_pressure_swing(
    head_mass_flow: Callable[[float], float],
    mean_pressure: float,
    warm_temperature: float,
    *,
    compressors: Sequence[Compressor],
    relief_width: float = 5e4,
) -> float:
    """Calculates the pressure swing that the supply equation and the relief valve set.

    The search for the supply equation's swing keeps the return pressure above a tenth of the mean pressure and the
    swing within 20 relief widths above the relief swing, where the smooth minimum is the relief swing.

    Args:
        head_mass_flow: Mass flow in kg/s into every head on the manifold at a pressure swing in Pa, about
            `mean_pressure`.
        mean_pressure: Mean of supply and return pressures in Pa.
        warm_temperature: Temperature in K of the gas the compressors draw in.
        compressors: Compressors on the manifold.
        relief_width: Width w in Pa of the relief valve's transition.

    Returns:
        Supply minus return pressure in Pa.
    """
    displacement_rate = sum(compressor.displacement_rate for compressor in compressors)
    relief_swing = min(compressor.relief_swing for compressor in compressors)

    def excess_flow(swing: float) -> float:
        supply = float(helium.density(warm_temperature, mean_pressure - swing / 2)) * displacement_rate
        return head_mass_flow(swing) - supply

    upper = min(relief_swing + 20 * relief_width, 1.8 * mean_pressure)
    supplied = upper if excess_flow(upper) <= 0 else brentq(excess_flow, 0.0, upper)
    return supplied - relief_width * _softplus((supplied - relief_swing) / relief_width)


def _gas_volume(model: HeadModel) -> float:
    """Gas volume in m^3 of a head's pores, heat-exchanger voids, pulse tubes, reservoirs, and dead volume."""
    pores = sum(layer.pore_volume for section in model.regenerators for layer in section.layers)
    tubes = sum(tube.volume + tube.reservoir_volume for tube in model.pulse_tubes)
    return pores + sum(model.heat_exchanger_voids) + tubes + model.dead_volume


def _softplus(x: float) -> float:
    return float(np.logaddexp(0.0, x))
