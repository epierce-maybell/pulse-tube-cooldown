"""Drive of one gas manifold, from `docs/adrs/coupled-capacity-maps.md`.

The heads and compressors on a manifold share its supply and return pressures. The mean pressure P is the one at which
the heads and the manifold's warm volume at the warm temperature T_0 hold the gas charge. The swing dP between supply
P + dP / 2 and return P - dP / 2 has the heads take the compressors' flow,

    sum of head flows = rho(T_0, P - dP / 2) * sum of displacement rates,

for helium's density rho, up to the smallest of the compressors' relief swings.
"""

from collections.abc import Callable, Sequence

from scipy.optimize import brentq

from pulse_tube_cooldown.library.gases import helium

__all__ = ["calc_mean_pressure", "calc_pressure_swing"]


def calc_mean_pressure(
    head_gas_mass: Callable[[float], float],
    charge: float,
    warm_volume: float,
    warm_temperature: float,
    pressure_bracket: tuple[float, float] = (1.0e4, 5.0e6),
) -> float:
    """Calculates the mean pressure at which a manifold holds its gas charge.

    Args:
        head_gas_mass: Helium mass in kg in every head on the manifold at a pressure in Pa.
        charge: Helium mass in kg on the manifold.
        warm_volume: Gas volume in m^3 of the manifold's lines and buffers.
        warm_temperature: Temperature in K of the warm volume.
        pressure_bracket: Pressures in Pa that bound the mean pressure.

    Returns:
        Mean pressure in Pa.
    """

    def excess_mass(pressure: float) -> float:
        return head_gas_mass(pressure) + warm_volume * float(helium.density(warm_temperature, pressure)) - charge

    return brentq(excess_mass, *pressure_bracket)


def calc_pressure_swing(
    head_flow: Callable[[float, float], float],
    mean_pressure: float,
    warm_temperature: float,
    displacement_rate: Sequence[float],
    relief_swing: Sequence[float],
) -> float:
    """Calculates the pressure swing at which a manifold's heads take its compressors' flow.

    Args:
        head_flow: Mass flow in kg/s into every head on the manifold at a supply and a return pressure in Pa.
        mean_pressure: Mean of supply and return pressures in Pa.
        warm_temperature: Temperature in K of the gas the compressors draw.
        displacement_rate: Volumetric displacement rate in m^3/s of each compressor.
        relief_swing: Relief swing in Pa of each compressor.

    Returns:
        Supply minus return pressure in Pa.
    """
    total_displacement_rate = sum(displacement_rate)

    def excess_flow(swing: float) -> float:
        return_pressure = mean_pressure - swing / 2
        compressor_flow = float(helium.density(warm_temperature, return_pressure)) * total_displacement_rate
        return head_flow(mean_pressure + swing / 2, return_pressure) - compressor_flow

    relief_cap = min(relief_swing)
    if excess_flow(relief_cap) <= 0:
        return relief_cap
    return brentq(excess_flow, 0.0, relief_cap)
