"""The compressors and heads sharing one gas supply."""

from dataclasses import dataclass

import numpy as np

from pulse_tube_cooldown.library.gases import helium

__all__ = ["Manifold"]


@dataclass(frozen=True, slots=True, kw_only=True)
class Manifold:
    """The compressors and heads sharing one gas supply, all heads identical and in the same state.

    Attributes:
        head_count: Number of heads on the manifold.
        compressor_count: Number of compressors on the manifold.
        displacement_rate: Volume each compressor displaces per unit time.
        volumetric_efficiency: Fraction of the displacement rate delivered at the return-line density.
        relief_amplitude: Largest pressure amplitude the relief valves allow.
        fill_pressure: Pressure of the whole manifold at rest, at the warm-flange temperature.
        supply_volume: Volume at supply pressure: supply lines and the compressors' high-pressure side.
        return_volume: Volume at return pressure: return lines and the compressors' low-pressure side.
    """

    head_count: int
    compressor_count: int
    displacement_rate: float
    volumetric_efficiency: float
    relief_amplitude: float
    fill_pressure: float
    supply_volume: float
    return_volume: float

    def pressures(
        self,
        *,
        warm_temperature: float,
        head_volume: float,
        head_mass: float,
        head_mass_slope: float,
        head_admittance: float,
    ) -> tuple[float, float]:
        """Return the mean pressure that conserves the helium charge and the pressure amplitude it sustains.

        Args:
            warm_temperature: Temperature of the supply and return gas, and of every volume at rest.
            head_volume: Gas volume of one head, its reservoirs included.
            head_mass: Gas mass of one head at the fill pressure.
            head_mass_slope: Rate of change of one head's gas mass with mean pressure.
            head_admittance: Magnitude of one head's mass admittance, mass-flow amplitude per pressure amplitude.
        """
        warm = helium.properties(warm_temperature, self.fill_pressure)
        density = float(warm.density)
        density_slope = density * float(warm.isothermal_compressibility)

        mass_slope = density_slope * (self.supply_volume + self.return_volume) + self.head_count * head_mass_slope
        line_imbalance = density_slope * (self.supply_volume - self.return_volume)
        head_excess = self.head_count * (head_mass - density * head_volume)

        delivery = np.pi * self.compressor_count * self.volumetric_efficiency * self.displacement_rate
        gain = delivery / (self.head_count * head_admittance)
        amplitude_at_fill = gain * density / (1.0 + gain * density_slope)
        amplitude_slope = gain * density_slope / (1.0 + gain * density_slope)

        excess_pressure = -(line_imbalance * amplitude_at_fill + head_excess) / (
            mass_slope + line_imbalance * amplitude_slope
        )
        amplitude = amplitude_at_fill + amplitude_slope * excess_pressure
        if amplitude > self.relief_amplitude:
            amplitude = self.relief_amplitude
            excess_pressure = -(line_imbalance * amplitude + head_excess) / mass_slope

        return self.fill_pressure + excess_pressure, amplitude
