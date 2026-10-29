"""He-4 properties from HEPAK, evaluated elementwise over temperature and pressure.

Each function takes temperature in K and pressure in Pa as arrays or scalars that broadcast together, and returns an
array of their broadcast shape in SI units. HEPAK raises `HepakError` for a state it cannot evaluate.
"""

import numpy as np
from hepak.properties import density_TP, heat_capacity_TP, properties_TP
from numpy.typing import ArrayLike, NDArray

__all__ = ["density", "isobaric_expansion_coefficient", "isobaric_specific_heat"]


def _expansion_coefficient(temperature: float, pressure: float) -> float:
    return properties_TP(temperature, pressure)["expansion_coefficient"] / temperature


_density = np.vectorize(lambda temperature, pressure: density_TP(temperature, pressure, basis="mass"), otypes=[float])
_specific_heat = np.vectorize(
    lambda temperature, pressure: heat_capacity_TP(temperature, pressure, basis="mass"), otypes=[float]
)
_expansion = np.vectorize(_expansion_coefficient, otypes=[float])


def density(temperature: ArrayLike, pressure: ArrayLike) -> NDArray[np.float64]:
    """Mass density in kg/m^3."""
    return _density(temperature, pressure)


def isobaric_specific_heat(temperature: ArrayLike, pressure: ArrayLike) -> NDArray[np.float64]:
    """Specific heat at constant pressure in J/(kg K)."""
    return _specific_heat(temperature, pressure)


def isobaric_expansion_coefficient(temperature: ArrayLike, pressure: ArrayLike) -> NDArray[np.float64]:
    """Isobaric expansion coefficient, -(1 / density) d(density)/dT at constant pressure, in 1/K."""
    return _expansion(temperature, pressure)
