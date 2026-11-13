"""He-4 properties from HEPAK, evaluated elementwise over temperature and pressure.

Each function takes temperature in K and pressure in Pa as arrays or scalars that broadcast together, and returns an
array of their broadcast shape in SI units. HEPAK raises `HepakError` for a state it cannot evaluate. HEPAK's thermal
conductivity starts at `CONDUCTIVITY_MINIMUM_TEMPERATURE`; below it, conductivity takes its value there at the same
pressure.
"""

from dataclasses import dataclass

import numpy as np
from hepak.properties import (
    _call_calc,
    _offset,
    density_TP,
    enthalpy_TP,
    entropy_TP,
    heat_capacity_TP,
    properties_TP,
    thermal_conductivity_TP,
)
from numpy.typing import ArrayLike, NDArray

__all__ = [
    "HeliumProperties",
    "density",
    "enthalpy",
    "heat_capacity_ratio",
    "isentropic_density",
    "isobaric_expansion_coefficient",
    "isobaric_specific_heat",
    "properties",
    "thermal_conductivity",
]

CONDUCTIVITY_MINIMUM_TEMPERATURE = 3.5


@dataclass(frozen=True, slots=True, kw_only=True)
class HeliumProperties:
    """Helium properties at a set of temperatures and pressures.

    Attributes:
        density: Mass density.
        isobaric_specific_heat: Specific heat at constant pressure.
        isobaric_expansion_coefficient: -(1 / density) d(density)/dT at constant pressure.
        isothermal_compressibility: (1 / density) d(density)/dp at constant temperature.
        heat_capacity_ratio: Ratio of isobaric to isochoric specific heat, gamma.
        thermal_conductivity: Thermal conductivity.
    """

    density: NDArray[np.float64]
    isobaric_specific_heat: NDArray[np.float64]
    isobaric_expansion_coefficient: NDArray[np.float64]
    isothermal_compressibility: NDArray[np.float64]
    heat_capacity_ratio: NDArray[np.float64]
    thermal_conductivity: NDArray[np.float64]


def _expansion_coefficient(temperature: float, pressure: float) -> float:
    return properties_TP(temperature, pressure)["expansion_coefficient"] / temperature


def _isentropic_density(temperature: float, reference_pressure: float, pressure: float) -> float:
    entropy = entropy_TP(temperature, reference_pressure, basis="mass")
    return _call_calc("P", pressure, "S", entropy)[_offset(3)]


_density = np.vectorize(lambda temperature, pressure: density_TP(temperature, pressure, basis="mass"), otypes=[float])
_enthalpy = np.vectorize(lambda temperature, pressure: enthalpy_TP(temperature, pressure, basis="mass"), otypes=[float])
_specific_heat = np.vectorize(
    lambda temperature, pressure: heat_capacity_TP(temperature, pressure, basis="mass"), otypes=[float]
)
_expansion = np.vectorize(_expansion_coefficient, otypes=[float])
_isentropic = np.vectorize(_isentropic_density, otypes=[float])
_ratio = np.vectorize(lambda temperature, pressure: properties_TP(temperature, pressure)["gamma"], otypes=[float])
_conductivity = np.vectorize(thermal_conductivity_TP, otypes=[float])


def density(temperature: ArrayLike, pressure: ArrayLike) -> NDArray[np.float64]:
    """Mass density in kg/m^3."""
    return _density(temperature, pressure)


def enthalpy(temperature: ArrayLike, pressure: ArrayLike) -> NDArray[np.float64]:
    """Specific enthalpy in J/kg."""
    return _enthalpy(temperature, pressure)


def isentropic_density(
    temperature: ArrayLike, reference_pressure: ArrayLike, pressure: ArrayLike
) -> NDArray[np.float64]:
    """Mass density in kg/m^3 at `pressure`, at the entropy of `temperature` and `reference_pressure`."""
    return _isentropic(temperature, reference_pressure, pressure)


def heat_capacity_ratio(temperature: ArrayLike, pressure: ArrayLike) -> NDArray[np.float64]:
    """Ratio of isobaric to isochoric specific heat, gamma."""
    return _ratio(temperature, pressure)


def isobaric_specific_heat(temperature: ArrayLike, pressure: ArrayLike) -> NDArray[np.float64]:
    """Specific heat at constant pressure in J/(kg K)."""
    return _specific_heat(temperature, pressure)


def isobaric_expansion_coefficient(temperature: ArrayLike, pressure: ArrayLike) -> NDArray[np.float64]:
    """Isobaric expansion coefficient, -(1 / density) d(density)/dT at constant pressure, in 1/K."""
    return _expansion(temperature, pressure)


def thermal_conductivity(temperature: ArrayLike, pressure: ArrayLike) -> NDArray[np.float64]:
    """Thermal conductivity in W/(m K)."""
    return _conductivity(np.maximum(temperature, CONDUCTIVITY_MINIMUM_TEMPERATURE), pressure)


def properties(temperature: ArrayLike, pressure: ArrayLike) -> HeliumProperties:
    """Return the properties in `HeliumProperties` from one HEPAK evaluation per state, in the broadcast shape."""
    temperature, pressure = np.broadcast_arrays(
        np.asarray(temperature, dtype=np.float64), np.asarray(pressure, dtype=np.float64)
    )
    states = [properties_TP(t, p) for t, p in zip(temperature.flat, pressure.flat, strict=True)]

    def collect(key: str) -> NDArray[np.float64]:
        return np.array([state[key] for state in states]).reshape(temperature.shape)

    conductivity = collect("thermal_conductivity")
    below_range = temperature < CONDUCTIVITY_MINIMUM_TEMPERATURE
    conductivity[below_range] = thermal_conductivity(CONDUCTIVITY_MINIMUM_TEMPERATURE, pressure[below_range])

    return HeliumProperties(
        density=collect("density"),
        isobaric_specific_heat=collect("cp"),
        isobaric_expansion_coefficient=collect("expansion_coefficient") / temperature,
        isothermal_compressibility=collect("isothermal_compressibility"),
        heat_capacity_ratio=collect("gamma"),
        thermal_conductivity=conductivity,
    )
