"""Temperature-property interpolants and primitives for the candidate experiment."""

import tomllib
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from hepak.properties import properties_TP
from numpy.typing import NDArray
from scipy.interpolate import PchipInterpolator

from pulse_tube_cooldown.core.types import Material


@dataclass(frozen=True, slots=True)
class ThermalProperty:
    """A temperature property and its first two temperature primitives."""

    value: PchipInterpolator
    primitive: PchipInterpolator
    second_primitive: PchipInterpolator

    @classmethod
    def from_values(cls, temperature: NDArray, values: NDArray) -> "ThermalProperty":
        """Interpolates values with a shape-preserving cubic."""
        value = PchipInterpolator(temperature, values, extrapolate=False)
        primitive = value.antiderivative()
        return cls(value, primitive, primitive.antiderivative())


@dataclass(frozen=True, slots=True)
class HeliumTable:
    """HEPAK properties at one prescribed mean pressure, in SI units."""

    temperature: NDArray
    density: ThermalProperty
    cp: ThermalProperty
    beta: ThermalProperty
    viscosity: ThermalProperty
    kinematic_viscosity: ThermalProperty

    @classmethod
    def build(cls, temperature: NDArray, *, pressure: float) -> "HeliumTable":
        """Evaluates HEPAK once per temperature and tabulates local properties."""
        rows = [properties_TP(float(t), pressure, basis="mass") for t in temperature]
        density = np.array([p["density"] for p in rows])
        cp = np.array([p["cp"] for p in rows])
        beta = np.array([p["expansion_coefficient"] for p in rows]) / temperature
        viscosity = np.array([p["viscosity"] for p in rows])
        return cls(
            temperature,
            *(ThermalProperty.from_values(temperature, v) for v in (density, cp, beta, viscosity, viscosity / density)),
        )


def temperature_knots(numerics: dict, *, material_directory: Path) -> NDArray:
    """Includes the material tables' resolved features in a common grid."""
    knots = np.geomspace(numerics["temperature_min"], numerics["temperature_max"], numerics["property_knots"])
    for path in material_directory.glob("*.toml"):
        with path.open("rb") as stream:
            spec = tomllib.load(stream)["specific_heat"]
        for key in ("table", "measured"):
            if key in spec:
                knots = np.union1d(knots, np.asarray(spec[key])[:, 0])
    return knots[(knots >= numerics["temperature_min"]) & (knots <= numerics["temperature_max"])]


def heat_capacity_table(material: Material, temperature: NDArray) -> ThermalProperty:
    """Tabulates the existing material continuation and preserves measured knots."""
    return ThermalProperty.from_values(temperature, material.specific_heat(temperature))


def mean_integral(property_: ThermalProperty, left: NDArray, right: NDArray) -> NDArray:
    """Averages a property on a linear-temperature segment."""
    delta = right - left
    close = np.abs(delta) < 1e-3
    result = np.empty_like(delta)
    result[~close] = (property_.primitive(right[~close]) - property_.primitive(left[~close])) / delta[~close]
    middle = (left[close] + right[close]) / 2
    half = delta[close] / (2 * np.sqrt(3))
    result[close] = (property_.value(middle - half) + property_.value(middle + half)) / 2
    return result
