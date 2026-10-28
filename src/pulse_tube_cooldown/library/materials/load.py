"""Loading materials from TOML files.

Each file describes one material by its name and two properties, specific heat and thermal conductivity:

    name = "OFHC Cu RRR100"

    [specific_heat]
    form = "log_poly"
    coefficients = [...]
    range = [4.0, 300.0]
    source = "Marquardt, Le & Radebaugh (2000), NIST Monograph 177"

    [thermal_conductivity]
    ...

Each property is a table read by `make_fit`. Within its range, a property takes the fitted value; outside, it
takes the value at the nearest end of the range and issues a RuntimeWarning naming the material, property, and side.
Python's default warning filter shows each distinct warning once.

Each material also has a conductivity integral, the integral of thermal conductivity from the low end of its range.
Within the range, it is a cubic Hermite spline through Simpson-rule integrals at knots, with the conductivity as its
slopes. Knots are log-spaced and include any table points. Outside the range, it extends linearly with the
conductivity at the nearest end and warns as the properties do.
"""

import os
from collections.abc import Sequence
from pathlib import Path

import numpy as np
from numpy.typing import NDArray
from scipy.interpolate import CubicHermiteSpline

from pulse_tube_cooldown.core.types import Material
from pulse_tube_cooldown.library.property_fits import PropertyFit, make_fit, warn_above, warn_below
from pulse_tube_cooldown.library.toml_files import read_named_toml

__all__ = ["load_material"]

_DATA_DIRECTORY = Path(__file__).parent / "data"


def load_material(name: str, /, *, directory: str | os.PathLike[str] = _DATA_DIRECTORY) -> Material:
    """Loads the material with the given name from a directory of material files.

    Args:
        name: The material's `name` field, such as "OFHC Cu RRR100".
        directory: The directory whose `*.toml` files are searched. Defaults to the library's material data.

    Returns:
        The material, with properties that evaluate elementwise over temperature arrays of any shape.

    Raises:
        ConfigurationError: If no file or more than one file defines the name.
    """
    data = read_named_toml(name, directory=directory, kind="materials")
    conductivity = make_fit(data["thermal_conductivity"], f"{name} thermal_conductivity")
    table = np.array(data["thermal_conductivity"].get("table", []), dtype=np.float64).reshape(-1, 2)

    return Material(
        name=name,
        specific_heat=make_fit(data["specific_heat"], f"{name} specific_heat"),
        thermal_conductivity=conductivity,
        conductivity_integral=_integrate(
            conductivity, data["thermal_conductivity"]["range"], table[:, 0], f"{name} conductivity_integral"
        ),
    )


def _integrate(
    conductivity: PropertyFit, valid_range: Sequence[float], table_temperature: NDArray[np.float64], label: str
) -> PropertyFit:
    """Builds the conductivity integral from the low end of the range, extended linearly outside the range."""
    lower, upper = valid_range
    knots = np.union1d(
        np.geomspace(lower, upper, 1025), table_temperature[(lower < table_temperature) & (table_temperature < upper)]
    )

    slope = conductivity(knots)
    midpoint_conductivity = conductivity((knots[:-1] + knots[1:]) / 2.0)
    simpson_integral = np.diff(knots) / 6.0 * (slope[:-1] + 4.0 * midpoint_conductivity + slope[1:])
    spline = CubicHermiteSpline(knots, np.concatenate(([0.0], np.cumsum(simpson_integral))), slope)

    def conductivity_integral(temperature: NDArray[np.float64]) -> NDArray[np.float64]:
        warn_below(temperature, valid_range, label, f"extending linearly from {lower:g} K")
        warn_above(temperature, valid_range, label, f"extending linearly from {upper:g} K")
        clipped = np.clip(temperature, lower, upper)
        return spline(clipped) + spline(clipped, 1) * (temperature - clipped)

    return conductivity_integral
