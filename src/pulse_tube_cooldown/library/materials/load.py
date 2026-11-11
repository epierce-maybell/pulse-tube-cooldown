"""Loading materials from TOML files.

Each file describes one material by its name, its density in kg/m^3, and two properties, specific heat and thermal
conductivity:

    name = "OFHC Cu RRR100"
    density = 8960.0

    [specific_heat]
    form = "log_poly"
    coefficients = [...]
    range = [4.0, 300.0]
    source = "Marquardt, Le & Radebaugh (2000), NIST Monograph 177"

    [thermal_conductivity]
    ...

Thermal conductivity is optional.

Each property is a table read by `make_fit`. Within its range, a property takes the fitted value; outside, it
takes the value at the nearest end of the range and issues a RuntimeWarning naming the material, property, and side.
Python's default warning filter shows each distinct warning once.

Each property also has an integral over temperature from the low end of its range: the specific-heat integral and,
where there is a thermal conductivity, the conductivity integral. Within the range, an integral is a cubic Hermite
spline through Simpson-rule integrals at knots, with the property as its slopes. Knots are log-spaced and include any
tabulated points. Outside the range, it extends linearly with the property at the nearest end and warns as the
properties do.
"""

import os
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

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
    specific_heat = make_fit(data["specific_heat"], f"{name} specific_heat")
    conductivity = conductivity_integral = None

    if "thermal_conductivity" in data:
        conductivity = make_fit(data["thermal_conductivity"], f"{name} thermal_conductivity")
        conductivity_integral = _integrate(conductivity, data["thermal_conductivity"], f"{name} conductivity_integral")

    return Material(
        name=name,
        density=float(data["density"]),
        specific_heat=specific_heat,
        specific_heat_integral=_integrate(specific_heat, data["specific_heat"], f"{name} specific_heat_integral"),
        thermal_conductivity=conductivity,
        conductivity_integral=conductivity_integral,
    )


def _integrate(fit: PropertyFit, spec: Mapping[str, Any], label: str) -> PropertyFit:
    """Builds a property's integral from the low end of its range, extended linearly outside the range."""
    valid_range: Sequence[float] = spec["range"]
    lower, upper = valid_range
    point_temperature = np.array(spec.get("table", spec.get("measured", [])), dtype=np.float64).reshape(-1, 2)[:, 0]
    if "low_temperature" in spec:
        extension = spec["low_temperature"]
        point_temperature = np.union1d(
            point_temperature, [*(row[0] for row in extension["table"]), extension["join_temperature"]]
        )
    knots = np.union1d(
        np.geomspace(lower, upper, 1025), point_temperature[(lower < point_temperature) & (point_temperature < upper)]
    )

    slope = fit(knots)
    midpoint_value = fit((knots[:-1] + knots[1:]) / 2.0)
    simpson_integral = np.diff(knots) / 6.0 * (slope[:-1] + 4.0 * midpoint_value + slope[1:])
    spline = CubicHermiteSpline(knots, np.concatenate(([0.0], np.cumsum(simpson_integral))), slope)

    def integral(temperature: NDArray[np.float64]) -> NDArray[np.float64]:
        warn_below(temperature, valid_range, label, f"extending linearly from {lower:g} K")
        warn_above(temperature, valid_range, label, f"extending linearly from {upper:g} K")
        clipped = np.clip(temperature, lower, upper)
        return spline(clipped) + spline(clipped, 1) * (temperature - clipped)

    return integral
