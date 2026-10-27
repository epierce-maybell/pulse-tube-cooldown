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

`form` selects a fit in `PROPERTY_FITS`, and the remaining keys besides `range` and `source` are that fit's arguments.
`range` is the fit's validity range in K, and `source` cites the data behind it. Within its range, a property takes the
fitted value; outside, it takes the value at the nearest end of the range and issues a RuntimeWarning naming the
material, property, and side. Python's default warning filter shows each distinct warning once.

Each material also has a conductivity integral, the integral of thermal conductivity from the low end of its range.
Within the range, it is a cubic Hermite spline through Simpson-rule integrals at knots, with the conductivity as its
slopes. Knots are log-spaced and include any table points. Outside the range, it extends linearly with the
conductivity at the nearest end and warns as the properties do.
"""

import inspect
import math
import os
import tomllib
import warnings
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray
from scipy.interpolate import CubicHermiteSpline

from pulse_tube_cooldown.core.errors import ConfigurationError
from pulse_tube_cooldown.core.types import Material
from pulse_tube_cooldown.library.property_fits import PROPERTY_FITS, PropertyFit

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
        ConfigurationError: If no file or more than one file defines the name, any file in the directory is not a
            readable material file, or the material's definition is malformed.
    """
    path, data = _find_material(name, Path(directory))

    specific_heat, specific_heat_range, _ = _parse_property(path, data, "specific_heat")
    conductivity, conductivity_range, table_temperature = _parse_property(path, data, "thermal_conductivity")

    return Material(
        name=name,
        specific_heat=_clamp(specific_heat, specific_heat_range, name, "specific_heat"),
        thermal_conductivity=_clamp(conductivity, conductivity_range, name, "thermal_conductivity"),
        conductivity_integral=_integrate(conductivity, conductivity_range, table_temperature, name),
    )


def _find_material(name: str, directory: Path) -> tuple[Path, dict[str, Any]]:
    """Returns the file defining `name` and its contents, raising ConfigurationError unless exactly one does."""
    files: dict[str, tuple[Path, dict[str, Any]]] = {}

    for path in sorted(directory.glob("*.toml")):
        data = _read_toml(path)
        material_name = data.get("name")

        if not isinstance(material_name, str):
            raise ConfigurationError(f"{path}: a material file needs a string name; got {material_name!r}.")

        if material_name in files:
            raise ConfigurationError(
                f"Material {material_name!r} is defined in both {files[material_name][0]} and {path}."
            )

        files[material_name] = (path, data)

    if name not in files:
        available = ", ".join(repr(material_name) for material_name in sorted(files)) or "none"
        raise ConfigurationError(f"No material named {name!r} in {directory}; available: {available}.")

    return files[name]


def _read_toml(path: Path) -> dict[str, Any]:
    """Parses a TOML file, raising ConfigurationError if it is invalid."""
    try:
        with path.open("rb") as file:
            return tomllib.load(file)
    except tomllib.TOMLDecodeError as error:
        raise ConfigurationError(f"{path}: invalid TOML: {error}") from error


def _parse_property(
    path: Path, data: dict[str, Any], key: str
) -> tuple[PropertyFit, tuple[float, float], NDArray[np.float64]]:
    """Builds the property's fit and returns it with its range and its table temperatures, if tabulated."""
    where = f"{path} [{key}]"
    spec = data.get(key)

    if not isinstance(spec, dict):
        raise ConfigurationError(f"{where}: missing, or not a table.")

    parameters = dict(spec)
    form = parameters.pop("form", None)
    lower, upper = _parse_range(where, parameters.pop("range", None))
    source = parameters.pop("source", None)

    if form not in PROPERTY_FITS:
        raise ConfigurationError(f"{where}: unknown form {form!r}; known forms are {', '.join(sorted(PROPERTY_FITS))}.")

    if not isinstance(source, str) or not source.strip():
        raise ConfigurationError(f"{where}: source must be a nonempty string citing the property data; got {source!r}.")

    make_fit = PROPERTY_FITS[form]

    try:
        inspect.signature(make_fit).bind(**parameters)
        fit = make_fit(**parameters)
    except (TypeError, ValueError) as error:
        raise ConfigurationError(f"{where}: invalid parameters for the {form} fit: {error}") from error

    table_temperature = np.array(parameters["table"], dtype=np.float64)[:, 0] if form == "table" else np.array([])

    return fit, (lower, upper), table_temperature


def _clamp(fit: PropertyFit, valid_range: tuple[float, float], name: str, key: str) -> PropertyFit:
    """Evaluates the fit at the nearest end of its range, with a warning, outside that range."""
    lower, upper = valid_range
    range_text = f"its range [{lower:g}, {upper:g}] K"
    below = f"{name} {key} evaluated below {range_text}; using its value at {lower:g} K."
    above = f"{name} {key} evaluated above {range_text}; using its value at {upper:g} K."

    def clamped_fit(temperature: NDArray[np.float64]) -> NDArray[np.float64]:
        _warn_outside_range(temperature, lower, upper, below, above)
        return fit(np.clip(temperature, lower, upper))

    return clamped_fit


def _integrate(
    conductivity: PropertyFit, valid_range: tuple[float, float], table_temperature: NDArray[np.float64], name: str
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

    range_text = f"its range [{lower:g}, {upper:g}] K"
    below = f"{name} conductivity_integral evaluated below {range_text}; extending linearly from {lower:g} K."
    above = f"{name} conductivity_integral evaluated above {range_text}; extending linearly from {upper:g} K."

    def conductivity_integral(temperature: NDArray[np.float64]) -> NDArray[np.float64]:
        _warn_outside_range(temperature, lower, upper, below, above)
        clipped = np.clip(temperature, lower, upper)
        return spline(clipped) + spline(clipped, 1) * (temperature - clipped)

    return conductivity_integral


def _warn_outside_range(temperature: NDArray[np.float64], lower: float, upper: float, below: str, above: str) -> None:
    """Warns, on behalf of the property's caller, for temperatures below or above the range."""
    if np.any(temperature < lower):
        warnings.warn(below, RuntimeWarning, stacklevel=3)
    if np.any(temperature > upper):
        warnings.warn(above, RuntimeWarning, stacklevel=3)


def _parse_range(where: str, value: object) -> tuple[float, float]:
    """Converts a range of [lower, upper] in K, raising ConfigurationError unless 0 < lower < upper."""
    if not (
        isinstance(value, list)
        and len(value) == 2
        and all(isinstance(bound, int | float) and not isinstance(bound, bool) for bound in value)
    ):
        raise ConfigurationError(f"{where}: range must be [lower, upper] in K; got {value!r}.")

    lower, upper = float(value[0]), float(value[1])

    if not 0 < lower < upper < math.inf:
        raise ConfigurationError(f"{where}: range must satisfy 0 < lower < upper; got {value!r}.")

    return lower, upper
