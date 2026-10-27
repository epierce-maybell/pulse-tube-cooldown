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

from pulse_tube_cooldown.core.errors import ConfigurationError
from pulse_tube_cooldown.core.types import Material
from pulse_tube_cooldown.library.property_fits import PROPERTY_FITS, PropertyFit

__all__ = ["load_material"]

_DATA_DIRECTORY = Path(__file__).parent / "data"

_PROPERTIES = ("specific_heat", "thermal_conductivity")


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

    unknown = data.keys() - {"name", *_PROPERTIES}
    if unknown:
        raise ConfigurationError(
            f"{path}: unknown keys {sorted(unknown)}; a material has name, specific_heat, and thermal_conductivity."
        )

    return Material(
        name=name,
        specific_heat=_load_property(path, data, name, "specific_heat"),
        thermal_conductivity=_load_property(path, data, name, "thermal_conductivity"),
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


def _load_property(path: Path, data: dict[str, Any], name: str, key: str) -> PropertyFit:
    """Builds the property's fit, evaluated at the nearest end of its range, with a warning, outside that range."""
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

    range_text = f"its range [{lower:g}, {upper:g}] K"
    below = f"{name} {key} evaluated below {range_text}; using its value at {lower:g} K."
    above = f"{name} {key} evaluated above {range_text}; using its value at {upper:g} K."

    def clamped_fit(temperature: NDArray[np.float64]) -> NDArray[np.float64]:
        if np.any(temperature < lower):
            warnings.warn(below, RuntimeWarning, stacklevel=2)
        if np.any(temperature > upper):
            warnings.warn(above, RuntimeWarning, stacklevel=2)
        return fit(np.clip(temperature, lower, upper))

    return clamped_fit


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
