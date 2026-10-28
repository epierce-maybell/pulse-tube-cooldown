"""Loading gases from TOML files.

Each file describes one gas by its name, five constants in SI units, and three temperature-dependent properties:

    name = "He-4"
    molar_mass = 4.002602e-3
    heat_capacity_ratio = 1.6666666666666667
    adsorption_energy = ...
    adsorption_constant = ...
    adsorption_capacity = ...

    [saturation_pressure]
    form = "wagner"
    critical_temperature = 5.1953
    ...
    range = [2.1768, inf]
    source = "..."

The properties are `accommodation_coeff`, `thermal_conductivity`, and `saturation_pressure`, each a table read by
`make_fit`. Outside its range, a property takes its value at the nearest end, with a RuntimeWarning. Saturation
pressure's range extends to infinity, where the Wagner form is infinite past the critical point.
"""

import os
from pathlib import Path

from pulse_tube_cooldown.core.types import Gas
from pulse_tube_cooldown.library.property_fits import make_fit
from pulse_tube_cooldown.library.toml_files import read_named_toml

__all__ = ["load_gas"]

_DATA_DIRECTORY = Path(__file__).parent / "data"


def load_gas(name: str, /, *, directory: str | os.PathLike[str] = _DATA_DIRECTORY) -> Gas:
    """Loads the gas with the given name from a directory of gas files.

    Args:
        name: The gas's `name` field, such as "He-4".
        directory: The directory whose `*.toml` files are searched. Defaults to the library's gas data.

    Returns:
        The gas, with properties that evaluate elementwise over temperature arrays of any shape.

    Raises:
        ConfigurationError: If no file or more than one file defines the name.
    """
    data = read_named_toml(name, directory=directory, kind="gases")

    return Gas(
        name=name,
        molar_mass=float(data["molar_mass"]),
        heat_capacity_ratio=float(data["heat_capacity_ratio"]),
        accommodation_coeff=make_fit(data["accommodation_coeff"], f"{name} accommodation_coeff"),
        thermal_conductivity=make_fit(data["thermal_conductivity"], f"{name} thermal_conductivity"),
        saturation_pressure=make_fit(data["saturation_pressure"], f"{name} saturation_pressure"),
        adsorption_energy=float(data["adsorption_energy"]),
        adsorption_constant=float(data["adsorption_constant"]),
        adsorption_capacity=float(data["adsorption_capacity"]),
    )
