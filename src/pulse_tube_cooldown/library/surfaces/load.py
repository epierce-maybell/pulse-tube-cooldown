"""Loading surfaces from TOML files.

Each file describes one surface finish by its name and its emissivity:

    name = "Polished OFHC Cu"

    [emissivity]
    form = "table"
    table = [[4.0, 0.006], [77.0, 0.015], [300.0, 0.03]]
    range = [4.0, 300.0]
    source = "..."

Emissivity is a table read by `make_fit`. Within its range, it takes the fitted value; outside, it takes the value at
the nearest end of the range and issues a RuntimeWarning naming the surface and side.
"""

import os
from pathlib import Path

from pulse_tube_cooldown.core.types import Surface
from pulse_tube_cooldown.library.property_fits import make_fit
from pulse_tube_cooldown.library.toml_files import read_named_toml

__all__ = ["load_surface"]

_DATA_DIRECTORY = Path(__file__).parent / "data"


def load_surface(name: str, /, *, directory: str | os.PathLike[str] = _DATA_DIRECTORY) -> Surface:
    """Loads the surface with the given name from a directory of surface files.

    Args:
        name: The surface's `name` field, such as "Polished OFHC Cu".
        directory: The directory whose `*.toml` files are searched. Defaults to the library's surface data.

    Returns:
        The surface, with an emissivity that evaluates elementwise over temperature arrays of any shape.

    Raises:
        ConfigurationError: If no file or more than one file defines the name.
    """
    data = read_named_toml(name, directory=directory, kind="surfaces")

    return Surface(name=name, emissivity=make_fit(data["emissivity"], f"{name} emissivity"))
