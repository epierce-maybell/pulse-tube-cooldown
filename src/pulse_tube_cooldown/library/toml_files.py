"""Reading TOML files: library files chosen by their `name` field, and data files named by path."""

import os
import tomllib
from pathlib import Path
from typing import Any

from pulse_tube_cooldown.core.errors import ConfigurationError

__all__ = ["read_named_toml", "read_toml"]


def read_named_toml(name: str, /, *, directory: str | os.PathLike[str], kind: str) -> dict[str, Any]:
    """Returns the contents of the one `*.toml` file in a directory whose `name` field is `name`.

    Args:
        name: The file's `name` field, such as "OFHC Cu RRR100".
        directory: The directory whose `*.toml` files are searched.
        kind: The plural noun for what the files describe, such as "materials", for error messages.

    Raises:
        ConfigurationError: If no file or more than one file has the name.
    """
    contents = [tomllib.loads(path.read_text(encoding="utf-8")) for path in Path(directory).glob("*.toml")]
    matches = [content for content in contents if content["name"] == name]

    if len(matches) != 1:
        raise ConfigurationError(f"Found {len(matches)} {kind} named {name!r} in {directory}; expected one.")

    return matches[0]


def read_toml(
    path: str | os.PathLike[str], /, *, root: str | os.PathLike[str] = Path(__file__).parents[3]
) -> dict[str, Any]:
    """Returns the contents of a TOML file.

    Args:
        path: The file's path, absolute or relative to `root`.
        root: The directory relative paths start from. Defaults to the repository root.
    """
    return tomllib.loads((Path(root) / path).read_text(encoding="utf-8"))
