"""Reading library files, chosen by their `name` field, from directories of TOML files."""

import os
import tomllib
from pathlib import Path
from typing import Any

from pulse_tube_cooldown.core.errors import ConfigurationError

__all__ = ["read_named_toml"]


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
