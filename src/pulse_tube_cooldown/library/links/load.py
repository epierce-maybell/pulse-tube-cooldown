"""Links that transfer heat between nodes, built from configuration.

A link's configuration is a table holding its `type`, a key of `LINK_TYPES`, and the fields that type's loader reads:

    name = "Rod"
    type = "conduction"
    nodes = ["CFP", "STP"]
    material = "304 SS"
    area = 1.0e-4
    length = 0.1
"""

from collections.abc import Callable, Mapping
from typing import Any

from pulse_tube_cooldown.core.errors import ConfigurationError
from pulse_tube_cooldown.core.types import Link
from pulse_tube_cooldown.library.links.conduction import load_conduction
from pulse_tube_cooldown.library.links.multi_layer_insulation import load_multi_layer_insulation
from pulse_tube_cooldown.library.links.radiation import load_radiation

__all__ = ["LINK_TYPES", "load_link"]

LINK_TYPES: Mapping[str, Callable[[Mapping[str, Any]], Link]] = {
    "conduction": load_conduction,
    "multi_layer_insulation": load_multi_layer_insulation,
    "radiation": load_radiation,
}


def load_link(config: Mapping[str, Any]) -> Link:
    """Builds a link from its configuration table with the loader registered for its type.

    Args:
        config: The link's `type` and the fields that type's loader reads.

    Returns:
        The link.

    Raises:
        ConfigurationError: The type is unknown, a key is missing, or the type's loader rejects the table.
    """
    link_type = config.get("type")

    if link_type not in LINK_TYPES:
        raise ConfigurationError(
            f"Link {config.get('name')!r} has unknown type {link_type!r}; known types are {', '.join(LINK_TYPES)}."
        )

    try:
        return LINK_TYPES[link_type](config)
    except KeyError as error:
        raise ConfigurationError(f"Link {config.get('name')!r} is missing {error}.") from error
