"""Thermal properties and heat balances of a graph at given node temperatures.

Temperature arrays hold one value in K per node, in node index order.
"""

from collections.abc import Sequence

import numpy as np
from numpy.typing import NDArray

from pulse_tube_cooldown.core.types import Node

__all__ = ["calc_heat_capacity"]


def calc_heat_capacity(nodes: Sequence[Node], temperature: NDArray[np.float64]) -> NDArray[np.float64]:
    """Calculates the heat capacity of each node at its temperature.

    Args:
        nodes: Nodes in index order.
        temperature: Temperature in K of each node.

    Returns:
        Heat capacity in J/K of each node.
    """
    return np.array(
        [
            sum(mass * material.specific_heat(temperature[i : i + 1])[0] for material, mass in node.composition)
            for i, node in enumerate(nodes)
        ],
        dtype=np.float64,
    )
