"""Heat flow across a pressed joint, such as a bolted flange or a clamp.

The joint's conductance in W/K scales with the force pressing it together and as a power of temperature,

    G(T) = g * F * (T / T_r)^n,

for clamping force F in N, conductance per unit force g in W/(K N) at reference temperature T_r in K, and temperature
exponent n. The heat flow is G at the mean of the two node temperatures times their difference. Dhuley, "Pressed
copper and gold-plated copper contacts at low temperatures," Cryogenics 101 (2019) 111-124, surveys g and n.
"""

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import NDArray

from pulse_tube_cooldown.core.types import EnergyBalanceTerms

__all__ = ["JointContact", "load_joint_contact"]


@dataclass(frozen=True, slots=True, kw_only=True)
class JointContact:
    """A pressed joint between two nodes.

    Attributes:
        name: Unique link identifier within a graph.
        nodes: Names of the two nodes the joint presses together.
        force: Clamping force in N.
        conductance_per_force: Conductance per unit clamping force in W/(K N) at `reference_temperature`.
        reference_temperature: Temperature in K at which `conductance_per_force` is given.
        temperature_exponent: Power of temperature the conductance scales with.
        source: Where the joint's coefficients come from.
    """

    name: str
    nodes: tuple[str, str]
    force: float
    conductance_per_force: float
    reference_temperature: float
    temperature_exponent: float
    source: str

    def calc_energy_balance_terms(self, temperature: NDArray[np.float64]) -> EnergyBalanceTerms:
        """Calculates the heat flow into each of the two nodes.

        Args:
            temperature: Temperatures in K of the two nodes, in the order of `nodes`.

        Returns:
            Heat flows in W, in the order of `nodes`. Positive values indicate heat flow into the node.
        """
        first, second = temperature
        mean_temperature = (first + second) / 2.0
        conductance = (
            self.conductance_per_force
            * self.force
            * (mean_temperature / self.reference_temperature) ** self.temperature_exponent
        )
        flow = conductance * (second - first)
        return EnergyBalanceTerms(heat_flow=np.array([flow, -flow], dtype=np.float64))


def load_joint_contact(config: Mapping[str, Any]) -> JointContact:
    """Builds a pressed-joint link from its configuration table.

    Args:
        config: The link's `name`, `type`, `nodes`, clamping `force` in N, `conductance_per_force` in W/(K N) at
            `reference_temperature` in K, `temperature_exponent`, and `source`.

    Raises:
        KeyError: A key is missing.
        ValueError: `nodes` does not hold exactly two entries.
    """
    first, second = config["nodes"]
    return JointContact(
        name=config["name"],
        nodes=(first, second),
        force=config["force"],
        conductance_per_force=config["conductance_per_force"],
        reference_temperature=config["reference_temperature"],
        temperature_exponent=config["temperature_exponent"],
        source=config["source"],
    )
