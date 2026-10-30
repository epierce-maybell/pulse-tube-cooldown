"""Heat flow through a multi-layer insulation blanket, from the Lockheed correlation.

The blanket's heat flux in W/m^2 is the sum of conduction through the spacers and radiation between the layers,

    q = C_s * n^b * T_m * (T_h - T_c) / (N + 1) + C_r * e * (T_h^4.67 - T_c^4.67) / N,

for warm and cold boundary temperatures T_h and T_c in K, their mean T_m, N layers packed n layers per cm, and the
layers' room-temperature emissivity e. The spacer coefficient C_s, layer density exponent b, radiation coefficient
C_r, and emissivity e depend on the blanket's construction. Keller, Cunnington, and Glassford, "Thermal Performance of
Multi-Layer Insulations," NASA CR-134477 (1974).
"""

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import NDArray

__all__ = ["MultiLayerInsulation", "load_multi_layer_insulation"]


@dataclass(frozen=True, slots=True, kw_only=True)
class MultiLayerInsulation:
    """A blanket wrapped around the inner node, facing the outer node.

    Attributes:
        name: Unique link identifier within a graph.
        nodes: Names of the inner (wrapped) and outer nodes, in that order.
        area: Blanket area in m^2.
        layers: Number of reflective layers.
        layer_density: Layers per cm of blanket thickness.
        spacer_coefficient: C_s in W m^-2 K^-2 (layers/cm)^-b.
        layer_density_exponent: b.
        radiation_coefficient: C_r in W m^-2 K^-4.67.
        emissivity: Room-temperature emissivity of a layer.
        source: Where the blanket's coefficients come from.
    """

    name: str
    nodes: tuple[str, str]
    area: float
    layers: int
    layer_density: float
    spacer_coefficient: float
    layer_density_exponent: float
    radiation_coefficient: float
    emissivity: float
    source: str

    def calc_heat_flow(self, temperature: NDArray[np.float64]) -> NDArray[np.float64]:
        """Calculates the heat flow into each of the two nodes.

        Args:
            temperature: Temperatures in K of the inner and outer nodes.

        Returns:
            Heat flows in W, in the order of `nodes`. Positive values indicate heat flow into the node.
        """
        t_inner, t_outer = temperature
        spacers = (
            self.spacer_coefficient
            * self.layer_density**self.layer_density_exponent
            * (t_outer + t_inner)
            / 2.0
            * (t_outer - t_inner)
            / (self.layers + 1)
        )
        radiation = self.radiation_coefficient * self.emissivity * (t_outer**4.67 - t_inner**4.67) / self.layers
        flow = self.area * (spacers + radiation)
        return np.array([flow, -flow], dtype=np.float64)


def load_multi_layer_insulation(config: Mapping[str, Any]) -> MultiLayerInsulation:
    """Builds a multi-layer insulation link from its configuration table.

    Args:
        config: The link's `name`, `type`, `nodes` as `[inner, outer]`, `area` in m^2, `layers`, `layer_density` in
            layers per cm, and the blanket's `spacer_coefficient`, `layer_density_exponent`, `radiation_coefficient`,
            `emissivity`, and `source`.

    Raises:
        KeyError: A key is missing.
        ValueError: `nodes` does not hold exactly two entries.
    """
    inner, outer = config["nodes"]
    return MultiLayerInsulation(
        name=config["name"],
        nodes=(inner, outer),
        area=config["area"],
        layers=config["layers"],
        layer_density=config["layer_density"],
        spacer_coefficient=config["spacer_coefficient"],
        layer_density_exponent=config["layer_density_exponent"],
        radiation_coefficient=config["radiation_coefficient"],
        emissivity=config["emissivity"],
        source=config["source"],
    )
