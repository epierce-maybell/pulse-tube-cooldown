"""Pulse-tube cooling as a thermal-network link over one gas manifold."""

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

import pulse_tube_cooldown.library.links.cooler.drive as drive
import pulse_tube_cooldown.library.links.cooler.head as head

__all__ = ["PulseTubeCooler"]


@dataclass(frozen=True, slots=True, kw_only=True)
class PulseTubeCooler:
    """One manifold's heads and compressors as a thermal-network link.

    Attributes:
        name: Unique link name in the graph.
        warm_flange: Node taking the heads' rejected heat and the warm gas temperature.
        heads: Units on this manifold, including their column descriptions and stage names.
        compressors: Compressors whose displacement rates add, under the smallest relief swing.
        charge: Total helium mass in kg, including the warm volume and every head.
        warm_volume: Lines' and buffers' gas volume in m^3, at the warm-flange temperature.
        waveform_factor: Pressure waveform factor kappa, absorbed into fitted conductances at its default of one.
        pressure_bracket: Pressures in Pa bounding the mean-pressure solve.
    """

    name: str
    warm_flange: str
    heads: tuple[head.PulseTubeHead, ...]
    compressors: tuple[drive.Compressor, ...]
    charge: float
    warm_volume: float
    waveform_factor: float = 1.0
    pressure_bracket: tuple[float, float] = (1.0e4, 2.5e6)

    @property
    def nodes(self) -> tuple[str, ...]:
        """Unique node names in head-chain order, starting at the shared warm flange."""
        return tuple(dict.fromkeys(node for unit in self.heads for node in unit.chain_nodes(self.warm_flange)))

    def calc_heat_flow(self, temperature: NDArray[np.float64]) -> NDArray[np.float64]:
        """Calculates heat flows in W into `nodes`, from their temperatures in K.

        All heads share the mean pressure set by the helium inventory and the swing set by the compressor drive.
        Contributions to shared nodes add, and the returned flows sum to zero.
        """
        node_index = {name: index for index, name in enumerate(self.nodes)}
        indices = [np.array([node_index[node] for node in unit.chain_nodes(self.warm_flange)]) for unit in self.heads]
        temperatures = [temperature[index] for index in indices]
        volumes = [unit.pore_volume for unit in self.heads]
        warm_temperature = float(temperature[node_index[self.warm_flange]])

        def gas_mass(pressure: float) -> float:
            return sum(
                head.calc_gas_mass(t, pore, unit.pulse_tube_volumes, pressure)
                for unit, t, pore in zip(self.heads, temperatures, volumes, strict=True)
            )

        def mass_flow(supply: float, back: float) -> float:
            return sum(
                float(head.calc_mass_flow(t, pore, unit.pulse_tube_volumes, supply, back, unit.frequency)[0])
                for unit, t, pore in zip(self.heads, temperatures, volumes, strict=True)
            )

        mean_pressure = drive.calc_mean_pressure(
            gas_mass, self.charge, self.warm_volume, warm_temperature, self.pressure_bracket
        )
        swing = drive.calc_pressure_swing(
            mass_flow,
            mean_pressure,
            warm_temperature,
            [compressor.displacement_rate for compressor in self.compressors],
            [compressor.relief_swing for compressor in self.compressors],
        )
        heat_flow = np.zeros(len(self.nodes), dtype=np.float64)
        for unit, t, pore, index in zip(self.heads, temperatures, volumes, indices, strict=True):
            flow = head.calc_mass_flow(
                t, pore, unit.pulse_tube_volumes, mean_pressure + swing / 2, mean_pressure - swing / 2, unit.frequency
            )
            first, second = (
                head.calc_enthalpy_flow(conductance, swing, self.waveform_factor) for conductance in unit.conductances
            )
            np.add.at(
                heat_flow,
                index,
                head.calc_heat_flow(
                    t,
                    unit.calc_heat_capacity(t),
                    flow,
                    (first, second),
                    unit.packing_constants,
                    mean_pressure,
                    unit.frequency,
                ),
            )
        return heat_flow
