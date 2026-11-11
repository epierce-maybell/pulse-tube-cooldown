"""A pulse-tube head as a reduced dynamical system: regenerator energies driven by stage and clamp temperatures."""

import itertools
from collections.abc import Callable, Iterator, Mapping
from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import NDArray
from scipy.interpolate import CubicHermiteSpline, PchipInterpolator, RBFInterpolator
from scipy.optimize import brentq

from pulse_tube_cooldown.core.errors import ConfigurationError
from pulse_tube_cooldown.core.types import EnergyBalanceTerms, Graph, Node
from pulse_tube_cooldown.engine.equations import calc_net_heat_flow
from pulse_tube_cooldown.library.links.pulse_tube.parameters import (
    HeadCalibration,
    HeadCoefficients,
    RegeneratorLayer,
    SharedParameters,
    load_shared_parameters,
    read_head_calibration,
)

__all__ = ["PulseTubeHead", "build_pulse_tube_head", "load_pulse_tube_head"]


@dataclass(frozen=True, slots=True, kw_only=True)
class Regenerator:
    """A regenerator's nodes from warm to cold, with their materials' specific energies.

    Attributes:
        nodes: The regenerator nodes.
        table_temperature: Temperatures in K at which `table_energy` is tabulated.
        table_energy: Each node's specific energy in J/kg at `table_temperature`, shape (nodes, knots).
        mass: Node masses in kg.
        midpoint: Node midpoints as fractions of the regenerator length from its warm end.
        positions: Quadrature positions in each node as fractions of the regenerator length, shape (nodes, order).
        weights: Quadrature weights summing to one.
    """

    nodes: tuple[Node, ...]
    table_temperature: NDArray[np.float64]
    table_energy: NDArray[np.float64]
    mass: NDArray[np.float64]
    midpoint: NDArray[np.float64]
    positions: NDArray[np.float64]
    weights: NDArray[np.float64]

    def specific_energy(self, temperature: NDArray[np.float64]) -> NDArray[np.float64]:
        """Returns each node's specific energy in J/kg at its temperature in K."""
        return np.array(
            [
                interpolate(value, self.table_temperature, energy)
                for value, energy in zip(temperature, self.table_energy, strict=True)
            ]
        )

    def target(self, profile: NDArray[np.float64]) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Returns each node's mean specific energy in J/kg over a profile at `positions`, and its temperature in K."""
        energy = np.array(
            [
                interpolate(temperature, self.table_temperature, table) @ self.weights
                for temperature, table in zip(profile, self.table_energy, strict=True)
            ]
        )
        temperature = np.array(
            [
                interpolate(value, table, self.table_temperature)
                for value, table in zip(energy, self.table_energy, strict=True)
            ]
        )
        return energy, temperature


@dataclass(frozen=True, slots=True, kw_only=True)
class ColdMap:
    """The load map as zero-clamp-load cooling power, continued outside its measured temperatures.

    Attributes:
        minimum_temperature: Lowest measured temperature of each stage, in K.
        maximum_temperature: Highest measured temperature of each stage, in K.
        tail_length: Length in K of each stage's exponential continuation.
        blend_width: Width `beta` of the warm blend, in span-normalized temperature.
        heater_power: Maps span-normalized stage temperatures, shape (n, 2), to heater powers in W.
        test_stand: The test stand the load map was measured on.
        test_stand_stages: Graph indices of the test stand's two stages.
        test_stand_temperature: Temperature in K of the test stand's other nodes.
    """

    minimum_temperature: NDArray[np.float64]
    maximum_temperature: NDArray[np.float64]
    tail_length: NDArray[np.float64]
    blend_width: float
    heater_power: RBFInterpolator
    test_stand: Graph
    test_stand_stages: NDArray[np.intp]
    test_stand_temperature: float

    def cooling_power(
        self, stage_temperature: NDArray[np.float64], warm_correction: NDArray[np.float64]
    ) -> NDArray[np.float64]:
        """Returns `q0 = H(B(T)) + L_stand(T) + w(T) * s` in W, given the warm correction `s` in W."""
        minimum, maximum, tail = self.minimum_temperature, self.maximum_temperature, self.tail_length
        span = maximum - minimum
        above = np.maximum(stage_temperature - maximum, 0.0)
        below = np.maximum(minimum - stage_temperature, 0.0)
        continued = (
            np.clip(stage_temperature, minimum, maximum)
            + tail * -np.expm1(-above / tail)
            - tail * -np.expm1(-below / tail)
        )
        heater_power = self.heater_power(((continued - minimum) / span)[None, :])[0]

        temperature = np.full(len(self.test_stand.nodes), self.test_stand_temperature)
        temperature[self.test_stand_stages] = stage_temperature
        test_stand_heating = calc_net_heat_flow(self.test_stand, temperature)[self.test_stand_stages]

        blend = -np.expm1(-np.sum(((above + below) / span) ** 2) / self.blend_width**2)
        return heater_power + test_stand_heating + blend * warm_correction


@dataclass(frozen=True, slots=True, kw_only=True)
class ClampResponse:
    """The changes with clamp load in stage temperatures, clamp temperature, and regenerator 2's target profile.

    Attributes:
        stage_shift: Maps clamp load in W to the stage temperature shifts `d_1` and `d_2` in K.
        clamp_increment: Maps clamp load in W to the clamp temperature increment `d_c` in K.
        profile_fraction: The unloaded profile fraction `F_c` at the clamp.
        unloaded_fraction: The unloaded profile fraction `F` at regenerator 2's quadrature positions.
        above_clamp: Whether each of regenerator 2's quadrature positions is at or above the clamp.
        loaded_shape: Maps clamp load in W to the profile shape `U` above the clamp and `V` below it, at regenerator
            2's quadrature positions.
    """

    stage_shift: PchipInterpolator
    clamp_increment: CubicHermiteSpline
    profile_fraction: float
    unloaded_fraction: NDArray[np.float64]
    above_clamp: NDArray[np.bool_]
    loaded_shape: CubicHermiteSpline

    def stage_shifts(self, clamp_load: float) -> NDArray[np.float64]:
        """Returns `d_1` and `d_2` in K, continued linearly beyond the tabulated loads."""
        return extend_linearly(self.stage_shift, clamp_load)

    def clamp_temperature(self, stage_temperature: NDArray[np.float64], clamp_load: float) -> float:
        """Returns `r_eq = F_c * T1 + (1 - F_c) * T2 + d_c` in K."""
        fraction = self.profile_fraction
        return (
            fraction * stage_temperature[0]
            + (1.0 - fraction) * stage_temperature[1]
            + extend_linearly(self.clamp_increment, clamp_load)
        )

    def regenerator_2_profile(self, stage_temperature: NDArray[np.float64], clamp_load: float) -> NDArray[np.float64]:
        """Returns regenerator 2's target temperatures in K at its quadrature positions."""
        stage_1, stage_2 = stage_temperature
        if clamp_load == 0.0:
            return self.unloaded_fraction * stage_1 + (1.0 - self.unloaded_fraction) * stage_2
        loads = self.loaded_shape.x
        shape = self.loaded_shape(np.clip(clamp_load, loads[0], loads[-1]))
        clamp = self.clamp_temperature(stage_temperature, clamp_load)
        return np.where(self.above_clamp, clamp + (stage_1 - clamp) * shape, stage_2 + (clamp - stage_2) * shape)


@dataclass(frozen=True, slots=True, kw_only=True)
class PulseTubeHead:
    """A pulse-tube head linking its stages, its clamp, and its regenerator nodes.

    Attributes:
        name: Unique link identifier within the graph.
        stages: Stage-1 and stage-2 node names.
        clamp: Clamp node name, or None on the test stand.
        contact_resistance: Resistance `R_contact` in K/W between the clamp node and the head-side clamp, or None on
            the test stand.
        coefficients: The fitted coefficients.
        warm_flange_temperature: Temperature `T0` in K at the warm end of regenerator 1.
        energy_scale: Normalization `E_scale` in J of excess regenerator-1 energy.
        cold_map: The zero-clamp-load cooling power.
        clamp_response: The clamp-load response.
        regenerator_1: Regenerator 1.
        regenerator_2: Regenerator 2.
        memory_node: Index in regenerator 2 of the node just below the clamp, which carries the clamp's thermal memory.
    """

    name: str
    stages: tuple[str, str]
    clamp: str | None
    contact_resistance: float | None
    coefficients: HeadCoefficients
    warm_flange_temperature: float
    energy_scale: float
    cold_map: ColdMap
    clamp_response: ClampResponse
    regenerator_1: Regenerator
    regenerator_2: Regenerator
    memory_node: int

    @property
    def internal_nodes(self) -> tuple[Node, ...]:
        """Regenerator nodes from warm to cold."""
        return (*self.regenerator_1.nodes, *self.regenerator_2.nodes)

    @property
    def nodes(self) -> tuple[str, ...]:
        """Node names: the stages, the clamp, if any, and the regenerator nodes."""
        clamp = () if self.clamp is None else (self.clamp,)
        return (*self.stages, *clamp, *(node.name for node in self.internal_nodes))

    def solve_clamp_load(
        self,
        stage_temperature: NDArray[np.float64],
        clamp_temperature: float,
        regenerator_2_temperature: NDArray[np.float64],
    ) -> float:
        """Solves `r_eq * theta_h / theta_eq,h + R_contact * load = TR` for the clamp load in W."""
        response, memory = self.clamp_response, self.memory_node

        def residual(clamp_load: float) -> float:
            _, target = self.regenerator_2.target(response.regenerator_2_profile(stage_temperature, clamp_load))
            head_side = (
                response.clamp_temperature(stage_temperature, clamp_load)
                * regenerator_2_temperature[memory]
                / target[memory]
            )
            return head_side + self.contact_resistance * clamp_load - clamp_temperature

        increment = response.clamp_increment
        lower = -response.clamp_temperature(stage_temperature, 0.0) / increment(increment.x[0], 1)
        upper = increment.x[-1]
        while residual(upper) < 0.0:
            upper *= 2.0
        return brentq(residual, lower, upper)

    def calc_energy_balance_terms(self, temperature: NDArray[np.float64]) -> EnergyBalanceTerms:
        """Calculates the stage and clamp heat flows and the regenerator nodes' energy rates.

        Args:
            temperature: Temperatures in K in the order of `nodes`.

        Returns:
            Heat flows in W in the order of `nodes`, positive into each node.
        """
        regenerator_1, regenerator_2 = self.regenerator_1, self.regenerator_2
        coefficients, warm = self.coefficients, self.warm_flange_temperature
        count_1, count_2 = len(regenerator_1.nodes), len(regenerator_2.nodes)
        stage_temperature = temperature[:2]
        regenerator_1_temperature = temperature[-count_1 - count_2 : -count_2]
        regenerator_2_temperature = temperature[-count_2:]
        if self.clamp is None:
            clamp_load = 0.0
        else:
            clamp_load = self.solve_clamp_load(stage_temperature, temperature[2], regenerator_2_temperature)

        positions = regenerator_1.positions
        target_energy_1, target_temperature_1 = regenerator_1.target(
            (1.0 - positions) * warm + positions * stage_temperature[0]
        )
        _, target_temperature_2 = regenerator_2.target(
            self.clamp_response.regenerator_2_profile(stage_temperature, clamp_load)
        )
        rate_1 = (
            regenerator_1.mass
            * coefficients.regenerator_1_exchange_coefficient
            * (target_temperature_1 - regenerator_1_temperature)
        )
        rate_2 = (
            regenerator_2.mass
            * coefficients.regenerator_2_exchange_coefficient
            * (target_temperature_2 - regenerator_2_temperature)
        )

        excess_energy = np.sum(
            regenerator_1.mass * (regenerator_1.specific_energy(regenerator_1_temperature) - target_energy_1)
        )
        multiplier = np.exp(-coefficients.energy_coupling_coefficient * excess_energy / self.energy_scale)

        shifted = stage_temperature - self.clamp_response.stage_shifts(clamp_load)
        amplitude = np.array([coefficients.stage_1_warm_amplitude, coefficients.stage_2_warm_amplitude])
        exponent = np.array([coefficients.stage_1_warm_exponent, coefficients.stage_2_warm_exponent])
        equilibrated = self.cold_map.cooling_power(shifted, amplitude * (shifted / warm) ** exponent)
        storage = np.array(
            [
                regenerator_1.midpoint @ rate_1 + (1.0 - regenerator_2.midpoint) @ rate_2,
                regenerator_2.midpoint @ rate_2,
            ]
        )
        cooling_power = np.array([multiplier, 1.0]) * equilibrated + storage

        clamp_heat_flow = () if self.clamp is None else (-clamp_load,)
        return EnergyBalanceTerms(heat_flow=np.concatenate((-cooling_power, clamp_heat_flow, rate_1, rate_2)))


def load_pulse_tube_head(
    config: Mapping[str, Any], *, load_graph: Callable[[Mapping[str, Any]], Graph]
) -> PulseTubeHead:
    """Loads a calibrated head from a link table.

    Args:
        config: `name`, `data` (a calibration file path relative to the working directory), `stages` in stage order,
            and `warm_flange_temperature` in K. A mounted head also has a `clamp` node and a `contact_resistance` in
            K/W.
        load_graph: Builds the calibration's test stand from its `fridge` table.

    Raises:
        ConfigurationError: Only one of `clamp` and `contact_resistance` is given, or the resistance is not positive.
    """
    clamp = config.get("clamp")
    resistance = config.get("contact_resistance")
    if (clamp is None) != (resistance is None) or (resistance is not None and resistance <= 0):
        raise ConfigurationError("A mounted head needs both a clamp node and a positive contact_resistance.")
    calibration = read_head_calibration(config["data"])
    first, second = config["stages"]
    return build_pulse_tube_head(
        calibration,
        load_graph(calibration.test_stand),
        name=config["name"],
        stages=(first, second),
        warm_flange_temperature=float(config["warm_flange_temperature"]),
        clamp=clamp,
        contact_resistance=resistance,
    )


def build_pulse_tube_head(
    calibration: HeadCalibration,
    test_stand: Graph,
    /,
    *,
    name: str,
    stages: tuple[str, str],
    warm_flange_temperature: float,
    clamp: str | None = None,
    contact_resistance: float | None = None,
    shared: SharedParameters | None = None,
) -> PulseTubeHead:
    """Builds a head from its calibration.

    Args:
        calibration: The head's calibration.
        test_stand: The calibration's test stand, built from `calibration.test_stand`.
        name: Unique link identifier within the graph.
        stages: Stage-1 and stage-2 node names in the graph the head joins.
        warm_flange_temperature: Temperature at the warm end of regenerator 1, in K.
        clamp: Clamp node name, or None on the test stand.
        contact_resistance: Resistance between the clamp node and the head-side clamp, in K/W.
        shared: Parameters shared by every head. Defaults to the packaged parameters.
    """
    shared = shared or load_shared_parameters()
    design = calibration.design
    node_names = (f"{name} regenerator node {number}" for number in itertools.count(1))
    table_temperature = np.geomspace(*shared.energy_temperature_range, shared.energy_knots)
    quadrature = np.polynomial.legendre.leggauss(shared.quadrature_order)
    regenerator_1 = build_regenerator(design.regenerator_1, node_names, table_temperature, quadrature)
    regenerator_2 = build_regenerator(
        design.regenerator_2, node_names, table_temperature, quadrature, split_at=design.clamp_fraction
    )

    temperature, heater_power = calibration.load_map.temperature, calibration.load_map.heater_power
    minimum, maximum = temperature.min(axis=0), temperature.max(axis=0)
    span = maximum - minimum
    cold_map = ColdMap(
        minimum_temperature=minimum,
        maximum_temperature=maximum,
        tail_length=span * np.asarray(shared.cold_tail_span_multiplier),
        blend_width=shared.cold_blend_width,
        heater_power=RBFInterpolator(
            (temperature - minimum) / span, heater_power, kernel="thin_plate_spline", degree=1
        ),
        test_stand=test_stand,
        test_stand_stages=np.array([test_stand.get_node_index(stage) for stage in calibration.test_stand_stages]),
        test_stand_temperature=calibration.warm_flange_temperature,
    )

    return PulseTubeHead(
        name=name,
        stages=stages,
        clamp=clamp,
        contact_resistance=contact_resistance,
        coefficients=calibration.coefficients,
        warm_flange_temperature=warm_flange_temperature,
        energy_scale=shared.energy_scale,
        cold_map=cold_map,
        clamp_response=build_clamp_response(shared, design.clamp_fraction, regenerator_2.positions),
        regenerator_1=regenerator_1,
        regenerator_2=regenerator_2,
        memory_node=int(np.argmax(regenerator_2.midpoint > design.clamp_fraction)),
    )


def build_regenerator(
    layers: tuple[RegeneratorLayer, ...],
    node_names: Iterator[str],
    table_temperature: NDArray[np.float64],
    quadrature: tuple[NDArray[np.float64], NDArray[np.float64]],
    *,
    split_at: float | None = None,
) -> Regenerator:
    """Builds a regenerator with one node per layer, splitting the layer that spans `split_at`.

    Args:
        layers: Layers from warm to cold.
        node_names: Supplies each node's name.
        table_temperature: Temperatures in K at which to tabulate specific energies.
        quadrature: Gauss-Legendre abscissas and weights on [-1, 1].
        split_at: Position, as a fraction of the regenerator length, at which to split a layer.
    """
    abscissa, weights = quadrature
    nodes, energies, spans = [], [], []
    start = 0.0
    for layer in layers:
        end = start + layer.length_fraction
        edges = [start, split_at, end] if split_at is not None and start < split_at < end else [start, end]
        for low, high in itertools.pairwise(edges):
            mass = float(layer.mass * (high - low) / (end - start))
            nodes.append(Node(next(node_names), ((layer.material, mass),)))
            energies.append(layer.material.specific_heat_integral(table_temperature))
            spans.append((low, high))
        start = end

    low, high = np.array(spans).T
    return Regenerator(
        nodes=tuple(nodes),
        table_temperature=table_temperature,
        table_energy=np.array(energies),
        mass=np.array([mass for node in nodes for _, mass in node.composition]),
        midpoint=(low + high) / 2.0,
        positions=low[:, None] + (abscissa + 1.0) * (high - low)[:, None] / 2.0,
        weights=weights / 2.0,
    )


def build_clamp_response(
    shared: SharedParameters, clamp_fraction: float, positions: NDArray[np.float64]
) -> ClampResponse:
    """Builds the clamp response from the published tables.

    Args:
        shared: Parameters shared by every head.
        clamp_fraction: Clamp position `x_R` as a fraction of regenerator 2's length.
        positions: Regenerator 2's quadrature positions.
    """
    profile = shared.unloaded_profile
    warm, cold = profile.reference_stage_temperature
    unloaded = PchipInterpolator(profile.position, (profile.temperature - cold) / (warm - cold))
    fraction = float(unloaded(clamp_fraction))

    table = shared.clamp_response
    load = table.load
    shift = np.column_stack(
        (
            table.stage_1_temperature - table.stage_1_temperature[0],
            table.stage_2_temperature - table.stage_2_temperature[0],
        )
    )
    increment = (
        table.clamp_temperature - table.clamp_temperature[0] - fraction * shift[:, 0] - (1.0 - fraction) * shift[:, 1]
    )
    increment_slope = PchipInterpolator(load, increment).derivative()(load)
    increment_slope[[0, -1]] = np.diff(increment)[[0, -1]] / np.diff(load)[[0, -1]]

    loaded = shared.loaded_profiles
    source_clamp = loaded.source_clamp_fraction
    source_warm, source_cold = loaded.reference_stage_temperature
    above = positions <= clamp_fraction
    source_position = np.where(
        above,
        source_clamp * positions / clamp_fraction,
        source_clamp + (1.0 - source_clamp) * (positions - clamp_fraction) / (1.0 - clamp_fraction),
    )
    unloaded_fraction = unloaded(positions)
    shapes = [np.where(above, (unloaded_fraction - fraction) / (1.0 - fraction), unloaded_fraction / fraction)]
    for temperature in loaded.temperature[1:]:
        at_clamp = np.interp(source_clamp, loaded.position, temperature)
        local = np.interp(source_position, loaded.position, temperature)
        shape = np.where(
            above, (local - at_clamp) / (source_warm - at_clamp), (local - source_cold) / (at_clamp - source_cold)
        )
        shapes.append(np.clip(shape, 0.0, 1.0))
    shape_slope = PchipInterpolator(loaded.load, np.array(shapes), axis=0).derivative()(loaded.load)
    shape_slope[[0, -1]] = 0.0

    return ClampResponse(
        stage_shift=PchipInterpolator(load, shift, axis=0),
        clamp_increment=CubicHermiteSpline(load, increment, increment_slope),
        profile_fraction=fraction,
        unloaded_fraction=unloaded_fraction,
        above_clamp=above,
        loaded_shape=CubicHermiteSpline(loaded.load, np.array(shapes), shape_slope, axis=0),
    )


def interpolate(
    value: float | NDArray[np.float64], x: NDArray[np.float64], y: NDArray[np.float64]
) -> np.float64 | NDArray[np.float64]:
    """Interpolates `y` against `x` linearly, continuing the end intervals' slopes beyond them."""
    below = np.minimum(value - x[0], 0.0)
    above = np.maximum(value - x[-1], 0.0)
    return np.interp(value, x, y) + below * (y[1] - y[0]) / (x[1] - x[0]) + above * (y[-1] - y[-2]) / (x[-1] - x[-2])


def extend_linearly(curve: PchipInterpolator | CubicHermiteSpline, x: float) -> NDArray[np.float64]:
    """Evaluates a curve, continued linearly along its endpoint tangents beyond its knots."""
    end = np.clip(x, curve.x[0], curve.x[-1])
    return curve(end) + (x - end) * curve(end, 1)
