"""Five-coefficient candidate, with matrix energies and a floating clamp.

Equations follow docs/proposals/pulse-tube-model-direction.md, 2026-11-10.
This experiment uses the existing head inventories without the retired head's
packing, valve, waveform, breakthrough, or stage-coupling coefficients.
"""

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from numpy.typing import NDArray
from scipy.interpolate import PchipInterpolator
from scipy.optimize import least_squares, root
from scipy.special import expit, logit

from pulse_tube_cooldown.library.heads.load import load_head_model
from pulse_tube_cooldown.library.heads.model import HeadModel
from simulations.exploration.pulse_tube_candidate.properties import (
    HeliumTable,
    ThermalProperty,
    heat_capacity_table,
    mean_integral,
    temperature_knots,
)


class PropertyLimitError(ValueError):
    """The thermal state lies outside the helium table's supported interval."""


class ProfileSolveError(RuntimeError):
    """The energy/profile or equilibrium solve did not meet its tolerance."""


@dataclass(frozen=True, slots=True)
class Parameters:
    """W_star [W], allocation a and b, and exchange G1 and G2 [W/K]."""

    power: float
    allocation: float
    branch_floor: float
    exchange_1: float
    exchange_2: float

    def as_array(self) -> NDArray:
        """Returns coefficients in the experiment configuration's order."""
        return np.array([self.power, self.allocation, self.branch_floor, self.exchange_1, self.exchange_2])


@dataclass(frozen=True, slots=True)
class Piece:
    """One uniform matrix inventory with a heat-capacity interpolant."""

    start: float
    end: float
    mass: float
    pore_volume: float
    heat: ThermalProperty
    alpha_beta: ThermalProperty
    alpha_cp: PchipInterpolator
    regenerator: int


@dataclass(frozen=True, slots=True)
class Panel:
    """A profile connection's intersection with one material piece."""

    piece: int
    connection: int
    left: NDArray
    right: NDArray
    length: float


@dataclass(frozen=True, slots=True)
class Evaluation:
    """Flows, matrix energy rates, and source diagnostics at one complete state."""

    energy_rate: NDArray
    heat_flow: NDArray
    flows: NDArray
    conductance: NDArray
    pressure_fraction: float
    demand: float
    pass_demand: float
    fill_demand: float
    profile_measure: float
    transmission: NDArray
    input_power_cap: float
    tube_power: NDArray
    wall_leak: NDArray


class Candidate:
    """Fixed geometry, material tables, and quadrature for one head design."""

    def __init__(self, design: HeadModel, helium: HeliumTable, *, source: dict, numerics: dict):
        """Prepares the specified material pieces and shared-grid interpolants."""
        self.design = design
        self.helium = helium
        self.source = source
        self.numerics = numerics
        self.temperature_min = float(helium.temperature[0])
        self.temperature_max = float(helium.temperature[-1])
        self.reference = numerics["energy_reference_temperature"]
        self.gauss_x, self.gauss_w = leggauss(numerics["quadrature_order"])
        self.flux_breaks = np.union1d(
            numerics["flux_temperature_breaks"],
            self.temperature_min
            * numerics["flux_temperature_ratio"]
            ** np.arange(
                1
                + int(np.log(self.temperature_max / self.temperature_min) / np.log(numerics["flux_temperature_ratio"]))
            ),
        )
        pieces = []
        self.positions = []
        self.node_codes = []
        self.connection_offsets = [0]
        for regenerator, section in enumerate(design.regenerators):
            position = 0.0
            centers = []
            for layer in section.layers:
                end = position + layer.length_fraction * section.length
                cuts = [position, end]
                clamp = design.clamp_fraction * section.length
                if regenerator == 1 and position < clamp < end:
                    cuts.insert(1, clamp)
                heat = heat_capacity_table(layer.material, helium.temperature)
                for left, right in zip(cuts[:-1], cuts[1:], strict=True):
                    fraction = (right - left) / (end - position)
                    mass = layer.mass * fraction
                    pore_volume = layer.pore_volume * fraction
                    matrix_capacity = mass * heat.value(helium.temperature)
                    gas_capacity = (
                        pore_volume * helium.density.value(helium.temperature) * helium.cp.value(helium.temperature)
                    )
                    alpha = matrix_capacity / (matrix_capacity + gas_capacity)
                    pieces.append(
                        Piece(
                            left,
                            right,
                            mass,
                            pore_volume,
                            heat,
                            ThermalProperty.from_values(
                                helium.temperature, alpha * helium.beta.value(helium.temperature)
                            ),
                            PchipInterpolator(helium.temperature, alpha * helium.cp.value(helium.temperature)),
                            regenerator,
                        )
                    )
                    centers.append(((left + right) / 2, len(pieces) - 1))
                position = end
            nodes = [(0.0, 8 if regenerator == 0 else 9), *centers, (section.length, 9 if regenerator == 0 else 10)]
            if regenerator == 1:
                nodes.append((design.clamp_fraction * section.length, 7))
            nodes.sort()
            self.positions.append(np.array([p for p, _ in nodes]))
            self.node_codes.append(np.array([code for _, code in nodes]))
            self.connection_offsets.append(self.connection_offsets[-1] + len(nodes) - 1)
        self.pieces = tuple(pieces)
        if len(pieces) != 7:
            raise ValueError("This experiment requires the seven material pieces in the proposal.")
        panels = []
        for i, _section in enumerate(design.regenerators):
            for connection, (left, right) in enumerate(zip(self.positions[i][:-1], self.positions[i][1:], strict=True)):
                for piece_index, piece in enumerate(pieces):
                    if piece.regenerator != i:
                        continue
                    start, end = max(left, piece.start), min(right, piece.end)
                    if start >= end:
                        continue
                    rows = []
                    for x in (start, end):
                        fraction = (x - left) / (right - left)
                        row = np.zeros(11)
                        row[self.node_codes[i][connection]] = 1 - fraction
                        row[self.node_codes[i][connection + 1]] = fraction
                        rows.append(row)
                    panels.append(Panel(piece_index, self.connection_offsets[i] + connection, *rows, end - start))
        self.panels = tuple(panels)
        self.panel_piece = np.array([panel.piece for panel in panels])
        self.panel_connection = np.array([panel.connection for panel in panels])
        self.panel_regenerator = np.array([pieces[panel.piece].regenerator for panel in panels])
        self.left = np.stack([panel.left for panel in panels])
        self.right = np.stack([panel.right for panel in panels])
        self.panel_lengths = np.array([panel.length for panel in panels])
        self.panel_mass = np.array(
            [pieces[p.piece].mass * p.length / (pieces[p.piece].end - pieces[p.piece].start) for p in panels]
        )
        self.panel_pore_volume = np.array(
            [pieces[p.piece].pore_volume * p.length / (pieces[p.piece].end - pieces[p.piece].start) for p in panels]
        )
        self.heat_primitive = np.stack([piece.heat.primitive.c for piece in pieces])
        self.heat_value = np.stack([piece.heat.value.c for piece in pieces])
        self.heat_second_primitive = np.stack([piece.heat.second_primitive.c for piece in pieces])
        self.heat_reference = np.array([piece.heat.primitive(self.reference) for piece in pieces])
        self.alpha_beta_value = np.stack([piece.alpha_beta.value.c for piece in pieces])
        self.alpha_beta_primitive = np.stack([piece.alpha_beta.primitive.c for piece in pieces])
        self.alpha_cp_value = np.stack([piece.alpha_cp.c for piece in pieces])
        self.wall_conductivity = np.stack(
            [
                PchipInterpolator(helium.temperature, section.wall_material.thermal_conductivity(helium.temperature)).c
                for section in design.regenerators
            ]
        )
        self.heat_inverse = tuple(
            PchipInterpolator(
                piece.heat.primitive(helium.temperature) - piece.heat.primitive(self.reference),
                helium.temperature,
                extrapolate=True,
            )
            for piece in pieces
        )

    def interpolate(self, coefficients: NDArray, rows: NDArray, temperature: NDArray) -> NDArray:
        """Evaluates shared-grid piecewise polynomials for different material rows."""
        indices = np.clip(
            np.searchsorted(self.helium.temperature, temperature, side="right") - 1, 0, len(self.helium.temperature) - 2
        )
        delta = temperature - self.helium.temperature[indices]
        result = coefficients[rows, 0, indices]
        for degree in range(1, coefficients.shape[1]):
            result = result * delta + coefficients[rows, degree, indices]
        return result

    def check_temperature(self, temperature: NDArray) -> None:
        """Rejects property evaluations outside the tabulated interval."""
        if (
            np.any(~np.isfinite(temperature))
            or np.any(temperature < self.temperature_min)
            or np.any(temperature > self.temperature_max)
        ):
            raise PropertyLimitError(
                f"Temperature range [{np.min(temperature):.8g}, {np.max(temperature):.8g}] K; "
                f"supported [{self.temperature_min:g}, {self.temperature_max:g}] K."
            )

    @staticmethod
    def complete_temperature(centers: NDArray, boundary: NDArray) -> NDArray:
        """Orders seven centers, clamp, warm flange, stage 1, and stage 2."""
        return np.concatenate((centers, boundary))

    def energies(self, centers: NDArray, boundary: NDArray) -> NDArray:
        """Integrates matrix energy on each material's physical extent."""
        temperature = self.complete_temperature(centers, boundary)
        ta, tb = self.left @ temperature, self.right @ temperature
        delta = tb - ta
        close = np.abs(delta) < 0.01
        average = np.empty_like(delta)
        rows = self.panel_piece[~close]
        average[~close] = (
            self.interpolate(self.heat_second_primitive, rows, tb[~close])
            - self.interpolate(self.heat_second_primitive, rows, ta[~close])
        ) / delta[~close] - self.heat_reference[rows]
        t = (ta[close, None] + tb[close, None]) / 2 + self.gauss_x * delta[close, None] / 2
        rows = np.broadcast_to(self.panel_piece[close, None], t.shape)
        values = self.interpolate(self.heat_primitive, rows, t) - self.heat_reference[rows]
        average[close] = values @ self.gauss_w / 2
        return np.bincount(self.panel_piece, self.panel_mass * average, minlength=7)

    def initial_centers(self, boundary: NDArray) -> NDArray:
        """Uses a linear initial profile through each regenerator's endpoints."""
        centers = np.empty(8)
        for j, piece in enumerate(self.pieces):
            i = piece.regenerator
            fraction = (piece.start + piece.end) / (2 * self.design.regenerators[i].length)
            centers[j] = boundary[i] * (1 - fraction) + boundary[i + 1] * fraction
        centers[7] = boundary[1] * (1 - self.design.clamp_fraction) + boundary[2] * self.design.clamp_fraction
        return centers

    def energy_jacobian(self, centers: NDArray, boundary: NDArray) -> NDArray:
        """Differentiates each integrated material energy with respect to profile knots."""
        temperature = self.complete_temperature(centers, boundary)
        ta, tb = self.left @ temperature, self.right @ temperature
        delta = tb - ta
        close = np.abs(delta) < 0.1
        left_derivative = np.empty_like(delta)
        right_derivative = np.empty_like(delta)
        rows = self.panel_piece[~close]
        average = (
            self.interpolate(self.heat_second_primitive, rows, tb[~close])
            - self.interpolate(self.heat_second_primitive, rows, ta[~close])
        ) / delta[~close]
        left_derivative[~close] = (average - self.interpolate(self.heat_primitive, rows, ta[~close])) / delta[~close]
        right_derivative[~close] = (self.interpolate(self.heat_primitive, rows, tb[~close]) - average) / delta[~close]
        t = (ta[close, None] + tb[close, None]) / 2 + self.gauss_x * delta[close, None] / 2
        rows = np.broadcast_to(self.panel_piece[close, None], t.shape)
        cp = self.interpolate(self.heat_value, rows, t)
        left_derivative[close] = cp @ (self.gauss_w * (1 - self.gauss_x) / 4)
        right_derivative[close] = cp @ (self.gauss_w * (1 + self.gauss_x) / 4)
        result = np.zeros((7, 11))
        np.add.at(
            result,
            self.panel_piece,
            self.panel_mass[:, None] * (left_derivative[:, None] * self.left + right_derivative[:, None] * self.right),
        )
        return result

    def evaluate(
        self,
        centers: NDArray,
        boundary: NDArray,
        parameters: Parameters,
        *,
        compressors: int,
        frequency: float,
        clamp_temperature: float | None = None,
        clamp_conductance: float = 0.0,
        other_demand: float = 0.0,
    ) -> Evaluation:
        """Calculates conservative flows and the common pressure-limited source."""
        temperature = self.complete_temperature(centers, boundary)
        self.check_temperature(temperature)
        ta, tb = self.left @ temperature, self.right @ temperature
        primitive_left = self.interpolate(self.alpha_beta_primitive, self.panel_piece, ta)
        logarithmic_increment = self.interpolate(self.alpha_beta_primitive, self.panel_piece, tb) - primitive_left
        regenerator = self.panel_regenerator
        prefix = np.zeros(len(self.panels))
        transmission = np.empty(2)
        for i in range(2):
            selected = np.flatnonzero(regenerator == i)
            increments = logarithmic_increment[selected]
            prefix[selected] = np.r_[0.0, np.cumsum(increments)[:-1]]
            transmission[i] = np.exp(np.sum(increments))
        viscosity_average = mean_integral(self.helium.kinematic_viscosity, ta, tb)
        selected = regenerator == 1
        profile_measure = float(
            np.sum(self.panel_lengths[selected] * viscosity_average[selected])
            / self.design.regenerators[1].length
            / self.helium.kinematic_viscosity.value(boundary[0])
        )
        branch_2 = parameters.power / (parameters.branch_floor + parameters.allocation * profile_measure)
        junction = parameters.power + branch_2
        input_power = junction / transmission[0]
        power_scale = np.array([input_power, branch_2])
        response_average = mean_integral(self.helium.mass_response, ta, tb)
        mass_response = self.panel_pore_volume @ response_average
        for i, tube in enumerate(self.design.pulse_tubes):
            response = mean_integral(self.helium.mass_response, np.array([boundary[i + 1]]), np.array([boundary[0]]))[0]
            mass_response += tube.volume * response
            mass_response += self.design.heat_exchanger_voids[i] * self.helium.mass_response.value(boundary[i + 1])
        mass_response += self.design.dead_volume * self.helium.mass_response.value(boundary[0])
        pressure = self.source["pressure_ceiling"]
        pass_demand = 2 * self.helium.density.value(boundary[0]) * input_power / (np.pi * pressure)
        fill_demand = 2 * frequency * pressure * mass_response
        demand = pass_demand + fill_demand
        z = min(1.0, compressors * self.source["delivery_per_compressor"] / (demand + other_demand))
        exchange = np.array([parameters.exchange_1, parameters.exchange_2])
        q_temperatures, q_weights, q_panels = [], [], []
        for index, panel in enumerate(self.panels):
            left, right = ta[index], tb[index]
            cuts_array = np.array([0.0, 1.0])
            if right != left:
                fractions = (self.flux_breaks - left) / (right - left)
                cuts_array = np.sort(np.r_[cuts_array, fractions[(0 < fractions) & (fractions < 1)]])
            fraction = (
                (cuts_array[:-1, None] + cuts_array[1:, None]) / 2
                + (cuts_array[1:, None] - cuts_array[:-1, None]) / 2 * self.gauss_x
            ).ravel()
            weights = ((cuts_array[1:, None] - cuts_array[:-1, None]) / 2 * self.gauss_w).ravel() * panel.length
            q_temperatures.append(left + fraction * (right - left))
            q_weights.append(weights)
            q_panels.append(np.full(len(weights), index, dtype=int))
        t, weights, panel_indices = np.concatenate(q_temperatures), np.concatenate(q_weights), np.concatenate(q_panels)
        pieces = self.panel_piece[panel_indices]
        regenerators = regenerator[panel_indices]
        lengths = np.array([section.length for section in self.design.regenerators])
        wall_areas = np.array([section.wall_area for section in self.design.regenerators])
        local_power_cap = power_scale[regenerators] * np.exp(
            prefix[panel_indices]
            + self.interpolate(self.alpha_beta_primitive, pieces, t)
            - primitive_left[panel_indices]
        )
        rho = self.helium.density.value(t)
        alpha_cp = self.interpolate(self.alpha_cp_value, pieces, t)
        dgrad = (
            2
            * z**2
            * (alpha_cp * rho * local_power_cap / pressure) ** 2
            / (exchange[regenerators] / lengths[regenerators])
        )
        conductivity = dgrad + wall_areas[regenerators] * self.interpolate(self.wall_conductivity, regenerators, t)
        phi = (1 - self.interpolate(self.alpha_beta_value, pieces, t) * t) * z**2 * local_power_cap
        connections = self.panel_connection[panel_indices]
        resistance = np.bincount(connections, weights / conductivity, minlength=self.connection_offsets[-1])
        forcing = np.bincount(connections, weights * phi / conductivity, minlength=self.connection_offsets[-1])
        difference = np.concatenate([temperature[codes[:-1]] - temperature[codes[1:]] for codes in self.node_codes])
        flows = (difference + forcing) / resistance
        # R1 has three connections; R2 has seven, with the clamp at connections 3 and 4.
        energy_rate = np.array(
            [
                flows[0] - flows[1],
                flows[1] - flows[2],
                flows[3] - flows[4],
                flows[4] - flows[5],
                flows[5] - flows[6],
                flows[7] - flows[8],
                flows[8] - flows[9],
            ]
        )
        qr = 0.0 if clamp_temperature is None else clamp_conductance * (centers[7] - clamp_temperature)
        wall_leak = np.array(
            [
                tube.wall_area
                / tube.length
                * (
                    tube.wall_material.conductivity_integral(boundary[0])
                    - tube.wall_material.conductivity_integral(boundary[i + 1])
                )
                for i, tube in enumerate(self.design.pulse_tubes)
            ]
        )
        tube_power = z**2 * np.array([parameters.power, branch_2 * transmission[1]])
        q = np.array(
            [
                flows[2] - flows[3] - tube_power[0] + wall_leak[0],
                qr,
                flows[9] - tube_power[1] + wall_leak[1],
            ]
        )
        return Evaluation(
            energy_rate,
            q,
            flows,
            1 / resistance,
            float(z),
            float(demand),
            float(pass_demand),
            float(fill_demand),
            profile_measure,
            transmission,
            float(input_power),
            tube_power,
            wall_leak,
        )

    def _encode(self, temperature: NDArray) -> NDArray:
        fraction = (temperature - self.temperature_min) / (self.temperature_max - self.temperature_min)
        return logit(np.clip(fraction, 1e-12, 1 - 1e-12))

    def _decode(self, coordinate: NDArray) -> NDArray:
        return self.temperature_min + (self.temperature_max - self.temperature_min) * expit(coordinate)

    def recover_profile(
        self,
        energy: NDArray,
        boundary: NDArray,
        parameters: Parameters,
        *,
        compressors: int,
        frequency: float,
        guess: NDArray | None = None,
    ) -> tuple[NDArray, Evaluation]:
        """Inverts the seven energies with a massless, unloaded clamp balance."""
        self.check_temperature(boundary)
        if guess is None:
            guess = np.array(
                [
                    inverse(e / piece.mass)
                    for inverse, e, piece in zip(self.heat_inverse, energy, self.pieces, strict=True)
                ]
                + [self.initial_centers(boundary)[7]]
            )
        capacity = np.array(
            [
                piece.mass * piece.heat.value(float(np.clip(t, self.temperature_min, self.temperature_max)))
                for piece, t in zip(self.pieces, guess[:7], strict=True)
            ]
        )

        cached_coordinate, cached_centers, cached_evaluation = None, None, None

        def at(coordinate: NDArray) -> tuple[NDArray, Evaluation]:
            nonlocal cached_coordinate, cached_centers, cached_evaluation
            if cached_coordinate is None or not np.array_equal(coordinate, cached_coordinate):
                cached_coordinate = coordinate.copy()
                cached_centers = self._decode(coordinate)
                cached_evaluation = self.evaluate(
                    cached_centers, boundary, parameters, compressors=compressors, frequency=frequency
                )
            return cached_centers, cached_evaluation

        def residual(coordinate: NDArray) -> NDArray:
            centers, evaluation = at(coordinate)
            clamp_scale = evaluation.conductance[6] + evaluation.conductance[7]
            return np.r_[
                (self.energies(centers, boundary) - energy) / capacity,
                (evaluation.flows[6] - evaluation.flows[7]) / clamp_scale,
            ]

        def jacobian(coordinate: NDArray) -> NDArray:
            centers, evaluation = at(coordinate)
            result = np.zeros((8, 8))
            result[:7] = self.energy_jacobian(centers, boundary)[:, :8] / capacity[:, None]
            clamp_scale = evaluation.conductance[6] + evaluation.conductance[7]
            result[7, 4] = evaluation.conductance[6] / clamp_scale
            result[7, 5] = evaluation.conductance[7] / clamp_scale
            result[7, 7] = -1
            fraction = expit(coordinate)
            return result * ((self.temperature_max - self.temperature_min) * fraction * (1 - fraction))

        solution = root(residual, self._encode(guess), jac=jacobian, tol=1e-9)
        centers = self._decode(solution.x)
        error = float(np.max(np.abs(residual(solution.x))))
        if error > self.numerics["inversion_temperature_tolerance"]:
            solution = least_squares(residual, solution.x, xtol=1e-12, ftol=1e-12, gtol=1e-12, max_nfev=60)
            centers = self._decode(solution.x)
            error = float(np.max(np.abs(residual(solution.x))))
        if error > self.numerics["inversion_temperature_tolerance"]:
            raise ProfileSolveError(f"Energy/profile inversion residual {error:.6g} K.")
        return centers, self.evaluate(centers, boundary, parameters, compressors=compressors, frequency=frequency)

    def equilibrate(
        self,
        boundary: NDArray,
        parameters: Parameters,
        *,
        compressors: int,
        frequency: float,
        guess: NDArray | None = None,
    ) -> tuple[NDArray, Evaluation]:
        """Relaxes all matrix energies and the unloaded clamp at held stages."""
        self.check_temperature(boundary)
        initial = self.initial_centers(boundary) if guess is None else guess

        def residual(coordinate: NDArray) -> NDArray:
            evaluation = self.evaluate(
                self._decode(coordinate), boundary, parameters, compressors=compressors, frequency=frequency
            )
            return np.r_[evaluation.energy_rate, evaluation.flows[6] - evaluation.flows[7]]

        solution = root(residual, self._encode(initial), tol=1e-9)
        error = float(np.max(np.abs(residual(solution.x))))
        if error > self.numerics["equilibrium_flow_tolerance"]:
            solution = least_squares(residual, solution.x, xtol=1e-12, ftol=1e-12, gtol=1e-12, max_nfev=150)
            error = float(np.max(np.abs(residual(solution.x))))
        if error > self.numerics["equilibrium_flow_tolerance"]:
            raise ProfileSolveError(f"Held-stage equilibrium residual {error:.6g} W.")
        centers = self._decode(solution.x)
        return centers, self.evaluate(centers, boundary, parameters, compressors=compressors, frequency=frequency)


def build_candidate(config: dict) -> Candidate:
    """Builds the experiment using the existing PT1.5 design and property data."""
    root = Path(__file__).resolve().parents[3]
    grid = temperature_knots(
        config["numerics"], material_directory=root / "src/pulse_tube_cooldown/library/materials/data"
    )
    helium = HeliumTable.build(grid, pressure=config["source"]["mean_pressure"])
    return Candidate(load_head_model("PT1.5"), helium, source=config["source"], numerics=config["numerics"])
