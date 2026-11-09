"""Fits SN09's five head parameters to its May 2026 test-stand load map.

At each load-map point the warm flange and both stages sit at their measured temperatures and each regenerator column
is solved to steady state with no load on the clamp. The residual at each stage is the heat flow into it from the head
and the stand's tube walls plus its heater power, zero where the parameters reproduce the point. Running this file
writes a candidate parameter file, diagnostics, and a figure beside this module.
The fit searches C_2 / a_2 and eta / a_2 on a selected equilibrium branch, with a nonnegative
linear fit for a_1, a_2, and C_1 at each pair of ratios.
Heat-flow residuals are in W.
"""

import argparse
import hashlib
import itertools
import json
import time
import tomllib
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np
import xarray as xr
from numpy.typing import NDArray
from scipy.optimize import least_squares, minimize_scalar, nnls

import pulse_tube_cooldown
from pulse_tube_cooldown.core.types import Link, Material
from pulse_tube_cooldown.library.gases import helium
from pulse_tube_cooldown.library.links.conduction import Conduction
from pulse_tube_cooldown.library.links.cooler.columns import load_column
from pulse_tube_cooldown.library.links.cooler.drive import calc_mean_pressure, calc_pressure_swing
from pulse_tube_cooldown.library.links.cooler.head import (
    calc_capacity_ratio,
    calc_enthalpy_flow,
    calc_gas_mass,
    calc_heat_flow,
    calc_mass_flow,
    calc_utilization,
)
from pulse_tube_cooldown.library.materials.load import load_material
from pulse_tube_cooldown.library.toml_files import read_named_toml
from pulse_tube_cooldown.measurements.test_stand import read_cooldown, read_load_map

HERE = Path(__file__).parent
SN09 = HERE.parents[2] / "data" / "heads" / "pt1_1"
STAND = HERE.parents[2] / "data/test_stand/s0009"
WARM_FLANGE, STAGE_1, STAGE_2 = "Warm flange", "Stage 1", "Stage 2"


@dataclass(frozen=True, slots=True, kw_only=True)
class Head:
    """A head model's columns and pulse tubes, as the network equations in `cooler/head.py` take them.

    Attributes:
        material: Matrix material of each node in the stage-1 and stage-2 columns, from the warm end.
        mass: Matrix mass in kg of each node in the two columns.
        pore_volume: Pore gas volume in m^3 of each node in the two columns.
        pulse_tube_volume: Gas volumes in m^3 of the stage-1 and stage-2 pulse tubes.
        clamp: Index along the head's chain of the stage-2 column node containing the clamp.
        cold_exchanger_coupling: Whether the matrix cold end couples to the plate through its gas stream.
    """

    material: tuple[tuple[Material, ...], tuple[Material, ...]]
    mass: tuple[NDArray[np.float64], NDArray[np.float64]]
    pore_volume: tuple[NDArray[np.float64], NDArray[np.float64]]
    pulse_tube_volume: tuple[float, float]
    clamp: int
    cold_exchanger_coupling: bool = True

    @property
    def stage_1(self) -> int:
        """Index of stage 1 along the chain."""
        return len(self.mass[0]) + 1

    @property
    def stage_2(self) -> int:
        """Index of stage 2 along the chain."""
        return len(self.mass[0]) + len(self.mass[1]) + 2

    def calc_heat_capacity(self, temperature: NDArray[np.float64]) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Matrix heat capacity in J/K of each column node, from temperatures in K along the chain."""
        node_temperature = (temperature[1 : self.stage_1], temperature[self.stage_1 + 1 : self.stage_2])
        return tuple(
            mass * np.array([material.specific_heat(np.array(t)) for material, t in zip(materials, temps, strict=True)])
            for materials, mass, temps in zip(self.material, self.mass, node_temperature, strict=True)
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class Parameters:
    """One head's fitted values and its valve frequency.

    Attributes:
        packing_constant: Packing constants a_1 and a_2 of the stage-1 and stage-2 columns.
        conductance: Effective conductances C_1 and C_2 of the stage valves in W/Pa^2.
        exchanger_effectiveness: Stage-2 cold heat exchanger's effectiveness eta.
        frequency: Valve frequency in Hz.
    """

    packing_constant: tuple[float, float]
    conductance: tuple[float, float]
    exchanger_effectiveness: float
    frequency: float


@dataclass(frozen=True, slots=True, kw_only=True)
class Drive:
    """The compressors and gas charge on the head's manifold.

    Attributes:
        charge: Helium mass in kg on the manifold.
        warm_volume: Gas volume in m^3 of the manifold's lines and buffers.
        displacement_rate: Volumetric displacement rate in m^3/s of each compressor.
        relief_swing: Relief swing in Pa of each compressor.
        waveform_factor: The drive's waveform factor kappa.
    """

    charge: float
    warm_volume: float
    displacement_rate: tuple[float, ...]
    relief_swing: tuple[float, ...]
    waveform_factor: float = 1.0


@dataclass(frozen=True, slots=True, kw_only=True)
class ColumnSolution:
    """A head at steady state between fixed plates.

    Attributes:
        temperature: Temperatures in K along the head's chain.
        heat_flow: Heat flow in W from the head into each element of the chain, zero at the column nodes.
        supply_pressure: Supply pressure in Pa.
        return_pressure: Return pressure in Pa.
        utilization: Utilization of the stage-1 and stage-2 columns.
        capacity_ratio: Pore-gas to matrix heat-capacity ratio at each node of the two columns.
    """

    temperature: NDArray[np.float64]
    heat_flow: NDArray[np.float64]
    supply_pressure: float
    return_pressure: float
    utilization: tuple[float, float]
    capacity_ratio: tuple[NDArray[np.float64], NDArray[np.float64]]


@dataclass(frozen=True, slots=True, kw_only=True)
class LoadMapFit:
    """A head's parameters fitted to a load map.

    Attributes:
        parameters: The fitted parameters.
        residual: Net heat flow in W into stage 1 and stage 2 at each point, with one row per point.
        columns: The steady-state head at each point.
    """

    parameters: Parameters
    residual: NDArray[np.float64]
    columns: tuple[ColumnSolution, ...]


def read_head(name: str, *, directory: Path = HERE / "head_model") -> Head:
    """Reads a head model file in the column format this experiment was fitted with.

    Args:
        name: The model's `name`, such as "PT1.1".
        directory: Folder holding head model files; defaults to this experiment's snapshot.
    """
    model = read_named_toml(name, directory=directory, kind="heads")
    columns = [
        load_column(model[column], head_name=name, stage=stage)
        for stage, column in enumerate(("stage_1_column", "stage_2_column"), start=1)
    ]
    ends = itertools.accumulate(layer.length_fraction for layer in columns[1])
    clamp_layer = next(index for index, end in enumerate(ends) if model["clamp_fraction"] < end)
    first, second = model["pulse_tube_volumes"]
    return Head(
        material=tuple(tuple(layer.node.composition[0][0] for layer in column) for column in columns),
        mass=tuple(np.array([layer.node.composition[0][1] for layer in column]) for column in columns),
        pore_volume=tuple(np.array([layer.pore_volume for layer in column]) for column in columns),
        pulse_tube_volume=(first, second),
        clamp=len(columns[0]) + 2 + clamp_layer,
    )


def build_walls(geometry: Path = SN09 / "geometry.csv", material: str = "304 SS") -> tuple[Conduction, ...]:
    """Builds the conduction links of a head's four tube walls from its geometry record.

    The stage-1 regenerator and pulse tube join the warm flange to stage 1, the stage-2 regenerator joins stage 1 to
    stage 2, and the stage-2 pulse tube joins the warm flange to stage 2.

    Args:
        geometry: A `geometry.csv` with each tube's outer radius and length in mm and `tube_wall_thickness` in mm.
        material: The tubes' material.
    """
    rows = np.genfromtxt(geometry, delimiter=",", names=True, dtype=None, encoding="utf-8", usecols=(0, 1))
    value = {row["quantity"]: row["value"] for row in rows}
    wall = float(value["tube_wall_thickness"]) * 1e-3
    tubes = {
        "stage1_regenerator": (WARM_FLANGE, STAGE_1),
        "stage1_pulse_tube": (WARM_FLANGE, STAGE_1),
        "stage2_regenerator": (STAGE_1, STAGE_2),
        "stage2_pulse_tube": (WARM_FLANGE, STAGE_2),
    }
    walls = []
    for tube, nodes in tubes.items():
        radius = float(value[f"{tube}_outer_radius"]) * 1e-3
        walls.append(
            Conduction(
                name=f"{tube} wall",
                nodes=nodes,
                material=load_material(material),
                area=np.pi * (radius**2 - (radius - wall) ** 2),
                length=float(value[f"{tube}_length"]) * 1e-3,
            )
        )
    return tuple(walls)


def column_indices(head: Head) -> NDArray[np.int64]:
    """Returns matrix-node indices along the chain."""
    return np.setdiff1d(np.arange(head.stage_2 + 1), [0, head.stage_1, head.stage_2])


def calc_pressures(
    head: Head, parameters: Parameters, drive: Drive, temperature: NDArray[np.float64]
) -> tuple[float, float]:
    """Solves the original drive equations, returning mean pressure and swing in Pa."""
    volumes = (head.pore_volume, head.pulse_tube_volume)
    mean = calc_mean_pressure(
        lambda pressure: calc_gas_mass(temperature, *volumes, pressure),
        drive.charge,
        drive.warm_volume,
        float(temperature[0]),
    )
    swing = calc_pressure_swing(
        lambda supply, back: calc_mass_flow(temperature, *volumes, supply, back, parameters.frequency)[0],
        mean,
        float(temperature[0]),
        drive.displacement_rate,
        drive.relief_swing,
    )
    return mean, swing


def evaluate_head(
    head: Head, parameters: Parameters, drive: Drive, temperature: NDArray[np.float64], mean: float, swing: float
) -> ColumnSolution:
    """Evaluates the library's head equations at explicit pressures."""
    mass_flow = calc_mass_flow(
        temperature, head.pore_volume, head.pulse_tube_volume, mean + swing / 2, mean - swing / 2, parameters.frequency
    )
    heat_capacity = head.calc_heat_capacity(temperature)
    capacity_ratio = calc_capacity_ratio(temperature, heat_capacity, head.pore_volume, mean)
    enthalpy = tuple(calc_enthalpy_flow(c, swing, drive.waveform_factor) for c in parameters.conductance)
    return ColumnSolution(
        temperature=temperature.copy(),
        heat_flow=(calc_heat_flow if head.cold_exchanger_coupling else calc_uncoupled_heat_flow)(
            temperature,
            heat_capacity,
            capacity_ratio,
            mass_flow,
            enthalpy,
            parameters.packing_constant,
            parameters.exchanger_effectiveness,
            mean,
            parameters.frequency,
        ),
        supply_pressure=mean + swing / 2,
        return_pressure=mean - swing / 2,
        utilization=calc_utilization(temperature, heat_capacity, mass_flow, mean, parameters.frequency),
        capacity_ratio=capacity_ratio,
    )


def calc_uncoupled_heat_flow(
    temperature,
    heat_capacity,
    capacity_ratio,
    mass_flow,
    enthalpy_flow,
    packing_constant,
    exchanger_effectiveness,
    mean_pressure,
    frequency,
):
    """Finite-capacity control: the column ends on the plate, with finite-capacity terms."""
    stage_1 = len(heat_capacity[0]) + 1
    capacity = np.r_[0, heat_capacity[0], 0, heat_capacity[1], 0]
    midpoint = (temperature[:-1] + temperature[1:]) / 2
    cp = helium.isobaric_specific_heat(midpoint, mean_pressure)
    packing = np.repeat(packing_constant, [stage_1, len(mass_flow) - stage_1])
    flow = 2 * packing * (mass_flow * cp) ** 2 / (frequency * (capacity[:-1] + capacity[1:]))
    flow *= temperature[:-1] - temperature[1:]
    boundary = midpoint[stage_1:].copy()
    boundary[[0, -1]] = temperature[[stage_1, -1]]
    density = helium.density(boundary, mean_pressure)
    eps = capacity_ratio[1]
    attenuation = np.log(density[1:] / density[:-1]) / (1 + eps)
    acoustic = enthalpy_flow[1] * np.exp(np.r_[np.cumsum(attenuation[::-1])[::-1], 0])
    boundary_eps = np.r_[eps[0], (eps[:-1] + eps[1:]) / 2, eps[-1]]
    expansion = helium.isobaric_expansion_coefficient(boundary, mean_pressure) * boundary
    flow[stage_1:] += (1 - expansion / (1 + boundary_eps)) * acoustic
    result = -np.diff(flow, prepend=0.0, append=0.0)
    result[stage_1] -= enthalpy_flow[0]
    result[-1] -= enthalpy_flow[1]
    result[0] += sum(enthalpy_flow)
    return result


def solve_columns(
    head: Head,
    parameters: Parameters,
    drive: Drive,
    *,
    warm_temperature: float,
    stage_temperature: tuple[float, float],
    guess: NDArray[np.float64] | None = None,
    tolerance: float = 1e-7,
    temperature_range: tuple[float, float] = (2.0, 300.0),
    max_nfev: int = 200,
) -> ColumnSolution:
    """Solves matrix balances, gas inventory, and the capped compressor drive together.

    Matrix temperatures lie in the supported material range. The internal equations
    are normalized to unit packing constants and retain C_2 / a_2 and eta / a_2. Mean
    pressure and swing are unknowns in units of 1 MPa; thermal residuals use a 100 W scale.
    """
    nodes = column_indices(head)
    normalized = Parameters(
        packing_constant=(1.0, 1.0),
        conductance=(0.0, parameters.conductance[1] / parameters.packing_constant[1]),
        exchanger_effectiveness=parameters.exchanger_effectiveness / parameters.packing_constant[1],
        frequency=parameters.frequency,
    )
    template = np.empty(head.stage_2 + 1)
    template[[0, head.stage_1, head.stage_2]] = [warm_temperature, *stage_temperature]
    if guess is None:
        guess = template.copy()
        for warm, cold in ((0, head.stage_1), (head.stage_1, head.stage_2)):
            fraction = (np.arange(warm + 1, cold) - warm - 0.5) / (cold - warm - 1)
            guess[warm + 1 : cold] = template[warm] * (template[cold] / template[warm]) ** fraction
        count = head.stage_2 - head.stage_1 - 1
        fraction = (np.arange(count) + 0.5) / count
        guess[head.stage_1 + 1 : head.stage_2] = stage_temperature[0] * (1 - fraction) + stage_temperature[1] * fraction
        guess[head.stage_2 - 1] = max(9.0, guess[head.stage_2 - 1])
    guess = guess.copy()
    guess[[0, head.stage_1, head.stage_2]] = template[[0, head.stage_1, head.stage_2]]
    mean, swing = calc_pressures(head, parameters, drive, guess)
    cap = min(drive.relief_swing)
    lower = np.r_[np.full(len(nodes), np.log(temperature_range[0])), (cap / 2 + 1e4) / 1e6, 1e-4]
    upper = np.r_[np.full(len(nodes), np.log(temperature_range[1])), 2.5, cap / 1e6]
    initial = np.r_[np.log(guess[nodes]), mean / 1e6, swing / 1e6]

    def residual(state: NDArray[np.float64]) -> NDArray[np.float64]:
        temperature = template.copy()
        temperature[nodes] = np.exp(state[:-2])
        pressure, delta = state[-2:] * 1e6
        column = evaluate_head(head, normalized, drive, temperature, pressure, delta)
        inventory = calc_gas_mass(temperature, head.pore_volume, head.pulse_tube_volume, pressure)
        inventory += drive.warm_volume * float(helium.density(warm_temperature, pressure))
        head_flow = calc_mass_flow(
            temperature,
            head.pore_volume,
            head.pulse_tube_volume,
            pressure + delta / 2,
            pressure - delta / 2,
            parameters.frequency,
        )[0]
        compressor_flow = float(helium.density(warm_temperature, pressure - delta / 2)) * sum(drive.displacement_rate)
        compressor = min((cap - delta) / cap, (compressor_flow - head_flow) / 0.005)
        return np.r_[column.heat_flow[nodes] / 100, (inventory - drive.charge) / drive.charge, compressor]

    solution = least_squares(
        residual,
        np.clip(initial, lower + 1e-12, upper - 1e-12),
        bounds=(lower, upper),
        x_scale="jac",
        diff_step=1e-5,
        max_nfev=max_nfev,
        xtol=1e-12,
        ftol=1e-12,
        gtol=1e-12,
    )
    temperature = template.copy()
    temperature[nodes] = np.exp(solution.x[:-2])
    mean, swing = calc_pressures(head, parameters, drive, temperature)
    result = evaluate_head(head, parameters, drive, temperature, mean, swing)
    balance = max(abs(result.heat_flow[nodes]))
    constraint = max(abs(residual(solution.x)[-2:]))
    if balance > tolerance:
        raise RuntimeError(
            f"Column balance at {stage_temperature}: {balance:.6g} W; drive residual {constraint:.3g}; "
            f"ratio {normalized.conductance[1]:.6g}; temperatures {temperature[nodes]}."
        )
    return result


def calc_growth_rate(
    head: Head, parameters: Parameters, drive: Drive, column: ColumnSolution, relative_step: float = 1e-5
) -> float:
    """Returns the largest real thermal eigenvalue in 1/s, including drive response."""
    nodes = column_indices(head)
    capacities = np.concatenate(head.calc_heat_capacity(column.temperature))
    jacobian = np.empty((len(nodes), len(nodes)))
    for index, node in enumerate(nodes):
        step = relative_step * column.temperature[node]
        flows = []
        for sign in (1, -1):
            temperature = column.temperature.copy()
            temperature[node] += sign * step
            mean, swing = calc_pressures(head, parameters, drive, temperature)
            flows.append(evaluate_head(head, parameters, drive, temperature, mean, swing).heat_flow[nodes])
        jacobian[:, index] = (flows[0] - flows[1]) / (2 * step * capacities)
    return float(np.linalg.eigvals(jacobian).real.max())


def read_experiment() -> tuple[dict, xr.Dataset, float, Drive]:
    """Reads the May measurements and converts the recorded fill to helium mass."""
    inputs = tomllib.loads((HERE / "inputs.toml").read_text(encoding="utf-8"))
    head = read_head(inputs["head"])
    load_map = read_load_map(STAND / inputs["load_map"])
    warm = float(read_cooldown(STAND / inputs["cooldown"]).TPR[0])
    charge = calc_gas_mass(
        np.full(head.stage_2 + 1, warm), head.pore_volume, head.pulse_tube_volume, inputs["fill_pressure"]
    )
    charge += inputs["warm_volume"] * float(helium.density(warm, inputs["fill_pressure"]))
    drive = Drive(
        charge=charge,
        warm_volume=inputs["warm_volume"],
        displacement_rate=tuple(inputs["displacement_rate"]),
        relief_swing=tuple(inputs["relief_swing"]),
        waveform_factor=inputs["waveform_factor"],
    )
    return inputs, load_map, warm, drive


def calc_link_heat_flow(links: Sequence[Link], temperature: Mapping[str, float]) -> dict[str, float]:
    """Sums the heat flow in W from links into each node they touch.

    Args:
        links: Links over named nodes.
        temperature: Temperature in K of each node the links touch.
    """
    heat_flow = dict.fromkeys(temperature, 0.0)
    for link in links:
        flows = link.calc_heat_flow(np.array([temperature[node] for node in link.nodes]))
        for node, flow in zip(link.nodes, flows, strict=True):
            heat_flow[node] += float(flow)
    return heat_flow


def calc_stage_heat_flow(
    head: Head,
    parameters: Parameters,
    drive: Drive,
    links: Sequence[Link],
    *,
    warm_temperature: float,
    stage_temperature: tuple[float, float],
    guess: NDArray[np.float64] | None = None,
) -> tuple[NDArray[np.float64], ColumnSolution]:
    """Calculates the heat flow into stage 1 and stage 2 from the head and the links at steady state.

    Args:
        head: The head model.
        parameters: The head's parameters.
        drive: The head's manifold.
        links: Links over the warm flange and stages, named by `WARM_FLANGE`, `STAGE_1`, and `STAGE_2`.
        warm_temperature: Warm-flange temperature in K.
        stage_temperature: Stage-1 and stage-2 temperatures in K.
        guess: Starting temperatures in K along the chain for the column solve.

    Returns:
        The heat flows in W into stage 1 and stage 2, and the head's steady state.
    """
    column = solve_columns(
        head, parameters, drive, warm_temperature=warm_temperature, stage_temperature=stage_temperature, guess=guess
    )
    plates = {WARM_FLANGE: warm_temperature, STAGE_1: stage_temperature[0], STAGE_2: stage_temperature[1]}
    link_flow = calc_link_heat_flow(links, plates)
    head_flow = column.heat_flow[[head.stage_1, head.stage_2]]
    return head_flow + np.array([link_flow[STAGE_1], link_flow[STAGE_2]]), column


def fit_at_ratio(
    head: Head,
    drive: Drive,
    links: Sequence[Link],
    load_map: xr.Dataset,
    *,
    warm_temperature: float,
    ratio: float,
    effectiveness_ratio: float,
    frequency: float,
    guesses: Sequence[NDArray[np.float64] | None] | None = None,
) -> LoadMapFit:
    """Fits three nonnegative linear coefficients after solving profiles at C_2 / a_2 and eta / a_2."""
    reference_conductance = 1e-10
    stages = np.column_stack([load_map.T1.values, load_map.T2.values])
    heater = np.column_stack([load_map.Q1.values, load_map.Q2.values])
    normalized = Parameters(
        packing_constant=(1.0, 1.0),
        conductance=(0.0, ratio),
        exchanger_effectiveness=effectiveness_ratio,
        frequency=frequency,
    )
    if guesses is None:
        guesses = [None] * len(stages)
    profiles, matrices, target = [], [], []
    for point, stage in enumerate(stages):
        column = solve_columns(
            head,
            normalized,
            drive,
            warm_temperature=warm_temperature,
            stage_temperature=tuple(stage),
            guess=guesses[point],
        )
        profiles.append(column)
        mean = (column.supply_pressure + column.return_pressure) / 2
        swing = column.supply_pressure - column.return_pressure
        basis = [
            Parameters(
                packing_constant=(1.0, 0.0), conductance=(0.0, 0.0), exchanger_effectiveness=0.0, frequency=frequency
            ),
            Parameters(
                packing_constant=(0.0, 1.0),
                conductance=(0.0, ratio),
                exchanger_effectiveness=effectiveness_ratio,
                frequency=frequency,
            ),
            Parameters(
                packing_constant=(0.0, 0.0),
                conductance=(reference_conductance, 0.0),
                exchanger_effectiveness=0.0,
                frequency=frequency,
            ),
        ]
        matrices.append(
            np.column_stack(
                [
                    evaluate_head(head, p, drive, column.temperature, mean, swing).heat_flow[
                        [head.stage_1, head.stage_2]
                    ]
                    for p in basis
                ]
            )
        )
        wall = calc_link_heat_flow(links, {WARM_FLANGE: warm_temperature, STAGE_1: stage[0], STAGE_2: stage[1]})
        target.append(-heater[point] - [wall[STAGE_1], wall[STAGE_2]])
    matrix = np.concatenate(matrices)
    target = np.concatenate(target)
    coefficients, _ = nnls(matrix, target)
    a1, a2, c1 = coefficients
    parameters = Parameters(
        packing_constant=(float(a1), float(a2)),
        conductance=(float(c1 * reference_conductance), float(ratio * a2)),
        exchanger_effectiveness=float(effectiveness_ratio * a2),
        frequency=frequency,
    )
    columns = tuple(
        evaluate_head(
            head,
            parameters,
            drive,
            column.temperature,
            (column.supply_pressure + column.return_pressure) / 2,
            column.supply_pressure - column.return_pressure,
        )
        for column in profiles
    )
    residual = (matrix @ coefficients - target).reshape(-1, 2)
    return LoadMapFit(parameters=parameters, residual=residual, columns=columns)


def fit_load_map(
    head: Head,
    drive: Drive,
    links: Sequence[Link],
    load_map: xr.Dataset,
    *,
    warm_temperature: float,
    initial: Parameters,
    ratio_bounds: tuple[float, float] | None = None,
    guesses: Sequence[NDArray[np.float64] | None] | None = None,
    verbose: bool = False,
) -> LoadMapFit:
    """Searches the two nonlinear ratios and solves three coefficients linearly.

    C_2 / a_2 and eta / a_2 determine the equilibrium profiles. Each evaluation starts
    from the nearest solved profiles in log-ratio coordinates. The uncoupled control
    searches C_2 / a_2 alone. Verbose output reports the ratios, eta, and stage RMS.
    """
    reference = initial.conductance[1] / initial.packing_constant[1]
    effectiveness_ratio = initial.exchanger_effectiveness / initial.packing_constant[1]
    seed = fit_at_ratio(
        head,
        drive,
        links,
        load_map,
        warm_temperature=warm_temperature,
        ratio=reference,
        effectiveness_ratio=effectiveness_ratio,
        frequency=initial.frequency,
        guesses=guesses,
    )
    bounds = ratio_bounds if ratio_bounds is not None else (reference / 10, reference * 10)
    evaluated = {(0.0, 0.0): seed}

    def objective(log_ratios):
        ratio = reference * np.exp(log_ratios[0])
        nearest = min(evaluated, key=lambda key: np.linalg.norm(np.array(key) - log_ratios))
        starting_profiles = [c.temperature for c in evaluated[nearest].columns]
        result = fit_at_ratio(
            head,
            drive,
            links,
            load_map,
            warm_temperature=warm_temperature,
            ratio=ratio,
            effectiveness_ratio=effectiveness_ratio * np.exp(log_ratios[1]) if head.cold_exchanger_coupling else 0.0,
            frequency=initial.frequency,
            guesses=starting_profiles,
        )
        evaluated[tuple(log_ratios)] = result
        if verbose:
            print(
                "FIT",
                ratio,
                result.parameters.exchanger_effectiveness,
                np.sqrt(np.mean(result.residual**2, axis=0)),
                flush=True,
            )
        return result.residual.ravel()

    if head.cold_exchanger_coupling:

        def jacobian(log_ratios):
            step = 1e-3
            return np.column_stack(
                [
                    (objective(log_ratios + np.eye(2)[j] * step) - objective(log_ratios - np.eye(2)[j] * step))
                    / (2 * step)
                    for j in range(2)
                ]
            )

        solution = least_squares(
            objective,
            np.zeros(2),
            bounds=([np.log(bounds[0] / reference), -8.0], [np.log(bounds[1] / reference), 8.0]),
            jac=jacobian,
            xtol=1e-7,
            ftol=1e-7,
            gtol=1e-7,
            max_nfev=60,
        )
    else:
        solution = minimize_scalar(
            lambda x: np.sum(objective([x, 0.0]) ** 2),
            bounds=tuple(np.log(np.array(bounds) / reference)),
            method="bounded",
            options={"xatol": 1e-7},
        )
        solution.x = np.array([solution.x, 0.0])
    if not solution.success:
        raise RuntimeError(f"Ratio search failed: {solution.message}")
    return evaluated[tuple(solution.x)]


def _write_parameters(path: Path, parameters: Parameters, sources: Mapping[str, str]) -> None:
    """Writes a head's parameter file."""
    a_1, a_2 = parameters.packing_constant
    c_1, c_2 = parameters.conductance
    path.write_text(
        f"# Written by simulations/validation/sn09_load_map/fit.py.\n"
        f"packing_constants = [{a_1:.12g}, {a_2:.12g}]\n"
        f"conductances = [{c_1:.12e}, {c_2:.12e}]  # W/Pa^2\n"
        f"exchanger_effectiveness = {parameters.exchanger_effectiveness:.12g}\n"
        f"frequency = {parameters.frequency}  # Hz\n"
        f"\n[sources]\n" + "".join(f'{key} = "{value}"\n' for key, value in sources.items()),
        encoding="utf-8",
    )


def write_results(
    head: Head,
    drive: Drive,
    load_map: xr.Dataset,
    result: LoadMapFit,
    inputs: Mapping,
    warm: float,
    elapsed: float,
    *,
    folder: Path = HERE,
) -> None:
    """Writes the fitted branch, physical checks, and power residuals beside the experiment."""
    rows = []
    for point, column in enumerate(result.columns):
        mean, swing = calc_pressures(head, result.parameters, drive, column.temperature)
        checked = evaluate_head(head, result.parameters, drive, column.temperature, mean, swing)
        growth = calc_growth_rate(head, result.parameters, drive, column)
        rows.append(
            dict(
                point=point,
                heater_power_W=[float(load_map.Q1[point]), float(load_map.Q2[point])],
                stage_temperature_K=[float(load_map.T1[point]), float(load_map.T2[point])],
                residual_W=result.residual[point].tolist(),
                TPR_K=float(load_map.TPR[point]),
                clamp_temperature_K=float(column.temperature[head.clamp]),
                temperature_K=column.temperature.tolist(),
                supply_pressure_Pa=column.supply_pressure,
                return_pressure_Pa=column.return_pressure,
                utilization=list(column.utilization),
                capacity_ratio=[ratio.tolist() for ratio in column.capacity_ratio],
                max_column_balance_W=float(max(abs(checked.heat_flow[column_indices(head)]))),
                pressure_check_relative=float(
                    max(
                        abs(
                            np.array([mean + swing / 2, mean - swing / 2])
                            / [column.supply_pressure, column.return_pressure]
                            - 1
                        )
                    )
                ),
                growth_rate_per_s=growth,
            )
        )
    if max(row["growth_rate_per_s"] for row in rows) >= 0:
        raise RuntimeError("The fitted branch contains a thermally unstable equilibrium.")
    if max(row["max_column_balance_W"] for row in rows) > 1e-6:
        raise RuntimeError("Independent drive recomputation left a column out of balance.")
    source_paths = [
        Path(__file__),
        HERE / "inputs.toml",
        HERE / inputs["search" if head.cold_exchanger_coupling else "uncoupled_search"]["profile_seed"],
        STAND / inputs["load_map"],
        STAND / inputs["cooldown"],
        SN09 / "geometry.csv",
    ]
    package = Path(pulse_tube_cooldown.__file__).parent
    source_paths += [
        package / "library/links/cooler/head.py",
        package / "library/links/cooler/drive.py",
        package / "library/links/cooler/columns.py",
        package / "library/gases/helium.py",
        package / "library/property_fits.py",
        HERE / "head_model/pt1_1.toml",
    ]
    source_paths += [
        package / "library/materials/data" / name
        for name in ["ss_304.toml", "pb.toml", "er3ni.toml", "hocu2.toml", "gd2o2s.toml"]
    ]
    source_paths += [package / "library/links/conduction.py", package / "library/materials/load.py"]
    root = HERE.parents[2]
    provenance = {str(path.relative_to(root)): hashlib.sha256(path.read_bytes()).hexdigest() for path in source_paths}
    parameters = result.parameters
    summary = dict(
        status="converged_stable_branch_fit",
        global_minimum_established=False,
        cold_exchanger_coupling=head.cold_exchanger_coupling,
        warm_temperature_K=warm,
        helium_charge_kg=drive.charge,
        elapsed_s=elapsed,
        ratio_C2_over_a2=parameters.conductance[1] / parameters.packing_constant[1],
        packing_constants=list(parameters.packing_constant),
        conductances_W_per_Pa2=list(parameters.conductance),
        exchanger_effectiveness=parameters.exchanger_effectiveness,
        frequency_Hz=parameters.frequency,
        sum_squared_residual_W2=float(np.sum(result.residual**2)),
        rms_residual_W=np.sqrt(np.mean(result.residual**2, axis=0)).tolist(),
        max_abs_residual_W=np.max(abs(result.residual), axis=0).tolist(),
        points=rows,
        sha256=provenance,
    )
    folder.mkdir(exist_ok=True)
    (folder / "results.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    sources = dict(inputs["sources"])
    sources["packing_constants"] = "May map; least squares in W on the cooling branch, see results.json"
    sources["conductances"] = sources["packing_constants"]
    sources["exchanger_effectiveness"] = (
        "May map; eta / a_2 and C_2 / a_2 jointly fitted"
        if head.cold_exchanger_coupling
        else "Uncoupled finite-capacity control"
    )
    _write_parameters(folder / "candidate_parameters.toml", parameters, sources)
    lines = [
        "SN09 May load-map fit: finite capacity"
        + (" and cold-exchanger coupling" if head.cold_exchanger_coupling else " only"),
        "",
        "A converged fit on a stable cooling branch; a global minimum has not been established.",
        "Residual = head heat into a stage + wall heat into it + its measured heater power.",
        f"Warm flange: {warm:.6f} K; helium charge: {drive.charge:.9g} kg.",
        f"Fit time: {elapsed:.3f} s.",
        f"Packing constants: {parameters.packing_constant}",
        f"Exchanger effectiveness eta: {parameters.exchanger_effectiveness:.9g}.",
        f"Valve conductances: {parameters.conductance} W/Pa^2; frequency: {parameters.frequency} Hz.",
        f"Sum of squared residuals: {summary['sum_squared_residual_W2']:.9g} W^2.",
        f"RMS stage residuals: {summary['rms_residual_W']} W.",
        "",
        " Q1/W Q2/W   T1/K   T2/K   residual1/W residual2/W  clamp/K TPR/K growth/(1/s)",
    ]
    for row in rows:
        lines.append(
            f"{row['heater_power_W'][0]:5.1f} {row['heater_power_W'][1]:4.1f} "
            f"{row['stage_temperature_K'][0]:7.3f} {row['stage_temperature_K'][1]:6.3f} "
            f"{row['residual_W'][0]:11.5f} {row['residual_W'][1]:11.5f} "
            f"{row['clamp_temperature_K']:8.3f} {row['TPR_K']:6.3f} {row['growth_rate_per_s']:12.5g}"
        )
    lines += [
        "",
        "All matrix nodes balance with pressures independently recomputed by the library drive solver.",
        "All largest real thermal eigenvalues are negative. The warm flange and stages are held fixed.",
        "These checks establish local stability, rather than reachability from room temperature.",
        "",
        "Baseline and contribution comparison: comparison.txt.",
        "The candidate remains in this experiment directory; production parameters have not been written.",
        "Older branch scans here predate the finite-capacity terms and exchanger coupling; baseline/ keeps that fit.",
        "Zero-a_2 fits are degenerate: their stage-2 matrix has no thermal restoring rate.",
        "",
        "Reproduce from the repository root:",
        "  .venv/Scripts/python simulations/validation/sn09_load_map/fit.py",
        "  .venv/Scripts/python simulations/validation/sn09_load_map/verify.py",
        "  .venv/Scripts/python -m pytest simulations/validation/sn09_load_map/test_fit.py -q -p no:cacheprovider",
        "Run explore.py with seed, scan, or mesh for the current equations.",
    ]
    report = "\n".join(lines) + "\n"
    (folder / "fit_report.txt").write_text(report, encoding="utf-8")
    print(report, flush=True)
    plot_results(load_map, result, folder=folder)


def plot_results(load_map: xr.Dataset, result: LoadMapFit, *, folder: Path = HERE) -> None:
    """Plots heater power, predicted available cooling, and stage residuals."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(10, 7), layout="constrained")
    for stage in range(2):
        temperature = load_map[f"T{stage + 1}"].values
        heater = load_map[f"Q{stage + 1}"].values
        prediction = heater - result.residual[:, stage]
        axes[0, stage].scatter(temperature, heater, marker="o", label="Measured heater")
        axes[0, stage].scatter(temperature, prediction, marker="x", label="Model cooling after walls")
        axes[0, stage].set(title=f"Stage {stage + 1}", xlabel="Measured temperature (K)", ylabel="Power (W)")
        axes[0, stage].legend(fontsize=8)
        axes[1, stage].axhline(0, color="0.5", linewidth=0.8)
        axes[1, stage].scatter(temperature, result.residual[:, stage])
        axes[1, stage].set(xlabel="Measured temperature (K)", ylabel="Net stage imbalance (W)")
    fig.suptitle("SN09 May map: fit on a stable cooling branch")
    fig.savefig(folder / "load_map_fit.png", dpi=180)
    plt.close(fig)


def main() -> None:
    """Fits the selected branch and writes all results inside this experiment directory."""
    inputs, load_map, warm, drive = read_experiment()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--uncoupled", action="store_true", help="Fit the uncoupled control, with the column ending on the plate."
    )
    args = parser.parse_args()
    head = replace(read_head(inputs["head"]), cold_exchanger_coupling=not args.uncoupled)
    search = inputs["uncoupled_search"] if args.uncoupled else inputs["search"]
    initial_data = search.get("initial", inputs["initial"])
    initial = Parameters(
        packing_constant=tuple(initial_data["packing_constant"]),
        conductance=tuple(initial_data["conductance"]),
        exchanger_effectiveness=initial_data["exchanger_effectiveness"],
        frequency=inputs["frequency"],
    )
    seeds = json.loads((HERE / search["profile_seed"]).read_text(encoding="utf-8"))
    start = time.perf_counter()
    result = fit_load_map(
        head,
        drive,
        build_walls(),
        load_map,
        warm_temperature=warm,
        initial=initial,
        ratio_bounds=tuple(search["ratio_bounds"]),
        guesses=[np.array(t) for t in seeds["temperatures"]],
        verbose=True,
    )
    write_results(
        head,
        drive,
        load_map,
        result,
        inputs,
        warm,
        time.perf_counter() - start,
        folder=HERE / "uncoupled" if args.uncoupled else HERE,
    )


if __name__ == "__main__":
    main()
