"""Fitting a head's coefficients to the temperature-crossing times of a test-stand cooldown."""

import logging
from dataclasses import asdict, fields, replace
from datetime import datetime
from typing import Any

import numpy as np
import xarray as xr
from numpy.typing import NDArray
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares

from pulse_tube_cooldown.core.errors import ConfigurationError, SimulationError
from pulse_tube_cooldown.core.types import Graph
from pulse_tube_cooldown.engine.equations import calc_temperature_rate_of_change
from pulse_tube_cooldown.library.links.pulse_tube.empirical.head import PulseTubeHead, build_pulse_tube_head
from pulse_tube_cooldown.library.links.pulse_tube.empirical.parameters import (
    FitParameters,
    HeadCalibration,
    HeadCoefficients,
    load_fit_parameters,
)

__all__ = ["first_crossing_time", "fit_pulse_tube_head", "measured_crossing_times"]


def fit_pulse_tube_head(
    initial: HeadCalibration,
    cooldown: xr.Dataset,
    /,
    *,
    test_stand: Graph,
    start: datetime | np.datetime64 | str,
    end: datetime | np.datetime64 | str,
    parameters: FitParameters | None = None,
) -> HeadCalibration:
    """Fits the coefficients to the relative times of the first downward temperature crossings of a cooldown.

    The simulated cooldown starts with every node at the calibration's warm-flange temperature.

    Args:
        initial: The calibration holding the initial coefficients.
        cooldown: T1 and T2 in K over datetime `time`, as returned by `read_cooldown`.
        test_stand: The calibration's test stand, built from `initial.test_stand`.
        start: Compressor start, the zero of elapsed time.
        end: The last reading compared.
        parameters: Fit settings. Defaults to the packaged settings.

    Returns:
        The calibration with the fitted coefficients.

    Raises:
        ConfigurationError: The cooldown misses a crossing.
        SimulationError: Integration or optimization fails, or the fitted head misses a crossing.
    """
    parameters = parameters or load_fit_parameters()
    warm = initial.warm_flange_temperature
    start, end = np.datetime64(start), np.datetime64(end)
    measured = measured_crossing_times(cooldown, start, end, parameters=parameters, stage_2_max_temperature=warm)
    targets = np.array([*parameters.stage_1_crossings, *parameters.stage_2_crossings])
    stage = np.repeat([0, 1], [len(parameters.stage_1_crossings), len(parameters.stage_2_crossings)])

    record = (cooldown.time.sel(time=slice(start, end)).values[-1] - start) / np.timedelta64(1, "s")
    duration = parameters.duration or max(
        parameters.record_duration_multiplier * record, parameters.crossing_duration_multiplier * measured.max()
    )
    time = np.linspace(0.0, duration, parameters.sample_points)
    head = build_pulse_tube_head(
        initial, test_stand, name="Head", stages=initial.test_stand_stages, warm_flange_temperature=warm
    )

    def simulate(
        coefficients: HeadCoefficients, solver: dict[str, Any]
    ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        trace = simulate_stage_temperatures(replace(head, coefficients=coefficients), test_stand, time, solver)
        crossings = np.array(
            [first_crossing_time(time, trace[:, index], target) for index, target in zip(stage, targets, strict=True)]
        )
        return crossings, trace

    def residual(coordinates: NDArray[np.float64]) -> NDArray[np.float64]:
        crossings, trace = simulate(from_coordinates(coordinates), parameters.solver)
        missing = np.isnan(crossings)
        crossings[missing] = duration * trace[-1, stage[missing]] / targets[missing]
        return (crossings - measured) / measured

    fit = least_squares(
        residual,
        to_coordinates(initial.coefficients),
        bounds=(to_coordinates(parameters.lower_bounds), to_coordinates(parameters.upper_bounds)),
        **parameters.optimizer,
    )
    if not fit.success:
        raise SimulationError(f"Head calibration optimization failed: {fit.message}")

    fitted = from_coordinates(fit.x)
    final_solver = parameters.solver | {
        key: min(parameters.solver[key], value) for key, value in parameters.final_solver.items()
    }
    crossings, _ = simulate(fitted, final_solver)
    if np.isnan(crossings).any():
        raise SimulationError("The fitted head misses a temperature crossing.")
    error = (crossings - measured) / measured
    logging.getLogger(__name__).info(
        "Head calibration objective %.6g; worst relative crossing-time error %.3g", error @ error, np.abs(error).max()
    )
    provenance = {**initial.provenance, "compressor_start": str(start), "comparison_end": str(end)}
    return replace(initial, coefficients=fitted, provenance=provenance)


def measured_crossing_times(
    cooldown: xr.Dataset,
    start: datetime | np.datetime64 | str,
    end: datetime | np.datetime64 | str,
    *,
    parameters: FitParameters,
    stage_2_max_temperature: float,
) -> NDArray[np.float64]:
    """Returns the first downward crossing times in s after `start`, stage-1 crossings first.

    Readings after `end`, nonfinite readings, stage-1 readings at or above `parameters.stage_1_max_temperature`, and
    stage-2 readings above `stage_2_max_temperature` are excluded.

    Raises:
        ConfigurationError: The cooldown misses a crossing.
    """
    start = np.datetime64(start)
    window = cooldown.sel(time=slice(start, np.datetime64(end)))
    time = (window.time.values - start) / np.timedelta64(1, "s")
    stage_1, stage_2 = window.T1.values, window.T2.values
    retained_1 = np.isfinite(stage_1) & (stage_1 < parameters.stage_1_max_temperature)
    retained_2 = np.isfinite(stage_2) & (stage_2 <= stage_2_max_temperature)
    traces = {1: (time[retained_1], stage_1[retained_1]), 2: (time[retained_2], stage_2[retained_2])}
    crossings = [(1, target) for target in parameters.stage_1_crossings] + [
        (2, target) for target in parameters.stage_2_crossings
    ]
    measured = np.array([first_crossing_time(*traces[stage], target) for stage, target in crossings])
    missing = [crossing for crossing, value in zip(crossings, measured, strict=True) if np.isnan(value)]
    if missing:
        raise ConfigurationError(f"The cooldown misses crossings at (stage, temperature in K): {missing}.")
    return measured


def first_crossing_time(time: NDArray[np.float64], temperature: NDArray[np.float64], target: float) -> float:
    """Returns the time at which `temperature` first falls to `target`, interpolated linearly, or NaN."""
    downward = np.flatnonzero((temperature[:-1] > target) & (temperature[1:] <= target))
    if not downward.size:
        return np.nan
    index = downward[0]
    return float(np.interp(target, temperature[index : index + 2][::-1], time[index : index + 2][::-1]))


def simulate_stage_temperatures(
    head: PulseTubeHead, test_stand: Graph, time: NDArray[np.float64], solver: dict[str, Any]
) -> NDArray[np.float64]:
    """Integrates a test-stand cooldown from a uniform warm-flange temperature, in logarithmic temperatures.

    The test stand's nodes other than the head's stages are held at the warm-flange temperature.

    Returns:
        Stage temperatures in K at `time`, shape (times, 2).

    Raises:
        SimulationError: The integration fails.
    """
    graph = Graph(nodes=(*test_stand.nodes, *head.internal_nodes), links=(*test_stand.links, head))
    stages = [test_stand.get_node_index(stage) for stage in head.stages]
    boundary = tuple(index for index in range(len(test_stand.nodes)) if index not in stages)
    warm = head.warm_flange_temperature

    def rate(elapsed: float, state: NDArray[np.float64]) -> NDArray[np.float64]:
        temperature = warm * np.exp(state)
        return calc_temperature_rate_of_change(graph, boundary, elapsed, temperature) / temperature

    solution = solve_ivp(rate, (0.0, time[-1]), np.zeros(len(graph.nodes)), dense_output=True, **solver)
    if not solution.success:
        raise SimulationError(f"Head calibration integration failed: {solution.message}")
    return warm * np.exp(solution.sol(time)[stages].T)


def to_coordinates(coefficients: HeadCoefficients) -> NDArray[np.float64]:
    """Returns the optimizer's coordinates: the logarithms of the positive coefficients, then the energy coupling."""
    values = asdict(coefficients)
    coupling = values.pop("energy_coupling_coefficient")
    return np.array([*np.log(list(values.values())), coupling])


def from_coordinates(coordinates: NDArray[np.float64]) -> HeadCoefficients:
    """Returns the coefficients at the optimizer's coordinates."""
    names = [field.name for field in fields(HeadCoefficients) if field.name != "energy_coupling_coefficient"]
    positive = dict(zip(names, map(float, np.exp(coordinates[:-1])), strict=True))
    return HeadCoefficients(**positive, energy_coupling_coefficient=float(coordinates[-1]))
