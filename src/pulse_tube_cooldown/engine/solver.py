"""Integration of a model's phases through time."""

from collections.abc import Callable
from functools import partial

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

from pulse_tube_cooldown.core.errors import ConfigurationError
from pulse_tube_cooldown.core.types import Model, Phase, PhaseSolution
from pulse_tube_cooldown.engine.equations import calc_temperature_rate_of_change

__all__ = ["integrate_phase", "run_procedure"]


def integrate_phase(model: Model, phase: Phase, time: float, temperature: NDArray[np.float64]) -> PhaseSolution:
    """Integrates node temperatures through one phase.

    Args:
        model: The model containing `phase`. Its `solver_options` require `max_time`, the time in s at which the run
            ends, and `save_every`, the spacing in s of samples from the phase's start. `solver` names the
            `scipy.integrate.solve_ivp` method and defaults to "Radau". Other options pass to `solve_ivp`.
        phase: The phase to integrate, using its graph in `model`.
        time: Time in s at which the phase starts.
        temperature: Temperature in K of each node at `time`.

    Returns:
        Samples every `save_every` from `time`, followed by the state where the phase ended.

    Raises:
        ConfigurationError: `max_time` or `save_every` is missing from the solver options.
    """
    options = dict(model.solver_options)
    missing = [name for name in ("max_time", "save_every") if name not in options]
    if missing:
        raise ConfigurationError(f"Solver options need {missing}; options are {sorted(options)}.")

    max_time = options.pop("max_time")
    save_every = options.pop("save_every")
    method = options.pop("solver", "Radau")

    graph = model.graphs[model.procedure.phases.index(phase)]
    conditions = {"until": phase.until, "stop_condition": model.procedure.stop_condition}
    ends = [end for end, condition in conditions.items() if condition is not None]

    solution = solve_ivp(
        partial(calc_temperature_rate_of_change, graph, model.fixed_boundary_nodes),
        (time, max_time),
        temperature,
        method=method,
        events=[_terminal(conditions[end]) for end in ends],
        dense_output=True,
        **options,
    )

    if solution.status == -1:
        end = "solver_failure"
    elif solution.status == 1:
        end = next(end for end, times in zip(ends, solution.t_events, strict=True) if times.size)
    else:
        end = "max_time"

    end_time = solution.t[-1]
    sample_time = np.append(np.arange(time, end_time, save_every), end_time)
    sample_temperature = np.vstack([solution.sol(sample_time[:-1]).T, solution.y[:, -1]])

    return PhaseSolution(
        phase=phase, end=end, time=sample_time, temperature=sample_temperature, message=solution.message
    )


def run_procedure(model: Model) -> list[PhaseSolution]:
    """Integrates a model's phases in order from time zero.

    Args:
        model: The model to run, starting from its `initial_temperature`.

    Returns:
        The solution of each integrated phase, in order. A phase whose `until` condition is zero or negative at its
        start ends there with a single sample.
    """
    solutions = []
    time, temperature = 0.0, model.initial_temperature
    for phase in model.procedure.phases:
        if phase.until is not None and phase.until(time, temperature) <= 0:
            solution = PhaseSolution(
                phase=phase,
                end="until",
                time=np.array([time]),
                temperature=temperature[np.newaxis],
                message="The `until` condition was met at the phase's start.",
            )
        else:
            solution = integrate_phase(model, phase, time, temperature)
        solutions.append(solution)
        if solution.end != "until":
            break
        time, temperature = solution.time[-1], solution.temperature[-1]

    return solutions


def _terminal(
    condition: Callable[[float, NDArray[np.float64]], float],
) -> Callable[[float, NDArray[np.float64]], float]:
    def event(time: float, temperature: NDArray[np.float64]) -> float:
        return condition(time, temperature)

    event.terminal = True
    return event
