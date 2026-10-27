"""Integration of a model's phases through time."""

from collections.abc import Callable, Sequence
from functools import partial

import numpy as np
import xarray as xr
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

from pulse_tube_cooldown.core.errors import ConfigurationError
from pulse_tube_cooldown.core.types import Graph, Model, Phase, PhaseSolution, Result
from pulse_tube_cooldown.engine.equations import calc_link_heat_flows, calc_temperature_rate_of_change

__all__ = ["build_result", "integrate_phase", "reconstruct_heat_flow", "run_procedure"]


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

    graph = _graph(model, phase)
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


def reconstruct_heat_flow(graph: Graph, temperature: NDArray[np.float64]) -> NDArray[np.float64]:
    """Calculates the heat flow from each link into each node at sampled temperatures.

    Args:
        graph: The phase's nodes and active links.
        temperature: Temperature in K of each node, with one row per sample.

    Returns:
        Heat flow in W, indexed by sample, link, and node, with links and nodes in `graph` index order. Positive values
        indicate heat flow into the node; entries are NaN where a link does not touch a node.
    """
    return np.stack([calc_link_heat_flows(graph, sample) for sample in temperature])


def build_result(model: Model, solutions: Sequence[PhaseSolution]) -> Result:
    """Assembles the history of a run into a Result.

    Args:
        model: The simulated model.
        solutions: The solution of each integrated phase, in order, as returned by `run_procedure`.

    Returns:
        The run's Result, with the final solution's `end` as its status. The link axis holds every link in any phase
        of the procedure, in order of first appearance, and a phase switch keeps the samples on both sides.
    """
    links = {}
    for phase in model.procedure.phases:
        for link in phase.links:
            links.setdefault(link.name, link)
    link_names = list(links)
    node_names = [node.name for node in model.nodes]

    heat_flow = []
    for solution in solutions:
        graph = _graph(model, solution.phase)
        flow = np.full((solution.time.size, len(link_names), len(node_names)), np.nan)
        flow[:, [link_names.index(link.name) for link in graph.links]] = reconstruct_heat_flow(
            graph, solution.temperature
        )
        heat_flow.append(flow)

    connected = np.zeros((len(link_names), len(node_names)), dtype=bool)
    for row, link in enumerate(links.values()):
        connected[row, [node_names.index(name) for name in link.nodes]] = True
    status = solutions[-1].end

    data = xr.Dataset(
        {
            "phase": ("time", [solution.phase.name for solution in solutions for _ in solution.time]),
            "temperature": (
                ("time", "node"),
                np.concatenate([solution.temperature for solution in solutions]),
                {"units": "K"},
            ),
            "heat_flow": (("time", "link", "node"), np.concatenate(heat_flow), {"units": "W"}),
            "connected": (("link", "node"), connected),
        },
        coords={
            "time": ("time", np.concatenate([solution.time for solution in solutions]), {"units": "s"}),
            "node": node_names,
            "link": link_names,
        },
        attrs={"status": status, "message": solutions[-1].message},
    )

    return Result(status=status, model=model, data=data)


def _graph(model: Model, phase: Phase) -> Graph:
    return model.graphs[model.procedure.phases.index(phase)]


def _terminal(
    condition: Callable[[float, NDArray[np.float64]], float],
) -> Callable[[float, NDArray[np.float64]], float]:
    def event(time: float, temperature: NDArray[np.float64]) -> float:
        return condition(time, temperature)

    event.terminal = True
    return event
