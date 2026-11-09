"""The canonical PT1.5 head fitted to the test-stand load maps.

Run from the repository root with `python -m simulations.validation.test_stand_head.fit`. It writes `fit.txt` beside
this file and `data/heads/pt1_5/parameters.toml`.

The head's a1, a2, R1, R2, and stage coupling c are fitted jointly to every map in `data/test_stand/fit_set.csv`, each
run at its own compressor count and valve frequency, on the test-stand template `data/test_stand/fridge.toml` with its
bare-shield radiation link. At each point the stages are held at the measured T1 and T2 and the matrix nodes are
solved. The residual at a stage is the net heat flow into it, from the cooler, the radiation link, and the heater.

Each stage's residual is weighted by 1 / (T dQ/dT), the heat that moves the stage's log temperature by one, with
dQ/dT the pooled power law for the run's compressor count in `simulations/exploration/loadmaps/functional_forms.txt`.
Each map carries equal weight, whatever its number of points.
"""

import dataclasses
import functools
import sys
import tomllib
import warnings
from collections.abc import Callable, Sequence
from concurrent.futures import Executor, ProcessPoolExecutor
from pathlib import Path

import numpy as np
import scipy.optimize
import xarray as xr
from numpy.typing import NDArray

import pulse_tube_cooldown as ptc
from pulse_tube_cooldown.core.errors import SimulationError
from pulse_tube_cooldown.core.types import Graph
from pulse_tube_cooldown.engine.equations import calc_energy_balance_terms
from pulse_tube_cooldown.engine.steady_state import solve_steady_state
from pulse_tube_cooldown.library.links.cooler import Cooler, load_cooler
from pulse_tube_cooldown.measurements import read_load_map_set
from simulations.validation.test_stand_steady_state.coupling import analyze_coupling

sys.path.insert(0, str(Path(__file__).parents[2] / "exploration" / "loadmaps"))
from characterize import cell_partials, grids_of  # noqa: E402
from runs import load_runs  # noqa: E402

__all__ = ["calc_stage_residual", "fit_head"]

TEMPLATE = "data/test_stand/fridge.toml"
TWO_COMPRESSOR_DRIVE = "data/test_stand/pulse_tube/drive_two_compressors.toml"


def calc_stage_residual(
    points: xr.Dataset,
    *,
    packing: tuple[float, float],
    resistance: tuple[float, float],
    stage_coupling: float,
    template: str = TEMPLATE,
    two_compressor_drive: str = TWO_COMPRESSOR_DRIVE,
    workers: int = 2,
) -> xr.Dataset:
    """Calculates the net heat flow into each stage at each load-map point, with the matrix at steady state.

    Args:
        points: Q1 and Q2 in W, and T1 and T2 in K, along `point`, with each point's `compressors` and
            `valve_frequency` in Hz as coordinates, as `read_load_map_set` returns them.
        packing: Stage-1 and stage-2 packing coefficients a.
        resistance: Stage-1 and stage-2 needle-valve resistances R in Pa s/m^3.
        stage_coupling: The design's stage-coupling coefficient c.
        template: Test-stand fridge template, with one head on one compressor.
        two_compressor_drive: Drive file for points on two compressors.
        workers: Processes solving points in parallel.

    Returns:
        Q1 and Q2 along `point`: the net heat flow in W into each stage from the cooler, the radiation link, and the
        heater power.
    """
    with _executor(points, workers) as executor:
        uncoupled, coupling_basis, _ = _solve_points(
            executor, points, (*packing, *resistance), None, template, two_compressor_drive
        )
    residual = uncoupled + stage_coupling * np.column_stack([coupling_basis, np.zeros_like(coupling_basis)])
    return xr.Dataset(
        {name: ("point", column, {"units": "W"}) for name, column in zip(("Q1", "Q2"), residual.T, strict=True)},
        coords=points.coords,
    )


def fit_head(
    points: xr.Dataset,
    *,
    packing: tuple[float, float],
    resistance: tuple[float, float],
    template: str = TEMPLATE,
    two_compressor_drive: str = TWO_COMPRESSOR_DRIVE,
    workers: int = 2,
    relative_tolerance: float = 1e-4,
    verbose: int = 0,
) -> tuple[tuple[float, float], tuple[float, float], float]:
    """Fits a head's packing coefficients, resistances, and stage coupling to load-map points.

    The fit minimizes the weighted stage residuals of the module description, searching a and R from the given
    starting values. With the stages held, the stage coupling c adds c E_2 x (1 - x) to stage 1's residual and leaves
    the matrix unchanged, so at each a and R it takes its weighted least-squares value. Each point's matrix solve
    starts from its solution at the best parameters so far, so the fit follows the cold matrix branch that
    `calc_stage_residual` starts on. Points without a matrix steady state at the starting parameters are left out of
    the fit.

    Args:
        points: Load-map points, as `calc_stage_residual` takes them, with each point's map `serial` as a coordinate.
        packing: Starting stage-1 and stage-2 packing coefficients a.
        resistance: Starting stage-1 and stage-2 resistances R in Pa s/m^3.
        template: Test-stand fridge template, with one head on one compressor.
        two_compressor_drive: Drive file for points on two compressors.
        workers: Processes solving points in parallel.
        relative_tolerance: Relative change in the parameters, and in the weighted sum of squares, at which the fit
            stops.
        verbose: Progress the optimizer prints: 0 for none, 1 at the end, 2 at each iteration.

    Returns:
        The fitted packing coefficients, resistances in Pa s/m^3, and stage coupling.
    """
    start = np.array([*packing, *resistance])

    with _executor(points, workers) as executor:
        uncoupled, _, states = _solve_points(executor, points, tuple(start), None, template, two_compressor_drive)
        solved = ~np.isnan(uncoupled).any(axis=1)
        points = points.isel(point=solved)
        best_cost, best_states = np.inf, [state for state, kept in zip(states, solved, strict=True) if kept]
        _, map_index, map_size = np.unique(points["serial"].to_numpy(), return_inverse=True, return_counts=True)
        stage_weight = 1 / (np.sqrt(map_size[map_index])[:, None] * _log_temperature_scale(points))

        def evaluate(log_ratio: NDArray[np.float64]) -> tuple[NDArray[np.float64], float]:
            """The weighted residuals at the best stage coupling, and that coupling."""
            nonlocal best_cost, best_states
            uncoupled, coupling_basis, states = _solve_points(
                executor, points, tuple(start * np.exp(log_ratio)), best_states, template, two_compressor_drive
            )
            stage_1 = uncoupled[:, 0] * stage_weight[:, 0]
            basis = coupling_basis * stage_weight[:, 0]
            stage_coupling = -np.sum(stage_1 * basis) / np.sum(basis**2)
            weighted = np.concatenate([stage_1 + stage_coupling * basis, uncoupled[:, 1] * stage_weight[:, 1]])
            if np.sum(weighted**2) < best_cost:
                best_cost, best_states = np.sum(weighted**2), states
            return weighted, stage_coupling

        result = scipy.optimize.least_squares(
            lambda log_ratio: evaluate(log_ratio)[0],
            np.zeros(len(start)),
            x_scale=0.1,
            diff_step=1e-3,
            xtol=relative_tolerance,
            ftol=relative_tolerance,
            verbose=verbose,
        )
        _, stage_coupling = evaluate(result.x)
    a1, a2, r1, r2 = start * np.exp(result.x)
    return (float(a1), float(a2)), (float(r1), float(r2)), float(stage_coupling)


def _log_temperature_scale(points: xr.Dataset) -> NDArray[np.float64]:
    """Heat in W that moves each stage's log temperature by one at each point: T dQ/dT from the pooled power laws."""
    t1, t2 = points["T1"].to_numpy(), points["T2"].to_numpy()
    two = points["compressors"].to_numpy() == 2
    stage_1 = np.where(two, 102.1 * t1**-0.08, 70.69 * t1**-0.01)
    stage_2 = np.where(two, 0.1958 * t2**2.34, 0.1004 * t2**2.38)
    return np.column_stack([stage_1, stage_2])


def _executor(points: xr.Dataset, workers: int) -> Executor:
    """A process pool for solving points in parallel."""
    return ProcessPoolExecutor(max_workers=min(workers, points.sizes["point"]))


def _solve_points(
    executor: Executor,
    points: xr.Dataset,
    parameters: tuple[float, ...],
    states: Sequence[NDArray[np.float64] | None] | None,
    template: str,
    two_compressor_drive: str,
) -> tuple[NDArray[np.float64], NDArray[np.float64], list[NDArray[np.float64]]]:
    """Solves every point in parallel at a1, a2, R1, and R2, as `_solve_point` does, stacking its results."""
    count = points.sizes["point"]
    tasks = [
        (
            int(points["compressors"][index]),
            float(points["valve_frequency"][index]),
            (float(points["T1"][index]), float(points["T2"][index])),
            (float(points["Q1"][index]), float(points["Q2"][index])),
            None if states is None else states[index],
        )
        for index in range(count)
    ]
    solve = functools.partial(
        _solve_point, parameters=parameters, template=template, two_compressor_drive=two_compressor_drive
    )
    solved = list(executor.map(solve, *zip(*tasks, strict=True)))
    uncoupled, coupling_basis, states = zip(*solved, strict=True)
    return np.array(uncoupled), np.array(coupling_basis), list(states)


def _solve_point(
    compressors: int,
    frequency: float,
    stage_temperature: tuple[float, float],
    heater_power: tuple[float, float],
    state: NDArray[np.float64] | None,
    *,
    parameters: tuple[float, ...],
    template: str,
    two_compressor_drive: str,
) -> tuple[NDArray[np.float64], float, NDArray[np.float64]]:
    """Solves one point's matrix with the stages held, at a1, a2, R1, and R2 and no stage coupling.

    Without a starting state, the solve starts on the cold matrix branch with stage 2 at 2.3 K and steps stage 2 up
    to its measured temperature by at most 0.1 K at a time.

    Returns:
        The stage residuals in W, the stage-1 coupling term per unit c, E_2 x (1 - x) in W, and the steady state.
        Where no steady state is found, the residuals and term are NaN and the state is the starting state.
    """
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    graph, boundary, stages, profile = _point_graph(
        compressors, frequency, parameters, 0.0, template, two_compressor_drive
    )
    t1, t2 = stage_temperature
    if state is None:
        cold_start = min(t2, 2.3)
        start = profile(t1, cold_start)
        steps = np.linspace(cold_start, t2, int(np.ceil((t2 - cold_start) / 0.1)) + 1)
    else:
        start = state
        steps = np.array([t2])
    solved = start.copy()
    try:
        for stage_2 in steps:
            solved[list(stages)] = (t1, stage_2)
            solved = solve_steady_state(graph, solved, fixed_boundary_nodes=(*boundary, *stages), tolerance=1e-6)
    except SimulationError:
        return np.full(2, np.nan), np.nan, start
    heat_flow = calc_energy_balance_terms(graph, solved).heat_flow[list(stages)]
    coupled, _, _, _ = _point_graph(compressors, frequency, parameters, 1.0, template, two_compressor_drive)
    coupling_basis = calc_energy_balance_terms(coupled, solved).heat_flow[stages[0]] - heat_flow[0]
    return heat_flow + heater_power, float(coupling_basis), solved


def _cell_coupling(
    compressors: int,
    frequency: float,
    stage_temperature: tuple[float, float],
    *,
    parameters: tuple[float, ...],
    stage_coupling: float,
    template: str,
    two_compressor_drive: str,
) -> NDArray[np.float64]:
    """dQ_i/dT_j in W/K, heat into stage i, at held stage temperatures, or NaN without a matrix steady state.

    Args:
        compressors: Compressor count.
        frequency: Valve frequency in Hz.
        stage_temperature: T1 and T2 in K.
        parameters: a1, a2, R1, and R2.
        stage_coupling: The stage coupling c.
        template: Test-stand fridge template.
        two_compressor_drive: Drive file for two compressors.
    """
    residual, _, state = _solve_point(
        compressors,
        frequency,
        stage_temperature,
        (0.0, 0.0),
        None,
        parameters=parameters,
        template=template,
        two_compressor_drive=two_compressor_drive,
    )
    if np.isnan(residual).any():
        return np.full((2, 2), np.nan)
    graph, boundary, _, _ = _point_graph(
        compressors, frequency, parameters, stage_coupling, template, two_compressor_drive
    )
    coupling = analyze_coupling(
        graph, state, stages=("Stage 1", "Stage 2"), fixed_boundary_nodes=boundary, tolerance=1e-6
    )
    return coupling.stage_derivative


def _point_graph(
    compressors: int,
    frequency: float,
    parameters: tuple[float, ...],
    stage_coupling: float,
    template: str,
    two_compressor_drive: str,
) -> tuple[Graph, tuple[int, ...], tuple[int, int], Callable[[float, float], NDArray[np.float64]]]:
    """The template's graph for a run, head parameters a1, a2, R1, and R2, and stage coupling c.

    Returns:
        The graph, the template's fixed boundary nodes, the stage indices, and a starting profile in K of every node
        for given T1 and T2.
    """
    graph, boundary, stages, profile = _template(template, two_compressor_drive, compressors)
    a1, a2, r1, r2 = parameters
    return (
        _with_head(graph, packing=(a1, a2), resistance=(r1, r2), frequency=frequency, stage_coupling=stage_coupling),
        boundary,
        stages,
        profile,
    )


@functools.cache
def _template(
    template: str, two_compressor_drive: str, compressors: int
) -> tuple[Graph, tuple[int, ...], tuple[int, int], Callable[[float, float], NDArray[np.float64]]]:
    """The template's graph on a compressor count, with its fixed boundary nodes, stage indices, and starting profile.

    The profile places each matrix node between its section's warm and cold ends at the fraction, in log
    temperature, that it holds on the cold matrix branch of `stage2_gate.txt` at T1 = 50 K and T2 = 3.4 K.
    """
    model = ptc.configure_model(template)
    graph = model.graphs[0]
    if compressors == 2:
        with open(template, "rb") as file:
            (config,) = (link for link in tomllib.load(file)["fridge"]["links"] if link["type"] == "cooler")
        cooler = load_cooler({**config, "drive": two_compressor_drive})
        graph = dataclasses.replace(
            graph, links=tuple(cooler if isinstance(link, Cooler) else link for link in graph.links)
        )
    first_stage, second_stage = (graph.get_node_index(name) for name in ("Stage 1", "Stage 2"))
    warm = float(model.initial_temperature[model.fixed_boundary_nodes[0]])
    (cooler,) = (link for link in graph.links if isinstance(link, Cooler))
    first, second = ([graph.get_node_index(name) for name in section] for section in cooler.heads[0].matrix_nodes)
    gate_first, gate_second = (217, 93), (45, 26, 13, 4.7, 3.9, 3.7, 3.6, 3.55, 3.5, 3.45, 3.42)
    first_fraction = np.log(np.divide(gate_first, 50.0)) / np.log(warm / 50.0)
    second_fraction = np.log(np.divide(gate_second, 3.4)) / np.log(50.0 / 3.4)

    def profile(t1: float, t2: float) -> NDArray[np.float64]:
        state = model.initial_temperature.copy()
        state[first] = t1 * (warm / t1) ** first_fraction
        state[second] = t2 * (t1 / t2) ** second_fraction
        return state

    return graph, model.fixed_boundary_nodes, (first_stage, second_stage), profile


def _with_head(
    graph: Graph,
    *,
    packing: tuple[float, float],
    resistance: tuple[float, float],
    frequency: float,
    stage_coupling: float,
) -> Graph:
    """The graph with every head's packing, resistance, valve frequency, and stage coupling replaced."""

    def replace(link):
        if not isinstance(link, Cooler):
            return link
        heads = tuple(
            dataclasses.replace(
                mounted,
                head=dataclasses.replace(
                    mounted.head,
                    model=dataclasses.replace(mounted.head.model, stage_coupling=stage_coupling),
                    packing=packing,
                    resistance=resistance,
                    frequency=frequency,
                ),
            )
            for mounted in link.heads
        )
        return dataclasses.replace(link, heads=heads)

    return dataclasses.replace(graph, links=tuple(replace(link) for link in graph.links))


def _grid_cell_couplings(
    executor: Executor,
    parameters: tuple[float, ...],
    stage_coupling: float,
    template: str,
    two_compressor_drive: str,
) -> list[dict]:
    """Measured and model couplings at each grid cell of `07_capacity_partials.png` on a known compressor count.

    Each cell carries its compressor count, mean T1 and T2 in K, and the measured and model dQ_i/dT_j in W/K as 2x2
    arrays, with heat into the stage positive. The model runs at the cell's valve frequency, or 1.0 Hz where unknown.
    """
    cells = [cell for run in grids_of(load_runs()) if run.compressors in (1, 2) for cell in cell_partials(run)]
    couple = functools.partial(
        _cell_coupling,
        parameters=parameters,
        stage_coupling=stage_coupling,
        template=template,
        two_compressor_drive=two_compressor_drive,
    )
    model = executor.map(
        couple,
        [cell["run"].compressors for cell in cells],
        [cell["run"].valve_frequency or 1.0 for cell in cells],
        [(cell["t1"], cell["t2"]) for cell in cells],
    )
    return [
        {
            "compressors": cell["run"].compressors,
            "t1": cell["t1"],
            "t2": cell["t2"],
            "measured": -np.array([[cell["dq1_dt1"], cell["dq1_dt2"]], [cell["dq2_dt1"], cell["dq2_dt2"]]]),
            "model": derivative,
        }
        for cell, derivative in zip(cells, model, strict=True)
    ]


def _report(
    points: xr.Dataset,
    residual: xr.Dataset,
    cells: list[dict],
    *,
    start: tuple[float, ...],
    fitted: tuple[float, ...],
) -> str:
    """The fit's text report: parameters, each run's residual, each point's residual, and the coupling comparison.

    `start` holds the starting a1, a2, R1, and R2, and `fitted` those and c.
    """
    scale = _log_temperature_scale(points)
    stage_residual = np.column_stack([residual["Q1"].to_numpy(), residual["Q2"].to_numpy()])
    solved = ~np.isnan(stage_residual).any(axis=1)
    names = ("a1", "a2", "R1 (Pa s/m^3)", "R2 (Pa s/m^3)", "c")
    lines = [
        "The canonical PT1.5 head fitted to data/test_stand/fit_set.csv, written by fit.py beside this",
        "file. Residuals are the net heat flow into each stage in W, heater power included, with the matrix at steady",
        "state and the stages held at the measured T1 and T2. Weighted residuals divide each by T dQ/dT, the heat that",
        "moves the stage's log temperature by one.",
        "",
        f"{'':14} {'start':>10} {'fitted':>10}",
        *(f"{name:14} {s:10.4g} {f:10.4g}" for name, s, f in zip(names, start, fitted, strict=False)),
        f"{names[-1]:14} {'solved':>10} {fitted[-1]:10.4g}",
        "",
        "Each run's rms residual in W over its solved points, and in parentheses its weighted rms. The fit leaves out",
        "points without a matrix steady state above the helium floor at its starting parameters.",
        "",
        f"{'serial':7} {'date':10} {'comp.':>5} {'f (Hz)':>6} {'points':>6} {'stage 1':>16} {'stage 2':>16}",
    ]
    serials = list(dict.fromkeys(points["serial"].to_numpy()))
    for serial in [s for s in serials if s != "S0028"] + [s for s in serials if s == "S0028"]:
        if serial == "S0028":
            lines += ["", "S0028, on two compressors, reported separately:"]
        mask = points["serial"].to_numpy() == serial
        first = np.flatnonzero(mask)[0]
        rms = np.sqrt(np.nanmean(stage_residual[mask] ** 2, axis=0))
        weighted = np.sqrt(np.nanmean((stage_residual[mask] / scale[mask]) ** 2, axis=0))
        lines.append(
            f"{serial:7} {str(points['date'][first].item()):10} {int(points['compressors'][first]):5d} "
            f"{float(points['valve_frequency'][first]):6.1f} {int((mask & solved).sum()):6d} "
            f"{rms[0]:7.3f} ({weighted[0]:6.4f}) {rms[1]:7.3f} ({weighted[1]:6.4f})"
        )
    lines += [
        "",
        "Each point's heater powers P1 and P2 in W, stage temperatures in K, and residuals in W. NaN marks a point",
        "without a matrix steady state above the helium floor at the fitted parameters.",
        "",
        f"{'serial':7} {'P1':>5} {'P2':>4} {'T1':>6} {'T2':>5} {'stage 1':>8} {'stage 2':>8}",
        *(
            f"{serial:7} {p1:5.1f} {p2:4.1f} {t1:6.2f} {t2:5.2f} {r1:+8.3f} {r2:+8.3f}"
            for serial, p1, p2, t1, t2, (r1, r2) in zip(
                points["serial"].to_numpy(),
                points["Q1"].to_numpy(),
                points["Q2"].to_numpy(),
                points["T1"].to_numpy(),
                points["T2"].to_numpy(),
                stage_residual,
                strict=True,
            )
        ),
        "",
        "Couplings dQi/dTj in W/K, heat into the stage positive, at figure 07's grid cells: the median of each T2",
        "bin's cells as measured / model. The model is solved at each cell's mean T1 and T2.",
    ]
    for compressors in (1, 2):
        group = [cell for cell in cells if cell["compressors"] == compressors]
        missing = sum(np.isnan(cell["model"]).any() for cell in group)
        lines += [
            "",
            f"{compressors} compressor{'s' if compressors == 2 else ''}: {len(group)} cells, {missing} without a model "
            "steady state",
            f"{'T2 (K)':9} {'cells':>5} {'dQ1/dT1':>15} {'dQ1/dT2':>15} {'dQ2/dT1':>15} {'dQ2/dT2':>15}",
        ]
        for low in np.arange(2.0, 5.0, 0.5):
            binned = [cell for cell in group if low <= cell["t2"] < low + 0.5]
            if not binned:
                continue
            measured = np.median([cell["measured"] for cell in binned], axis=0)
            model = np.nanmedian([cell["model"] for cell in binned], axis=0)
            entries = " ".join(
                f"{measured[i, j]:+7.3f}/{model[i, j]:+7.3f}" for i, j in ((0, 0), (0, 1), (1, 0), (1, 1))
            )
            lines.append(f"{low:3.1f}-{low + 0.5:3.1f}   {len(binned):5d} {entries}")
    return "\n".join(lines) + "\n"


def main(
    packing: tuple[float, float] = (0.05, 0.02),
    resistance: tuple[float, float] = (1e10, 4.5e10),
    fit_set: str = "data/test_stand/fit_set.csv",
    parameters_file: str = "data/heads/pt1_5/parameters.toml",
    frequency: float = 1.0,
    workers: int = 2,
) -> None:
    """Fits the head, writes its parameter file, and writes `fit.txt` beside this file.

    Args:
        packing: Starting stage-1 and stage-2 packing coefficients.
        resistance: Starting stage-1 and stage-2 resistances in Pa s/m^3.
        fit_set: The load-map set table.
        parameters_file: The head parameter file written with the fitted a and R.
        frequency: Valve frequency in Hz written to the parameter file.
        workers: Processes solving points in parallel.
    """
    points = read_load_map_set(fit_set)
    fitted_packing, fitted_resistance, fitted_coupling = fit_head(
        points, packing=packing, resistance=resistance, workers=workers, verbose=2
    )
    residual = calc_stage_residual(
        points,
        packing=fitted_packing,
        resistance=fitted_resistance,
        stage_coupling=fitted_coupling,
        workers=workers,
    )
    with ProcessPoolExecutor(max_workers=workers) as executor:
        cells = _grid_cell_couplings(
            executor, (*fitted_packing, *fitted_resistance), fitted_coupling, TEMPLATE, TWO_COMPRESSOR_DRIVE
        )

    report = Path(__file__).with_name("fit.txt")
    source = f"{report.resolve().relative_to(Path.cwd().resolve()).as_posix()}: fitted to {fit_set}"
    Path(parameters_file).write_text(
        f"packing = [{fitted_packing[0]:.4g}, {fitted_packing[1]:.4g}]\n"
        f"resistance = [{fitted_resistance[0]:.4g}, {fitted_resistance[1]:.4g}]  # Pa s/m^3\n"
        f"frequency = {frequency}  # Hz\n"
        "\n[sources]\n"
        f'packing = "{source}"\n'
        f'resistance = "{source}"\n'
        "frequency = \"data/test_stand/fit_set.csv: the standard runs' 1.0 Hz, such as S0009's\"\n",
        encoding="utf-8",
    )
    report.write_text(
        _report(
            points,
            residual,
            cells,
            start=(*packing, *resistance),
            fitted=(*fitted_packing, *fitted_resistance, fitted_coupling),
        ),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
