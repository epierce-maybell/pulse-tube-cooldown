"""Run, jointly calibrate, and plot the five-coefficient candidate experiment."""

import argparse
import hashlib
import json
import os
import time
import tomllib
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from numpy.typing import NDArray
from scipy.integrate import BDF, Radau
from scipy.optimize import least_squares, root

from pulse_tube_cooldown.core.types import Node
from pulse_tube_cooldown.library.links.radiation import Radiation, load_radiation
from pulse_tube_cooldown.library.materials.load import load_material
from pulse_tube_cooldown.measurements import read_cooldown
from simulations.exploration.pulse_tube_candidate.model import (
    Candidate,
    Parameters,
    ProfileSolveError,
    PropertyLimitError,
    build_candidate,
)
from simulations.validation.head_comparison.comparison import (
    HeadData,
    Prediction,
    cell_derivatives,
    comparison_style,
    load_population,
    plot_comparison,
)


@dataclass(frozen=True, slots=True)
class Stand:
    """The configured bare-head stage inventories and shield radiation."""

    stages: tuple[Node, Node]
    radiation: Radiation

    @classmethod
    def load(cls, path: Path) -> "Stand":
        """Reads the existing template without constructing the retired cooler."""
        with path.open("rb") as stream:
            config = tomllib.load(stream)
        nodes = {
            spec["name"]: Node(
                spec["name"], tuple((load_material(material), mass) for material, mass in spec["composition"])
            )
            for spec in config["fridge"]["nodes"]
            if spec["name"] in ("Stage 1", "Stage 2")
        }
        radiation = next(load_radiation(spec) for spec in config["fridge"]["links"] if spec["type"] == "radiation")
        return cls((nodes["Stage 1"], nodes["Stage 2"]), radiation)

    def external_heat(self, temperature: NDArray, warm: float) -> NDArray:
        """Returns known stage heat inputs in W, excluding applied heaters."""
        return np.array([self.radiation.calc_energy_balance_terms(np.array([temperature[0], warm])).heat_flow[0], 0.0])


def initial_temperature(head: HeadData) -> float:
    """Uses the last finite pre-start TPR reading as the room proxy."""
    raw = read_cooldown(head.cooldown.source)
    with (Path(__file__).resolve().parents[2] / "validation/head_comparison/records.toml").open("rb") as stream:
        records = tomllib.load(stream)
    start = np.datetime64(next(row["cooldown"]["start"] for row in records["heads"] if row["serial"] == head.serial))
    pre = raw.TPR.where(raw.time < start, drop=True).dropna("time")
    return float(pre[-1] if pre.size else raw.TPR.dropna("time")[0])


def simulate(
    candidate: Candidate,
    stand: Stand,
    head: HeadData,
    parameters: Parameters,
    *,
    warm: float,
) -> tuple[xr.Dataset, dict]:
    """Integrates seven matrix energies and two plate log temperatures with BDF."""
    uniform = np.full(8, warm)
    boundary = np.full(3, warm)
    energy = candidate.energies(uniform, boundary)
    state = np.r_[energy, 0.0, 0.0]
    last_profile = uniform.copy()
    maximum_clamp_residual = 0.0
    maximum_inversion_error = 0.0
    rhs_calls = 0

    def rhs(_time: float, value: NDArray) -> NDArray:
        nonlocal last_profile, maximum_clamp_residual, maximum_inversion_error, rhs_calls
        stage = warm * np.exp(value[7:])
        boundary = np.r_[warm, stage]
        centers, evaluated = candidate.recover_profile(
            value[:7],
            boundary,
            parameters,
            guess=last_profile,
        )
        last_profile = centers
        maximum_clamp_residual = max(maximum_clamp_residual, abs(evaluated.flows[6] - evaluated.flows[7]))
        maximum_inversion_error = max(
            maximum_inversion_error, float(np.max(np.abs(candidate.energies(centers, boundary) - value[:7])))
        )
        heat = evaluated.heat_flow[[0, 2]] + stand.external_heat(stage, warm)
        capacity = np.array([node.calc_heat_capacity(float(t)) for node, t in zip(stand.stages, stage, strict=True)])
        rhs_calls += 1
        return np.r_[evaluated.energy_rate, heat / (capacity * stage)]

    final_time = float(head.cooldown.data.time[-1])
    numerics = candidate.numerics
    atol = np.r_[
        np.full(7, numerics["energy_absolute_tolerance"]), np.full(2, numerics["temperature_absolute_tolerance"] / warm)
    ]
    samples = np.r_[np.arange(0, final_time, numerics["sample_interval"]), final_time]
    states = [state.copy()]
    sample_index = 1
    solver_type = {"BDF": BDF, "Radau": Radau}[numerics["solver"]]
    solver = solver_type(
        rhs, 0.0, state, final_time, rtol=numerics["relative_tolerance"], atol=atol, max_step=numerics["maximum_step"]
    )
    rejected_domain_steps = 0
    integrator_fallbacks = 0
    started = time.perf_counter()
    accepted_time, accepted_state = 0.0, state.copy()
    while solver.status == "running":
        attempted_step = solver.h_abs
        try:
            message = solver.step()
        except (PropertyLimitError, ProfileSolveError) as error:
            rejected_domain_steps += 1
            next_step = min(attempted_step, numerics["maximum_step"]) / 4
            if next_step < 1e-4:
                raise ProfileSolveError(
                    f"Integration stopped at {accepted_time:.6g} s after bounded step retries: {error}"
                ) from error
            last_profile = candidate.initial_centers(np.r_[warm, warm * np.exp(accepted_state[7:])])
            solver = solver_type(
                rhs,
                accepted_time,
                accepted_state,
                final_time,
                rtol=numerics["relative_tolerance"],
                atol=atol,
                first_step=next_step,
                max_step=numerics["maximum_step"],
            )
            continue
        if solver.status == "failed":
            if solver_type is BDF:
                integrator_fallbacks += 1
                solver_type = Radau
                solver = Radau(
                    rhs,
                    accepted_time,
                    accepted_state,
                    final_time,
                    rtol=numerics["relative_tolerance"],
                    atol=atol,
                    max_step=numerics["maximum_step"],
                )
                continue
            raise ProfileSolveError(
                f"{solver_type.__name__} failed at {solver.t:.6g} s; "
                f"stages {(warm * np.exp(accepted_state[7:])).tolist()} K; {message}"
            )
        dense = solver.dense_output()
        while sample_index < len(samples) and samples[sample_index] <= solver.t:
            states.append(dense(samples[sample_index]))
            sample_index += 1
        accepted_time, accepted_state = solver.t, solver.y.copy()
    states = np.asarray(states)
    data = xr.Dataset(
        {
            **{name: ("time", warm * np.exp(states[:, 7 + i]), {"units": "K"}) for i, name in enumerate(("T1", "T2"))},
            "matrix_energy": (("time", "piece"), states[:, :7], {"units": "J"}),
        },
        coords={"time": ("time", samples, {"units": "s"}), "piece": np.arange(1, 8)},
    )
    diagnostics = {
        "seconds": time.perf_counter() - started,
        "rhs_calls": rhs_calls,
        "domain_step_retries": rejected_domain_steps,
        "integrator_fallbacks": integrator_fallbacks,
        "maximum_clamp_residual_W": float(maximum_clamp_residual),
        "maximum_energy_inversion_error_J": maximum_inversion_error,
    }
    print(
        f"{head.serial} cooldown: {diagnostics['seconds']:.2f} s, {rhs_calls} RHS calls, "
        f"{rejected_domain_steps} domain retries; final stages "
        f"{data.T1[-1].item():.3f}, {data.T2[-1].item():.3f} K",
        flush=True,
    )
    return data, diagnostics


def solve_map(
    candidate: Candidate,
    stand: Stand,
    head: HeadData,
    parameters: Parameters,
    *,
    warm: float,
) -> tuple[xr.Dataset, dict]:
    """Solves all material and plate balances at each measured heater pair."""
    measured = head.load_map.data
    predicted = []
    residuals = []
    previous = None
    for index in range(measured.sizes["point"]):
        heater = np.array([measured[name][index].item() for name in ("Q1", "Q2")])
        stage_guess = np.array([measured[name][index].item() for name in ("T1", "T2")])
        centers = candidate.initial_centers(np.r_[warm, stage_guess])
        initial = np.r_[candidate._encode(centers), candidate._encode(stage_guess)]

        def residual(coordinate: NDArray, heater: NDArray = heater) -> NDArray:
            centers, stage = candidate._decode(coordinate[:8]), candidate._decode(coordinate[8:])
            evaluated = candidate.evaluate(
                centers,
                np.r_[warm, stage],
                parameters,
            )
            stage_heat = evaluated.heat_flow[[0, 2]] + stand.external_heat(stage, warm) + heater
            return np.r_[evaluated.energy_rate, evaluated.flows[6] - evaluated.flows[7], stage_heat]

        solution = root(residual, initial if previous is None else previous, tol=1e-9)
        error = float(np.max(np.abs(residual(solution.x))))
        if error > candidate.numerics["equilibrium_flow_tolerance"]:
            solution = least_squares(residual, initial, max_nfev=200, xtol=1e-11, ftol=1e-11, gtol=1e-11)
            error = float(np.max(np.abs(residual(solution.x))))
        if error > candidate.numerics["equilibrium_flow_tolerance"]:
            raise ProfileSolveError(f"Map point {index}, heaters {heater}: balance residual {error:.6g} W.")
        previous = solution.x
        predicted.append(candidate._decode(solution.x[8:]))
        residuals.append(error)
    predicted = np.array(predicted)
    result = measured[["Q1", "Q2"]].copy()
    for i, name in enumerate(("T1", "T2")):
        result[name] = ("point", predicted[:, i], {"units": "K"})
    return result, {"maximum_balance_residual_W": float(max(residuals))}


def predict_one(candidate: Candidate, stand: Stand, head: HeadData, parameters: Parameters, *, warm: float):
    """Calculates both curves and a fully equilibrated map using one parameter set."""
    load_map, map_diagnostics = solve_map(candidate, stand, head, parameters, warm=warm)
    cooldown, cooldown_diagnostics = simulate(candidate, stand, head, parameters, warm=warm)
    return Prediction(cooldown=cooldown, load_map=load_map), {"cooldown": cooldown_diagnostics, "map": map_diagnostics}


def residual_blocks(head: HeadData, prediction: Prediction, config: dict) -> dict[str, NDArray]:
    """Returns eight equally weighted blocks, including four fitted Jacobians."""
    blocks = {}
    for name in ("T1", "T2"):
        measured = head.cooldown.data[name]
        valid = np.isfinite(measured.values)
        if name == "T1":
            valid &= measured.values < head.stage_1_max_temperature
        t = head.cooldown.data.time.values[valid]
        observed = measured.values[valid]
        estimate = np.interp(t, prediction.cooldown.time, prediction.cooldown[name])
        # Time quadrature prevents logger density from changing the fitting objective.
        time_weight = np.gradient(t)
        blocks[f"cooldown_{name}"] = np.log(estimate / observed) * np.sqrt(time_weight / np.sum(time_weight))
        observed_map = head.load_map.data[name].values
        blocks[f"map_{name}"] = (prediction.load_map[name].values / observed_map - 1) / np.sqrt(len(observed_map))
    observed_cells = cell_derivatives(head.load_map.data)
    predicted_cells = cell_derivatives(prediction.load_map)
    for index, name in enumerate(("dQ1dT1", "dQ1dT2", "dQ2dT1", "dQ2dT2")):
        observed = observed_cells[name].values
        scale = max(float(np.sqrt(np.mean(observed**2))), config["fit"]["jacobian_scale_floor"][index])
        blocks[name] = (predicted_cells[name].values - observed) / (scale * np.sqrt(len(observed)))
    return blocks


def transformed_parameters(values: NDArray) -> Parameters:
    """Uses logarithms for positive parameters and log(1+a) for nonnegative a."""
    result = np.exp(values)
    result[1] -= 1
    return Parameters(*result)


def parameter_coordinates(values: NDArray) -> NDArray:
    """Maps the five physical coefficients to optimization coordinates."""
    values = np.array(values, dtype=float)
    values[1] += 1
    return np.log(values)


def fit_one(candidate: Candidate, stand: Stand, head: HeadData, config: dict, *, destination: Path):
    """Jointly fits cooldown, map, and Jacobians, saving every evaluation."""
    warm = initial_temperature(head)
    config_fit = config["fit"]
    if config_fit.get("resume", False) and (destination / f"{head.serial}_best.json").exists():
        previous = json.loads((destination / f"{head.serial}_best.json").read_text(encoding="utf-8"))
        config_fit = {**config_fit, "initial": previous["parameters"]}
    start = parameter_coordinates(config_fit["initial"])
    best = None
    log_path = destination / f"{head.serial}_evaluations.jsonl"
    evaluation_number = 0
    if log_path.exists():
        evaluation_number = max(json.loads(line)["evaluation"] for line in log_path.read_text().splitlines())
    source_paths = [
        Path(__file__).with_name(name)
        for name in (
            "model.py",
            "properties.py",
            "run.py",
            "settings.toml",
        )
    ]
    source_paths += [
        Path(__file__).resolve().parents[2] / "validation/head_comparison/records.toml",
        head.cooldown.source,
        head.load_map.source,
        Path(__file__).resolve().parents[3] / "data/test_stand/test_stand.toml",
    ]
    manifest = {
        "starting_evaluation": evaluation_number + 1,
        "config": config,
        "starting_parameters": config_fit["initial"],
        "warm_temperature_K": warm,
        "input_hashes": {str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in source_paths},
    }
    with (destination / f"{head.serial}_runs.jsonl").open("a", encoding="utf-8") as stream:
        stream.write(json.dumps(manifest) + "\n")

    def objective(coordinate: NDArray) -> NDArray:
        nonlocal best, evaluation_number
        evaluation_number += 1
        parameters = transformed_parameters(coordinate)
        started = time.perf_counter()
        try:
            prediction, diagnostics = predict_one(candidate, stand, head, parameters, warm=warm)
            blocks = residual_blocks(head, prediction, config)
            residual = np.concatenate(list(blocks.values()))
            score = float(residual @ residual)
            record = {
                "evaluation": evaluation_number,
                "parameters": parameters.as_array().tolist(),
                "score": score,
                "blocks": {k: float(v @ v) for k, v in blocks.items()},
                "diagnostics": diagnostics,
            }
            if best is None or score < best[0]:
                best = (score, parameters, prediction, diagnostics, record["blocks"])
                write_prediction(destination, head.serial, prediction)
                (destination / f"{head.serial}_best.json").write_text(
                    json.dumps(
                        {
                            **record,
                            "warm_temperature_K": warm,
                        },
                        indent=2,
                    )
                    + "\n",
                    encoding="utf-8",
                )
        except (ProfileSolveError, PropertyLimitError) as error:
            measured_count = sum(
                np.isfinite(head.cooldown.data[name].values).sum()
                if name == "T2"
                else np.sum(
                    np.isfinite(head.cooldown.data[name].values)
                    & (head.cooldown.data[name].values < head.stage_1_max_temperature)
                )
                for name in ("T1", "T2")
            )
            count = (
                measured_count
                + 2 * head.load_map.data.sizes["point"]
                + 4 * cell_derivatives(head.load_map.data).sizes["cell"]
            )
            residual = np.full(count, 100.0 + np.linalg.norm(coordinate - start))
            record = {
                "evaluation": evaluation_number,
                "parameters": parameters.as_array().tolist(),
                "error": str(error),
                "score": float(residual @ residual),
            }
        record["seconds"] = time.perf_counter() - started
        with log_path.open("a", encoding="utf-8") as stream:
            stream.write(json.dumps(record) + "\n")
        print(
            f"{head.serial} evaluation {evaluation_number}: {record['score']:.6g}, "
            f"{record['seconds']:.2f} s" + (f"; {record['error']}" if "error" in record else ""),
            flush=True,
        )
        return residual

    result = least_squares(
        objective,
        start,
        bounds=(parameter_coordinates(config_fit["lower"]), parameter_coordinates(config_fit["upper"])),
        diff_step=config_fit["relative_step"],
        max_nfev=config_fit["maximum_evaluations"],
        xtol=1e-4,
        ftol=1e-4,
        gtol=1e-4,
        verbose=1,
    )
    if best is None:
        raise ProfileSolveError(f"No supported complete prediction was obtained for {head.serial}.")
    (destination / f"{head.serial}_optimizer.json").write_text(
        json.dumps(
            {
                "message": result.message,
                "success": bool(result.success),
                "evaluations": int(result.nfev),
                "optimality": float(result.optimality),
                "active_mask": result.active_mask.tolist(),
                "jacobian_singular_values": np.linalg.svd(result.jac, compute_uv=False).tolist(),
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    return best[2]


def write_prediction(destination: Path, serial: str, prediction: Prediction) -> None:
    """Writes the exact curves and equilibrium points supplied to the plotter."""
    for suffix, data in (("cooldown", prediction.cooldown), ("map", prediction.load_map)):
        target = destination / f"{serial}_{suffix}.nc"
        temporary = target.with_suffix(".next.nc")
        data.to_netcdf(temporary)
        temporary.replace(target)


def fit_worker(head: HeadData, config: dict, destination: Path) -> tuple[str, Prediction]:
    """Fits one head in a process with its own numerical tables."""
    with warnings.catch_warnings():
        warnings.simplefilter("once", RuntimeWarning)
        candidate = build_candidate(config)
        stand = Stand.load(Path(__file__).resolve().parents[3] / "data/test_stand/test_stand.toml")
        return head.serial, fit_one(candidate, stand, head, config, destination=destination)


def main() -> None:
    """Runs the candidate experiment serially without acceptance-test collection."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--head", action="append", help="Head serial; default is all six comparison records")
    parser.add_argument("--fit", action="store_true", help="Jointly fit all eight residual blocks")
    parser.add_argument("--plot-saved", action="store_true", help="Plot previously saved predictions")
    parser.add_argument("--maximum-evaluations", type=int)
    parser.add_argument("--relative-step", type=float, help="Relative finite-difference step in optimizer coordinates")
    parser.add_argument("--solver", choices=("BDF", "Radau"))
    parser.add_argument("--parameters", type=Path, help="Replay the five coefficients from a saved best-fit JSON")
    parser.add_argument("--resume", action="store_true", help="Start each fit at its saved best coefficients")
    parser.add_argument("--output", type=Path, help="Result directory; default results")
    parser.add_argument(
        "--workers", type=int, default=1, help="Independent fitting processes; user authorized up to 8 cores"
    )
    args = parser.parse_args()
    here = Path(__file__).resolve().parent
    root_directory = here.parents[2]
    with (here / "settings.toml").open("rb") as stream:
        config = tomllib.load(stream)
    if args.maximum_evaluations is not None:
        config["fit"]["maximum_evaluations"] = args.maximum_evaluations
    if args.relative_step is not None:
        config["fit"]["relative_step"] = args.relative_step
    if args.solver is not None:
        config["numerics"]["solver"] = args.solver
    if args.parameters is not None:
        config["fit"]["initial"] = json.loads(args.parameters.read_text(encoding="utf-8"))["parameters"]
    config["fit"]["resume"] = args.resume
    heads = load_population(here.parents[1] / "validation/head_comparison/records.toml", root=root_directory)
    if args.head:
        heads = tuple(head for head in heads if head.serial in args.head)
    destination = args.output or here / "results"
    destination.mkdir(exist_ok=True, parents=True)
    predictions = {}
    # Spawned fitting workers inherit one-thread numerical-library settings.
    for variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        os.environ[variable] = "1"
    if args.fit and args.workers > 1:
        with ProcessPoolExecutor(max_workers=min(args.workers, len(heads))) as pool:
            pending = [pool.submit(fit_worker, head, config, destination) for head in heads]
            for future in as_completed(pending):
                serial, prediction = future.result()
                predictions[serial] = prediction
    with warnings.catch_warnings():
        warnings.simplefilter("once", RuntimeWarning)
        candidate = build_candidate(config)
        stand = Stand.load(root_directory / "data/test_stand/test_stand.toml")
        for head in heads:
            if head.serial in predictions:
                continue
            if args.plot_saved:
                predictions[head.serial] = Prediction(
                    cooldown=xr.load_dataset(destination / f"{head.serial}_cooldown.nc"),
                    load_map=xr.load_dataset(destination / f"{head.serial}_map.nc"),
                )
            elif args.fit:
                predictions[head.serial] = fit_one(candidate, stand, head, config, destination=destination)
            else:
                prediction, diagnostics = predict_one(
                    candidate, stand, head, Parameters(*config["fit"]["initial"]), warm=initial_temperature(head)
                )
                predictions[head.serial] = prediction
                write_prediction(destination, head.serial, prediction)
                print(json.dumps({"head": head.serial, "diagnostics": diagnostics}, indent=2), flush=True)
    plt.switch_backend("Agg")
    figure = plot_comparison(heads, predictions=predictions, candidate_name="Five-coefficient candidate")
    stem = "comparison" if len(heads) == 6 else "comparison_" + "_".join(head.serial for head in heads)
    with comparison_style():
        for extension in ("png", "svg"):
            path = destination / f"{stem}.{extension}"
            figure.savefig(path, dpi=180)
            print(path, flush=True)
    plt.close(figure)


if __name__ == "__main__":
    main()
