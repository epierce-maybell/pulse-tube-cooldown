"""Fit a two-temperature stand model with bounded, coupled cooling powers.

For each stage, u_i = g_i*(T_i - B_i)/A_i and
P_i = A_i*tanh(u_i - r_i*tanh(u_other)). B is in K, g in W/K,
A in W, and r is dimensionless. P is positive cooling power.
Only the two external stage temperatures are integrated. Heat capacity and
radiation use the existing stand template. The fitted net head powers include
head losses. The unloaded clamp is implicit in this stand-only experiment.
"""

import argparse
import hashlib
import json
import time
import tomllib
import warnings
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Protocol

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from numpy.typing import NDArray
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator
from scipy.optimize import least_squares, root

from pulse_tube_cooldown.library.heads.load import load_head_model
from simulations.exploration.pulse_tube_candidate.run import Stand, initial_temperature
from simulations.exploration.pulse_tube_reduced.five_coefficient import FiveCoefficientLaw, calibrate_map
from simulations.validation.head_comparison.comparison import (
    HeadData,
    Prediction,
    cell_derivatives,
    load_population,
)


class PowerLaw(Protocol):
    """Net head cooling powers evaluated at the two stage temperatures."""

    def power(self, temperature: NDArray) -> NDArray:
        """Returns cooling powers in W, with stages on the last axis."""
        ...


@dataclass(frozen=True, slots=True)
class CoolingLaw:
    """Two offsets, slopes, limits, and cross-stage strengths."""

    base: NDArray
    slope: NDArray
    limit: NDArray
    coupling: NDArray

    @classmethod
    def from_array(cls, values: NDArray) -> "CoolingLaw":
        """Reads coefficients in settings.toml order."""
        return cls(*np.asarray(values).reshape(4, 2))

    def power(self, temperature: NDArray) -> NDArray:
        """Returns net head cooling power in W, with stages on the last axis."""
        scaled = self.slope * (temperature - self.base) / self.limit
        return self.limit * np.tanh(scaled - self.coupling * np.tanh(scaled[..., ::-1]))


@dataclass(frozen=True, slots=True)
class Experiment:
    """Stand properties and numerical settings for one warm boundary."""

    capacity: PchipInterpolator
    radiation: PchipInterpolator
    warm: float
    numerics: dict

    @classmethod
    def build(cls, stand: Stand, warm: float, numerics: dict) -> "Experiment":
        """Tabulates the configured stand properties once for the fitting loop."""
        grid = np.geomspace(
            numerics["property_min_temperature"], numerics["property_max_temperature"], numerics["property_knots"]
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            capacity = np.column_stack(
                [
                    sum(mass * material.specific_heat(grid) for material, mass in node.composition)
                    for node in stand.stages
                ]
            )
            radiation = np.array([stand.external_heat(np.array([t, t]), warm)[0] for t in grid])
        return cls(PchipInterpolator(grid, capacity), PchipInterpolator(grid, radiation), warm, numerics)

    def cooldown(self, law: PowerLaw, times: NDArray) -> xr.Dataset:
        """Integrates C_i(T_i)*dT_i/dt = radiation_i - P_i."""

        def rhs(_time: float, temperature: NDArray) -> NDArray:
            capacities = np.diag(self.capacity(temperature))
            return (np.array([self.radiation(temperature[0]), 0.0]) - law.power(temperature)) / capacities

        solution = solve_ivp(
            rhs,
            (0.0, float(times[-1])),
            np.full(2, self.warm),
            method="Radau",
            t_eval=times,
            rtol=self.numerics["relative_tolerance"],
            atol=self.numerics["absolute_tolerance"],
            max_step=self.numerics["maximum_step"],
        )
        if not solution.success:
            raise RuntimeError(solution.message)
        return xr.Dataset(
            {name: ("time", solution.y[i], {"units": "K"}) for i, name in enumerate(("T1", "T2"))},
            coords={"time": ("time", solution.t, {"units": "s"})},
            attrs={"complete": 1, "state_count": 2},
        )

    def load_map(self, law: PowerLaw, measured: xr.Dataset) -> xr.Dataset:
        """Balances each measured heater pair using the same cooling law."""
        stages = []
        for index in range(measured.sizes["point"]):
            heater = np.array([measured[name][index].item() for name in ("Q1", "Q2")])

            def residual(temperature: NDArray, heater: NDArray = heater) -> NDArray:
                return law.power(temperature) - heater - np.array([self.radiation(temperature[0]), 0.0])

            guess = np.array([measured[name][index].item() for name in ("T1", "T2")])
            solution = root(residual, guess, tol=1e-10)
            error = np.max(np.abs(residual(solution.x)))
            if error > 1e-7:
                raise RuntimeError(f"Load-map balance error: {error} W.")
            stages.append(solution.x)
        result = measured[["Q1", "Q2"]].copy()
        for i, name in enumerate(("T1", "T2")):
            result[name] = ("point", np.asarray(stages)[:, i], {"units": "K"})
        return result


def fitting_samples(head: HeadData, interval: float) -> tuple[NDArray, NDArray]:
    """Selects observations at uniform time spacing, retaining logger gaps."""
    times = np.r_[np.arange(0, float(head.cooldown.data.time[-1]), interval), float(head.cooldown.data.time[-1])]
    observed = head.cooldown.data.interp(time=times)
    temperatures = np.column_stack([observed[name].values for name in ("T1", "T2")])
    temperatures[temperatures[:, 0] >= head.stage_1_max_temperature, 0] = np.nan
    return times, temperatures


def draw(
    head: HeadData,
    prediction: Prediction,
    destination: Path,
    *,
    inventory: str,
    candidate_name: str = "Two temperatures",
) -> None:
    """Writes readable S0035 cooldown, map, and response comparisons."""
    figure, axes = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
    for axis, name in zip(axes, ("T1", "T2"), strict=True):
        observed = head.cooldown.data[name]
        if name == "T1":
            observed = observed.where(observed < head.stage_1_max_temperature)
        axis.plot(observed.time / 3600, observed, color="black", label="Measured")
        axis.plot(prediction.cooldown.time / 3600, prediction.cooldown[name], label=candidate_name)
        baseline = (
            Path(__file__).resolve().parents[1] / "pulse_tube_candidate/results/prescribed_pressure/S0035_cooldown.nc"
        )
        if baseline.exists():
            with xr.open_dataset(baseline) as data:
                axis.plot(data.time / 3600, data[name], color="0.65", linestyle="--", label="Internal states (partial)")
        axis.set(xlabel="Time (h)", ylabel=f"{name} (K)", title=f"{head.serial}: {name}")
        axis.grid(alpha=0.2)
    axes[0].set_ylim(30, 220)
    axes[1].set_yscale("log")
    axes[1].legend()
    figure.suptitle(f"Heat capacity: {inventory}")
    figure.savefig(destination / "cooldown.png", dpi=180)
    plt.close(figure)

    figure, axes = plt.subplots(2, 3, figsize=(13, 7), constrained_layout=True)
    axis = axes[0, 0]
    for data, style, label in (
        (head.load_map.data, "ko-", "Measured"),
        (prediction.load_map, "C0s--", candidate_name),
    ):
        for i, power in enumerate(np.unique(data.Q2)):
            row = data.where(data.Q2 == power, drop=True).sortby("Q1")
            axis.plot(row.T1, row.T2, style, label=label if i == 0 else None)
    axis.set(xlabel="T1 (K)", ylabel="T2 (K)", title="Cold load map")
    axis.legend()
    measured = cell_derivatives(head.load_map.data)
    predicted = cell_derivatives(prediction.load_map)
    for axis, name, coordinate in zip(
        [axes[0, 1], axes[0, 2], axes[1, 1], axes[1, 2]],
        ("dQ1dT1", "dQ1dT2", "dQ2dT1", "dQ2dT2"),
        ("T1", "T2", "T1", "T2"),
        strict=True,
    ):
        axis.plot(measured[coordinate], measured[name], "ko", label="Measured")
        axis.plot(predicted[coordinate], predicted[name], "C0s", label=candidate_name)
        axis.set(xlabel=f"Cell mean {coordinate} (K)", ylabel="W/K", title=name)
        axis.grid(alpha=0.2)
    axes[1, 0].axis("off")
    axes[1, 0].text(
        0,
        0.85,
        "Heat-positive derivatives\n\nCross-couplings evaluated from map cells.\n"
        "They are not independent fit targets.\n\nSame parameters for map and cooldown.\n"
        "Heat capacity: " + inventory.replace(" + ", "\n+ ") + ".\nNo independent regenerator states.",
        va="top",
    )
    figure.savefig(destination / "map_response.png", dpi=180)
    plt.close(figure)


def main() -> None:
    """Runs a serial fit and independently tightens integration tolerances."""
    plt.switch_backend("Agg")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lump-regenerators", action="store_true")
    parser.add_argument("--law", choices=("bounded", "five"), default="bounded")
    args = parser.parse_args()
    here = Path(__file__).resolve().parent
    repository = here.parents[2]
    config = tomllib.loads((here / "settings.toml").read_text())
    law_type = CoolingLaw
    candidate_name = "Eight coefficients"
    source_paths = [Path(__file__), here / "settings.toml"]
    if args.law == "five":
        config.update(tomllib.loads((here / "five_coefficient.toml").read_text()))
        law_type = FiveCoefficientLaw
        candidate_name = "Five coefficients"
        source_paths += [here / "five_coefficient.py", here / "five_coefficient.toml"]
    head = next(
        h for h in load_population(repository / config["records"], root=repository) if h.serial == config["serial"]
    )
    stand = Stand.load(repository / config["stand"])
    inventory = "stand"
    if args.lump_regenerators:
        design = load_head_model(config["head_design"])
        stages = tuple(
            replace(
                node, composition=node.composition + tuple((layer.material, layer.mass) for layer in section.layers)
            )
            for node, section in zip(stand.stages, design.regenerators, strict=True)
        )
        stand = replace(stand, stages=stages)
        inventory = "stand + stage-lumped regenerators"
    experiment = Experiment.build(stand, initial_temperature(head), config["numerics"])
    times, observed = fitting_samples(head, config["fit"]["sample_interval"])
    destination = here / "results"
    if args.law == "five":
        destination /= "five_coefficient"
    if args.lump_regenerators:
        destination /= "lumped_regenerators"
    destination.mkdir(exist_ok=True, parents=True)
    start = np.array(config["law"]["initial"])
    if args.law == "five":
        temperatures = np.column_stack([head.load_map.data[name] for name in ("T1", "T2")])
        powers = np.column_stack([head.load_map.data[name] for name in ("Q1", "Q2")])
        powers[:, 0] += experiment.radiation(temperatures[:, 0])
        start = calibrate_map(temperatures, powers)
        map_calibration = law_type.from_array(start)
        calibration_prediction = Prediction(
            cooldown=experiment.cooldown(map_calibration, times),
            load_map=experiment.load_map(map_calibration, head.load_map.data),
        )
        calibration_directory = destination / "map_calibration"
        calibration_directory.mkdir(exist_ok=True)
        calibration_prediction.cooldown.to_netcdf(calibration_directory / "cooldown.nc")
        calibration_prediction.load_map.to_netcdf(calibration_directory / "load_map.nc")
        (calibration_directory / "parameters.json").write_text(json.dumps(start.tolist()) + "\n")
        draw(
            head,
            calibration_prediction,
            calibration_directory,
            inventory=inventory,
            candidate_name="Five coefficients: map calibration",
        )
    calls = 0
    started = time.perf_counter()

    def objective(values: NDArray) -> NDArray:
        nonlocal calls
        law = law_type.from_array(values)
        cooldown = experiment.cooldown(law, times)
        load_map = experiment.load_map(law, head.load_map.data)
        blocks = []
        for i, name in enumerate(("T1", "T2")):
            valid = np.isfinite(observed[:, i])
            blocks.append(
                np.log(cooldown[name].values[valid] / observed[valid, i])
                / (config["fit"]["cooldown_log_scales"][i] * np.sqrt(valid.sum()))
            )
            error = (load_map[name].values - head.load_map.data[name].values) / config["fit"]["map_temperature_scales"][
                i
            ]
            blocks.append(error / np.sqrt(error.size))
        calls += 1
        residual = np.concatenate(blocks)
        if calls % 100 == 0:
            print(
                f"{calls} calls; squared scaled residual {residual @ residual:.6g}; "
                f"{time.perf_counter() - started:.1f} s",
                flush=True,
            )
        return residual

    def jacobian(values: NDArray) -> NDArray:
        steps = config["fit"]["difference_fraction"] * np.asarray(config["law"]["initial"])
        directions = np.diag(steps)
        return np.column_stack(
            [
                (objective(values + direction) - objective(values - direction)) / (2 * step)
                for direction, step in zip(directions, steps, strict=True)
            ]
        )

    fit = least_squares(
        objective,
        start,
        bounds=(config["law"]["lower"], config["law"]["upper"]),
        x_scale="jac",
        max_nfev=config["fit"]["maximum_evaluations"],
        jac=jacobian,
        ftol=1e-7,
        xtol=1e-7,
        gtol=1e-7,
    )
    law = law_type.from_array(fit.x)
    prediction = Prediction(
        cooldown=experiment.cooldown(law, times), load_map=experiment.load_map(law, head.load_map.data)
    )
    refined = replace(
        experiment,
        numerics={
            **experiment.numerics,
            "relative_tolerance": 1e-10,
            "absolute_tolerance": 1e-11,
            "maximum_step": 15.0,
        },
    )
    verification = refined.cooldown(law, times)
    prediction.cooldown.to_netcdf(destination / "cooldown.nc")
    prediction.load_map.to_netcdf(destination / "load_map.nc")
    report = {
        "serial": head.serial,
        "inventory": inventory,
        "law": args.law,
        "fitted_coefficient_count": len(fit.x),
        "starting_parameters": start.tolist(),
        "config": config,
        "parameters": fit.x.tolist(),
        "optimizer": {
            "success": bool(fit.success),
            "message": fit.message,
            "optimality": float(fit.optimality),
            "evaluations": fit.nfev,
            "calls": calls,
            "active_bounds": fit.active_mask.tolist(),
            "scaled_squared_residual": float(2 * fit.cost),
        },
        "elapsed_seconds": time.perf_counter() - started,
        "cooldown_rmse_K": {
            name: float(np.sqrt(np.nanmean((prediction.cooldown[name].values - observed[:, i]) ** 2)))
            for i, name in enumerate(("T1", "T2"))
        },
        "map_rmse_K": {
            name: float(np.sqrt(np.mean((prediction.load_map[name].values - head.load_map.data[name].values) ** 2)))
            for name in ("T1", "T2")
        },
        "refined_maximum_difference_K": {
            name: float(np.max(np.abs(verification[name] - prediction.cooldown[name]))) for name in ("T1", "T2")
        },
        "map_cells_measured": cell_derivatives(head.load_map.data).to_dict(),
        "map_cells_predicted": cell_derivatives(prediction.load_map).to_dict(),
        "input_sha256": {
            str(path.relative_to(repository)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in (
                *source_paths,
                repository / config["stand"],
                repository / config["records"],
                head.cooldown.source,
                head.load_map.source,
                repository / "src/pulse_tube_cooldown/library/materials/load.py",
                repository / "src/pulse_tube_cooldown/library/heads/data/pt1_5.toml",
            )
        },
        "scope": (
            "Stand-only feasibility screen; clamp and drive dependence unresolved. "
            "Cooldown follows a hose change relative to the load map."
        ),
    }
    (destination / "fit.json").write_text(json.dumps(report, indent=2) + "\n")
    draw(head, prediction, destination, inventory=inventory, candidate_name=candidate_name)
    print(
        json.dumps(
            {k: report[k] for k in ("optimizer", "cooldown_rmse_K", "map_rmse_K", "refined_maximum_difference_K")},
            indent=2,
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
