"""Screen fixed inlet power and fixed splitting with algebraic profiles."""

import json
import time
import tomllib
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import least_squares, root

from simulations.exploration.pulse_tube_candidate.model import Parameters, build_candidate
from simulations.exploration.pulse_tube_candidate.run import Stand
from simulations.validation.head_comparison.comparison import cell_derivatives, load_population


def main():
    """Fit cold heat balances and check the resulting powers along the data."""
    here = Path(__file__).resolve().parent
    repository = here.parents[2]
    audit = json.loads((here / "results/algebraic_transport_audit.json").read_text())
    config = tomllib.loads((here.parent / "pulse_tube_candidate/settings.toml").read_text())
    candidate = build_candidate(config)
    stand = Stand.load(repository / "data/test_stand/test_stand.toml")
    head = next(
        h
        for h in load_population(repository / "simulations/validation/head_comparison/records.toml", root=repository)
        if h.serial == "S0035"
    )
    warm = audit["rows"][0]["T_K"][0]
    temperatures = np.column_stack([head.load_map.data[k].values for k in ("T1", "T2")])
    heaters = np.column_stack([head.load_map.data[k].values for k in ("Q1", "Q2")])
    guesses = [np.array(row["centers_K"]) for row in audit["cold_checks"] if "map_index" in row]
    controls = audit["global_controls"]
    initial = np.array(
        [controls["constant_inlet_power_W"], controls["constant_branch_floor"], *audit["parameters"][3:]]
    )
    bounds = np.array(config["fit"]["lower"])[[0, 2, 3, 4]], np.array(config["fit"]["upper"])[[0, 2, 3, 4]]
    calls = 0
    best = None
    start = time.monotonic()

    def evaluate(centers, stages, values):
        inlet, ratio, exchange_1, exchange_2 = values
        boundary = np.r_[warm, stages]
        nominal = Parameters(inlet, 0, ratio, exchange_1, exchange_2)
        transmission = candidate.evaluate(centers, boundary, nominal).transmission[0]
        coefficients = Parameters(inlet * transmission / (1 + 1 / ratio), 0, ratio, exchange_1, exchange_2)
        return candidate.evaluate(centers, boundary, coefficients)

    def solve(stages, values, guess):
        def balance(coordinates):
            evaluation = evaluate(candidate._decode(coordinates), stages, values)
            return np.r_[evaluation.energy_rate, evaluation.flows[6] - evaluation.flows[7]]

        solution = root(balance, candidate._encode(guess), tol=1e-10)
        error = float(np.max(abs(balance(solution.x))))
        if error > config["numerics"]["equilibrium_flow_tolerance"]:
            solution = least_squares(balance, solution.x, max_nfev=100, gtol=1e-11, ftol=1e-11, xtol=1e-11)
            error = float(np.max(abs(balance(solution.x))))
        if error > config["numerics"]["equilibrium_flow_tolerance"]:
            raise RuntimeError(f"Profile residual {error:g} W at {stages.tolist()} K")
        centers = candidate._decode(solution.x)
        return centers, evaluate(centers, stages, values)

    scales = np.max(heaters, axis=0)

    def residual(coordinates):
        nonlocal calls, best
        calls += 1
        values = np.exp(coordinates)
        heat = []
        for i, stages in enumerate(temperatures):
            centers, evaluation = solve(stages, values, guesses[i])
            guesses[i] = centers
            heat.append(evaluation.heat_flow[[0, 2]] + stand.external_heat(stages, warm) + heaters[i])
        normalized = np.array(heat) / scales
        score = float(np.sum(normalized**2))
        if best is None or score < best["score"]:
            best = {"score": score, "parameters": values.tolist(), "heat_balance_residual_W": np.array(heat).tolist()}
        if calls % 10 == 0:
            print(json.dumps({"calls": calls, "seconds": time.monotonic() - start, "best": best}), flush=True)
        return normalized.ravel()

    report = {
        "description": (
            "Fixed inlet power, fixed branch ratio, two exchange conductances; algebraic profiles. "
            "Cold heat-balance fit only."
        ),
        "initial_parameters": initial.tolist(),
        "parameter_order": ["inlet_W", "stage_1_to_stage_2_branch_ratio", "exchange_1_W_per_K", "exchange_2_W_per_K"],
        "heat_residual_scales_W": scales.tolist(),
    }
    try:
        optimization = least_squares(
            residual,
            np.log(initial),
            bounds=np.log(bounds),
            max_nfev=config["fit"]["maximum_evaluations"],
            diff_step=1e-4,
            ftol=1e-8,
            xtol=1e-8,
            gtol=1e-8,
        )
        report["optimizer"] = {
            "success": bool(optimization.success),
            "message": optimization.message,
            "optimality": float(optimization.optimality),
            "evaluations": optimization.nfev,
        }
    except RuntimeError as error:
        report["fit_error"] = str(error)
    report["best"] = best
    if best is not None:
        values = np.array(best["parameters"])
        heat = np.array(best["heat_balance_residual_W"]) - heaters
        coordinates = list(map(tuple, heaters))

        def edges(field, low, high):
            corners = dict(zip(coordinates, field, strict=True))
            return np.column_stack(
                (
                    (corners[high, 0] - corners[low, 0] + corners[high, 0.5] - corners[low, 0.5]) / 2,
                    (corners[low, 0.5] - corners[low, 0] + corners[high, 0.5] - corners[high, 0]) / 2,
                )
            )

        report["secant_J_on_measured_cells_W_per_K"] = [
            (edges(heat, low, high) @ np.linalg.inv(edges(temperatures, low, high))).tolist()
            for low, high in ((0, 20), (20, 25))
        ]
        report["trajectory_power_checks"] = []
        for row in [*audit["rows"], *audit["cold_checks"]]:
            if "time_s" not in row or "error" in row:
                continue
            check = {"time_s": row["time_s"], "T_K": row["T_K"]}
            try:
                _, evaluation = solve(np.array(row["T_K"]), values, np.array(row["centers_K"]))
                check["Q_W"] = evaluation.heat_flow[[0, 2]].tolist()
            except RuntimeError as error:
                check["error"] = str(error)
            report["trajectory_power_checks"].append(check)
        predicted = head.load_map.data[["Q1", "Q2"]].copy()
        stage_temperatures = []
        balance_errors = []
        for i, (stages, heater) in enumerate(zip(temperatures, heaters, strict=True)):
            centers, _ = solve(stages, values, guesses[i])

            def map_balance(coordinates, heater=heater):
                temperature = candidate._decode(coordinates)
                evaluation = evaluate(temperature[:8], temperature[8:], values)
                return np.r_[
                    evaluation.energy_rate,
                    evaluation.flows[6] - evaluation.flows[7],
                    evaluation.heat_flow[[0, 2]] + stand.external_heat(temperature[8:], warm) + heater,
                ]

            solution = root(map_balance, candidate._encode(np.r_[centers, stages]), tol=1e-10)
            error = float(np.max(abs(map_balance(solution.x))))
            if error > config["numerics"]["equilibrium_flow_tolerance"]:
                raise RuntimeError(f"Load-map balance residual {error:g} W")
            stage_temperatures.append(candidate._decode(solution.x)[8:])
            balance_errors.append(error)
        for i, name in enumerate(("T1", "T2")):
            predicted[name] = ("point", np.array(stage_temperatures)[:, i])
        predicted.to_netcdf(here / "results/constant_input_four_coefficient_map.nc")
        report["predicted_map_T_K"] = np.array(stage_temperatures).tolist()
        report["map_rmse_K"] = [
            float(np.sqrt(np.mean((predicted[name] - head.load_map.data[name]) ** 2))) for name in ("T1", "T2")
        ]
        report["map_maximum_balance_residual_W"] = max(balance_errors)
        names = ("dQ1dT1", "dQ1dT2", "dQ2dT1", "dQ2dT2")
        observed_slopes, fitted_slopes = cell_derivatives(head.load_map.data), cell_derivatives(predicted)
        report["map_cell_J_W_per_K"] = {
            name: np.column_stack([slopes[key].values for key in names]).reshape(-1, 2, 2).tolist()
            for name, slopes in (("measured", observed_slopes), ("predicted", fitted_slopes))
        }
        plt.switch_backend("Agg")
        figure, axes = plt.subplots(2, 2, figsize=(10, 6.5), constrained_layout=True)
        for axis, name in zip(axes.flat, names, strict=True):
            axis.plot([0, 1], observed_slopes[name], "ko-", label="Measured")
            axis.plot([0, 1], fitted_slopes[name], "s--", color="#b54b31", label="Four-coefficient screen")
            axis.set_xticks([0, 1], ["0-20 W", "20-25 W"])
            axis.set(xlabel="Stage-1 heater interval", ylabel="W/K", title=name)
            axis.grid(alpha=0.2)
            axis.axhline(0, color="0.7", linewidth=0.7)
        axes[0, 0].legend(fontsize=9)
        figure.suptitle(
            "S0035: fixed inlet power and fixed split\nSlopes from solved load-map cells; cold-point fit only"
        )
        figure.savefig(here / "results/constant_input_four_coefficient_jacobians.png", dpi=170)
        plt.close(figure)
        report["warm_energy_audit"] = []
        report["rate_half_window_s"] = 120
        for row in report["trajectory_power_checks"]:
            if "Q_W" not in row or row["time_s"] == 0:
                continue
            center_time = row["time_s"]
            selected = head.cooldown.data.where(
                abs(head.cooldown.data.time - center_time) <= report["rate_half_window_s"], drop=True
            ).dropna("time")
            fits = np.array(
                [
                    np.polynomial.polynomial.polyfit(selected.time.values - center_time, selected[name].values, 2)
                    for name in ("T1", "T2")
                ]
            )
            stages, rates = fits[:, 0], fits[:, 1]
            capacity = np.array(
                [
                    node.calc_heat_capacity(float(temperature))
                    + sum(layer.mass * float(layer.material.specific_heat(temperature)) for layer in section.layers)
                    for node, section, temperature in zip(
                        stand.stages, candidate.design.regenerators, stages, strict=True
                    )
                ]
            )
            measured_power = stand.external_heat(stages, warm) - capacity * rates
            source_row = next(
                item
                for item in [*audit["rows"], *audit["cold_checks"]]
                if item.get("time_s") == center_time and "centers_K" in item
            )
            _, evaluation = solve(stages, values, np.array(source_row["centers_K"]))
            predicted_power = -evaluation.heat_flow[[0, 2]]
            report["warm_energy_audit"].append(
                {
                    "time_s": center_time,
                    "temperature_K": stages.tolist(),
                    "rate_K_per_s": rates.tolist(),
                    "capacity_J_per_K": capacity.tolist(),
                    "implied_head_cooling_W": measured_power.tolist(),
                    "predicted_head_cooling_W": predicted_power.tolist(),
                    "ratio": (predicted_power / measured_power).tolist(),
                }
            )
        report["warm_energy_assumptions"] = [
            "Each regenerator inventory is lumped at its corresponding stage temperature, "
            "as in the prior reduced screens.",
            "Configured stand radiation and no cooldown heater inputs.",
            "Cold map and cooldown use the same head parameters; the hose change remains a recorded run difference.",
        ]
    report["seconds"] = time.monotonic() - start
    destination = here / "results/constant_input_four_coefficient_screen.json"
    destination.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report), flush=True)


if __name__ == "__main__":
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        main()
