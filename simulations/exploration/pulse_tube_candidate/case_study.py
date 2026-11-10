"""Audit S0035 equilibrium branches and nearby failed parameter evaluations."""

import json
import tomllib
import warnings
from pathlib import Path

import numpy as np
import xarray as xr
from scipy.integrate import solve_ivp
from scipy.linalg import eigvals
from scipy.optimize import root

from simulations.exploration.pulse_tube_candidate.model import Parameters, build_candidate
from simulations.exploration.pulse_tube_candidate.run import Stand
from simulations.validation.head_comparison.comparison import cell_derivatives, load_population


def main() -> None:
    """Continue the audited roots and check their full thermal stability."""
    here = Path(__file__).resolve().parent
    results = here / "results"
    output = results / "diagnostics"
    output.mkdir(exist_ok=True)
    config = tomllib.loads((here / "settings.toml").read_text())
    study = tomllib.loads((here / "case_study.toml").read_text())
    audit = json.loads((results / "S0035_branch_audit_refined.json").read_text())
    best = json.loads((results / "S0035_best.json").read_text())
    parameters = Parameters(*best["parameters"])
    warm = best["warm_temperature_K"]
    head = next(
        head
        for head in load_population(here.parents[1] / "validation/head_comparison/records.toml", root=here.parents[2])
        if head.serial == "S0035"
    )
    candidate = build_candidate(config)
    stand = Stand.load(here.parents[2] / "data/test_stand/test_stand.toml")

    def flows(temperature, heaters, coefficients):
        centers, stage = temperature[:8], temperature[8:]
        evaluated = candidate.evaluate(
            centers,
            np.r_[warm, stage],
            coefficients,
        )
        return np.r_[
            evaluated.energy_rate,
            evaluated.flows[6] - evaluated.flows[7],
            evaluated.heat_flow[[0, 2]] + stand.external_heat(stage, warm) + heaters,
        ]

    def solve(guess, heaters, coefficients):
        solution = root(
            lambda coordinate: flows(candidate._decode(coordinate), heaters, coefficients),
            candidate._encode(guess),
            tol=1e-10,
        )
        temperature = candidate._decode(solution.x)
        error = float(np.max(np.abs(flows(temperature, heaters, coefficients))))
        if error > config["numerics"]["equilibrium_flow_tolerance"]:
            raise RuntimeError(f"Continuation residual {error:g} W at {heaters} W")
        return temperature, error

    def stability(temperature, heaters, coefficients):
        mass = np.zeros((10, 10))
        energy_derivative = candidate.energy_jacobian(temperature[:8], np.r_[warm, temperature[8:]])
        mass[:7] = energy_derivative[:, [*range(8), 9, 10]]
        for index, node in enumerate(stand.stages):
            mass[8 + index, 8 + index] = node.calc_heat_capacity(float(temperature[8 + index]))
        derivative = np.empty((10, 10))
        for column in range(10):
            step = np.zeros(10)
            step[column] = 1e-4
            derivative[:, column] = (
                flows(temperature + step, heaters, coefficients) - flows(temperature - step, heaters, coefficients)
            ) / (2 * step[column])
        rates = eigvals(derivative, mass)
        finite = rates[np.isfinite(rates)]
        return {"finite_mode_count": len(finite), "maximum_growth_rate_per_s": float(finite.real.max())}

    reports = []
    stable_seeds = [
        entry
        for entry in audit["roots"]
        if "full_growth_rates_real_per_s" in entry and max(entry["full_growth_rates_real_per_s"]) < 0
    ]
    stable_seeds = [
        min(stable_seeds, key=lambda entry: entry["stages_K"][1]),
        max(stable_seeds, key=lambda entry: entry["stages_K"][1]),
    ]
    for seed in stable_seeds:
        initial, error = solve(np.r_[seed["centers_K"], seed["stages_K"]], np.array(audit["heaters_W"]), parameters)
        branch = {"anchor_stages_K": initial[8:].tolist(), "points": []}
        for index in range(head.load_map.data.sizes["point"]):
            heaters = np.array([head.load_map.data[name][index].item() for name in ("Q1", "Q2")])
            guess = initial.copy()
            point = {"heaters_W": heaters.tolist()}
            try:
                for fraction in np.linspace(0, 1, 21)[1:]:
                    intermediate = np.array(audit["heaters_W"]) * (1 - fraction) + heaters * fraction
                    guess, error = solve(guess, intermediate, parameters)
                point.update(
                    {
                        "temperature_K": guess.tolist(),
                        "balance_residual_W": error,
                        **stability(guess, heaters, parameters),
                    }
                )
            except RuntimeError as error:
                point["error"] = str(error)
            branch["points"].append(point)
        reports.append(branch)

    rows = [json.loads(line) for line in (results / "S0035_evaluations.jsonl").read_text().splitlines()]
    trials = []
    for row in rows:
        if row["evaluation"] not in study["failed_evaluations"]:
            continue
        trial = {"evaluation": row["evaluation"], "original_error": row.get("error"), "roots": []}
        coefficients = Parameters(*row["parameters"])
        for branch in reports:
            seed = branch["points"][audit["point"]]
            if "temperature_K" not in seed:
                continue
            try:
                temperature, error = solve(np.array(seed["temperature_K"]), np.array(audit["heaters_W"]), coefficients)
                trial["roots"].append(
                    {
                        "stages_K": temperature[8:].tolist(),
                        "balance_residual_W": error,
                        **stability(temperature, np.array(audit["heaters_W"]), coefficients),
                    }
                )
            except RuntimeError as error:
                trial["roots"].append({"error": str(error)})
        trials.append(trial)
    initial_point = reports[0]["points"][study["load_step_initial_point"]]
    initial = np.array(initial_point["temperature_K"])
    heaters = np.array(audit["heaters_W"])
    profile = initial[:8].copy()
    state = np.r_[candidate.energies(profile, np.r_[warm, initial[8:]]), np.log(initial[8:] / warm)]

    def rhs(_time, value):
        nonlocal profile
        stages = warm * np.exp(value[7:])
        profile, evaluated = candidate.recover_profile(
            value[:7],
            np.r_[warm, stages],
            parameters,
            guess=profile,
        )
        capacities = np.array([node.calc_heat_capacity(float(t)) for node, t in zip(stand.stages, stages, strict=True)])
        heat = evaluated.heat_flow[[0, 2]] + stand.external_heat(stages, warm) + heaters
        return np.r_[evaluated.energy_rate, heat / (capacities * stages)]

    solution = solve_ivp(
        rhs,
        (0, study["load_step_duration"]),
        state,
        method="Radau",
        rtol=config["numerics"]["relative_tolerance"],
        atol=np.r_[
            np.full(7, config["numerics"]["energy_absolute_tolerance"]),
            np.full(2, config["numerics"]["temperature_absolute_tolerance"] / warm),
        ],
        max_step=config["numerics"]["maximum_step"],
    )
    endpoint = np.r_[profile, warm * np.exp(solution.y[7:, -1])]
    trajectory = {
        "success": bool(solution.success),
        "message": solution.message,
        "duration_s": float(solution.t[-1]),
        "initial_heaters_W": initial_point["heaters_W"],
        "final_heaters_W": heaters.tolist(),
        "final_stages_K": endpoint[8:].tolist(),
    }
    if solution.success:
        endpoint, error = solve(endpoint, heaters, parameters)
        trajectory.update({"relaxed_stages_K": endpoint[8:].tolist(), "balance_residual_W": error})
    xr.Dataset(
        {"T1": ("time", warm * np.exp(solution.y[7])), "T2": ("time", warm * np.exp(solution.y[8]))},
        coords={"time": solution.t},
    ).to_netcdf(output / "S0035_load_step.nc")
    branch_map = head.load_map.data.copy()
    for column, name in enumerate(("T1", "T2")):
        branch_map[name] = ("point", [point["temperature_K"][8 + column] for point in reports[0]["points"]])
    branch_map.to_netcdf(output / "S0035_lower_branch_map.nc")
    responses = cell_derivatives(branch_map)
    report = {
        "parameters": best["parameters"],
        "study": study,
        "numerics": config["numerics"],
        "branches": reports,
        "previously_failed_trials": trials,
        "load_step": trajectory,
        "lower_branch_responses": {
            name: responses[name].values.tolist() for name in ("dQ1dT1", "dQ1dT2", "dQ2dT1", "dQ2dT2")
        },
    }
    (output / "S0035_continuation.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    with warnings.catch_warnings():
        warnings.simplefilter("once", RuntimeWarning)
        main()
