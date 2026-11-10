"""Replay and refine a saved candidate independently of acceptance tests."""

import argparse
import json
import time
import tomllib
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from hepak.properties import properties_TP
from scipy.integrate import quad

from simulations.exploration.pulse_tube_candidate.model import Parameters, build_candidate
from simulations.exploration.pulse_tube_candidate.run import Stand, initial_temperature, predict_one, simulate
from simulations.validation.head_comparison.comparison import load_population


def audit(config: dict) -> dict:
    """Reports property accuracy, energy inversion, and conservative balances."""
    candidate = build_candidate(config)
    sample = np.geomspace(2.11, 390, 137)
    exact = [properties_TP(float(t), config["source"]["mean_pressure"]) for t in sample]
    property_errors = {}
    for name, table, values in (
        ("density", candidate.helium.density, [p["density"] for p in exact]),
        ("cp", candidate.helium.cp, [p["cp"] for p in exact]),
        ("beta", candidate.helium.beta, [p["expansion_coefficient"] / t for t, p in zip(sample, exact, strict=True)]),
        ("viscosity", candidate.helium.viscosity, [p["viscosity"] for p in exact]),
        (
            "mass_response",
            candidate.helium.mass_response,
            [p["density"] * p["isothermal_compressibility"] for p in exact],
        ),
    ):
        values = np.asarray(values)
        property_errors[name] = float(np.max(np.abs(table.value(sample) / values - 1)))
    parameters = Parameters(*config["fit"]["initial"])
    boundary = np.array([295.0, 45.0, 3.0])
    start = time.perf_counter()
    centers, evaluation = candidate.equilibrate(boundary, parameters, compressors=1, frequency=1.0)
    equilibrium_seconds = time.perf_counter() - start
    energy = candidate.energies(centers, boundary)
    recovered, recovered_evaluation = candidate.recover_profile(
        energy, boundary, parameters, compressors=1, frequency=1.0
    )
    residual = (
        np.sum(evaluation.heat_flow)
        + np.sum(evaluation.energy_rate)
        - (evaluation.flows[0] - np.sum(evaluation.tube_power) + np.sum(evaluation.wall_leak))
    )
    off = candidate.evaluate(centers, boundary, parameters, compressors=0, frequency=1.0)
    off_expected = []
    temperature = candidate.complete_temperature(centers, boundary)
    for i, codes in enumerate(candidate.node_codes):
        section = candidate.design.regenerators[i]
        for k in range(len(codes) - 1):
            left, right = temperature[codes[k]], temperature[codes[k + 1]]
            resistance = (
                quad(
                    lambda x, section=section, left=left, right=right: (
                        1 / section.wall_material.thermal_conductivity(np.array([left + x * (right - left)]))[0]
                    ),
                    0,
                    1,
                    epsabs=1e-11,
                    epsrel=1e-10,
                )[0]
                * (candidate.positions[i][k + 1] - candidate.positions[i][k])
                / section.wall_area
            )
            off_expected.append((left - right) / resistance)
    # Adaptive quadrature independently evaluates the specified discrete closure.
    off_error = float(np.max(np.abs(off.flows - off_expected)))
    jacobian = candidate.energy_jacobian(centers, boundary)
    finite_difference = np.empty((7, 8))
    for index in range(8):
        perturbation = np.zeros(8)
        perturbation[index] = 1e-4
        finite_difference[:, index] = (
            candidate.energies(centers + perturbation, boundary) - candidate.energies(centers - perturbation, boundary)
        ) / 2e-4
    timings = []
    for _ in range(10):
        start = time.perf_counter()
        candidate.evaluate(centers, boundary, parameters, compressors=1, frequency=1.0)
        timings.append(time.perf_counter() - start)
    refined_config = {
        **config,
        "numerics": {
            **config["numerics"],
            "property_knots": 2 * config["numerics"]["property_knots"] - 1,
            "quadrature_order": 2 * config["numerics"]["quadrature_order"],
        },
    }
    refined = build_candidate(refined_config)
    refined_centers, refined_evaluation = refined.equilibrate(
        boundary, parameters, compressors=1, frequency=1.0, guess=centers
    )
    return {
        "helium_relative_errors_against_direct_hepak": property_errors,
        "parameters": parameters.as_array().tolist(),
        "boundary_K": boundary.tolist(),
        "centers_and_clamp_K": centers.tolist(),
        "energy_J": energy.tolist(),
        "flow_W": evaluation.flows.tolist(),
        "Q_W": evaluation.heat_flow.tolist(),
        "z": evaluation.pressure_fraction,
        "demand_kg_per_s": evaluation.demand,
        "s2": evaluation.profile_measure,
        "transmission": evaluation.transmission.tolist(),
        "equilibrium_residual_W": float(np.max(np.abs(evaluation.energy_rate))),
        "clamp_balance_W": float(evaluation.flows[6] - evaluation.flows[7]),
        "energy_balance_residual_W": float(residual),
        "energy_inversion_center_error_K": float(np.max(np.abs(recovered - centers))),
        "energy_inversion_flow_error_W": float(np.max(np.abs(recovered_evaluation.flows - evaluation.flows))),
        "drive_off_conduction_error_W": off_error,
        "energy_jacobian_maximum_error_J_per_K": float(np.max(np.abs(jacobian[:, :8] - finite_difference))),
        "refined_center_change_K": float(np.max(np.abs(refined_centers - centers))),
        "refined_Q_change_W": (refined_evaluation.heat_flow - evaluation.heat_flow).tolist(),
        "equilibrium_seconds": equilibrium_seconds,
        "median_evaluation_seconds": float(np.median(timings)),
    }


def reconstruct(candidate, head, parameters, cooldown, *, warm):
    """Recovers accepted profiles and source/heat budgets along a saved trace."""
    rows = []
    guess = np.full(8, warm)
    for index in range(cooldown.sizes["time"]):
        boundary = np.array([warm, cooldown.T1[index].item(), cooldown.T2[index].item()])
        guess, evaluated = candidate.recover_profile(
            cooldown.matrix_energy[index].values,
            boundary,
            parameters,
            compressors=head.cooldown.compressors,
            frequency=head.cooldown.valve_frequency,
            guess=guess,
        )
        rows.append(
            {
                "z": evaluated.pressure_fraction,
                "Dpass": evaluated.pass_demand,
                "Dfill": evaluated.fill_demand,
                "s2": evaluated.profile_measure,
                "W1": evaluated.tube_power[0],
                "W2c": evaluated.tube_power[1],
                "H1h": evaluated.flows[0],
                "H1c": evaluated.flows[2],
                "H2h": evaluated.flows[3],
                "H2c": evaluated.flows[9],
                "Q1": evaluated.heat_flow[0],
                "Q2": evaluated.heat_flow[2],
                "clamp": guess[7],
                "minimum_matrix_temperature": float(np.min(guess[:7])),
                "balance_residual": float(
                    np.sum(evaluated.heat_flow)
                    + np.sum(evaluated.energy_rate)
                    - evaluated.flows[0]
                    + np.sum(evaluated.tube_power)
                    - np.sum(evaluated.wall_leak)
                ),
            }
        )
    return xr.Dataset({key: ("time", [row[key] for row in rows]) for key in rows[0]}, coords={"time": cooldown.time})


def main():
    """Compares solvers or checks a fitted head under refined numerics."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--head")
    parser.add_argument("--audit", action="store_true", help="Check properties, conservation, and energy inversion")
    parser.add_argument("--evaluation", type=int, help="Replay a logged coefficient set, including rejected trials")
    parser.add_argument("--solver", choices=("BDF", "Radau"), default="Radau")
    parser.add_argument("--refine", action="store_true")
    parser.add_argument(
        "--cooldown-only", action="store_true", help="Verify transients independently of cold root selection"
    )
    args = parser.parse_args()
    here = Path(__file__).resolve().parent
    with (here / "settings.toml").open("rb") as stream:
        config = tomllib.load(stream)
    result_directory = here / "results"
    if args.audit:
        with warnings.catch_warnings():
            warnings.simplefilter("once", RuntimeWarning)
            result = audit(config)
        (result_directory / "numerical_audit.json").write_text(json.dumps(result, indent=2) + "\n")
        print(json.dumps(result, indent=2))
        return
    if args.head is None:
        parser.error("--head is required for a solver replay")
    head = next(
        head
        for head in load_population(here.parents[1] / "validation/head_comparison/records.toml", root=here.parents[2])
        if head.serial == args.head
    )
    if args.evaluation is None:
        record = json.loads((result_directory / f"{head.serial}_best.json").read_text(encoding="utf-8"))
    else:
        records = [
            json.loads(line)
            for line in (result_directory / f"{head.serial}_evaluations.jsonl").read_text().splitlines()
        ]
        record = next(row for row in records if row["evaluation"] == args.evaluation)
    parameters = Parameters(*record["parameters"])
    baseline = baseline_map = None
    if args.evaluation is None:
        baseline = xr.load_dataset(result_directory / f"{head.serial}_cooldown.nc")
        baseline_map = xr.load_dataset(result_directory / f"{head.serial}_map.nc")
    config["numerics"]["solver"] = args.solver
    if args.refine:
        config["numerics"]["property_knots"] = 2 * config["numerics"]["property_knots"] - 1
        config["numerics"]["quadrature_order"] *= 2
        config["numerics"]["relative_tolerance"] /= 4
        config["numerics"]["maximum_step"] /= 2
        config["numerics"]["sample_interval"] /= 2
    warm = initial_temperature(head)
    with warnings.catch_warnings():
        warnings.simplefilter("once", RuntimeWarning)
        candidate = build_candidate(config)
        stand = Stand.load(here.parents[2] / "data/test_stand/test_stand.toml")
        if args.cooldown_only:
            cooldown, transient_diagnostics = simulate(candidate, stand, head, parameters, warm=warm)
            prediction = None
            diagnostics = {"cooldown": transient_diagnostics}
        else:
            prediction, diagnostics = predict_one(candidate, stand, head, parameters, warm=warm)
            cooldown = prediction.cooldown
        budgets = reconstruct(candidate, head, parameters, cooldown, warm=warm)
    tag = f"{head.serial}_{args.solver}" + ("_refined" if args.refine else "")
    if args.cooldown_only:
        tag += "_cooldown_only"
    if args.evaluation is not None:
        tag += f"_evaluation_{args.evaluation}"
    output = result_directory / "verification"
    output.mkdir(exist_ok=True)
    cooldown.to_netcdf(output / f"{tag}_cooldown.nc")
    if prediction is not None:
        prediction.load_map.to_netcdf(output / f"{tag}_map.nc")
    budgets.to_netcdf(output / f"{tag}_budgets.nc")
    report = {
        "parameters": parameters.as_array().tolist(),
        "diagnostics": diagnostics,
        "numerics": config["numerics"],
        "source_at_start": {key: budgets[key][0].item() for key in ("z", "Dpass", "Dfill")},
        "source_at_end": {key: budgets[key][-1].item() for key in ("z", "Dpass", "Dfill")},
        "maximum_balance_residual_W": float(np.max(np.abs(budgets.balance_residual))),
    }
    if args.evaluation is None:
        report["maximum_temperature_change_K"] = {
            name: float(np.max(np.abs(cooldown[name].interp(time=baseline.time) - baseline[name])))
            for name in ("T1", "T2")
        }
        if prediction is not None:
            report["maximum_map_change_K"] = {
                name: float(np.max(np.abs(prediction.load_map[name] - baseline_map[name]))) for name in ("T1", "T2")
            }
    (output / f"{tag}.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2), flush=True)
    plt.switch_backend("Agg")
    figure, axes = plt.subplots(3, 1, figsize=(9, 8), sharex=True)
    axes[0].plot(budgets.time / 3600, budgets.z)
    axes[0].set_ylabel("Pressure fraction z")
    axes[1].plot(budgets.time / 3600, 1e3 * budgets.Dpass, label="Through-flow demand")
    axes[1].plot(budgets.time / 3600, 1e3 * budgets.Dfill, label="Filling demand")
    axes[1].axhline(
        1e3 * head.cooldown.compressors * config["source"]["delivery_per_compressor"],
        color="k",
        ls=":",
        label="Available compressor delivery",
    )
    axes[1].set_ylabel("Ceiling demand (g/s)")
    axes[1].legend()
    for name in ("W1", "W2c"):
        axes[2].plot(budgets.time / 3600, budgets[name], label=name)
    axes[2].set(xlabel="Time since compressor start (h)", ylabel="Delivered power (W)")
    axes[2].legend()
    figure.suptitle(head.serial + " source diagnostics")
    figure.tight_layout()
    figure.savefig(output / f"{tag}_source.png", dpi=150)
    plt.close(figure)


if __name__ == "__main__":
    main()
