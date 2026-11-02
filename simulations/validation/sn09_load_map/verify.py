"""Checks saved fit provenance, drive consistency, local identifiability, and gas quadrature."""

import argparse
import hashlib
import json
from dataclasses import replace
from unittest.mock import patch

import fit
import numpy as np
from scipy.special import roots_legendre

from pulse_tube_cooldown.library.links.cooler import head as head_equations


def perturbed(parameters, index, factor):
    """Multiplies one of a_1, a_2, C_1, C_2, eta by a positive factor."""
    values = [*parameters.packing_constant, *parameters.conductance, parameters.exchanger_effectiveness]
    if index == 2:
        values[index] += 1e-10 * np.log(factor)
    else:
        values[index] *= factor
    return replace(
        parameters,
        packing_constant=tuple(values[:2]),
        conductance=tuple(values[2:4]),
        exchanger_effectiveness=values[4],
    )


def main(folder=fit.HERE):
    """Writes local sensitivity and independent numerical consistency checks."""
    summary = json.loads((folder / "results.json").read_text())
    root = fit.HERE.resolve().parents[2]
    mismatches = [
        name
        for name, digest in summary["sha256"].items()
        if hashlib.sha256((root / name).read_bytes()).hexdigest() != digest
    ]
    if mismatches:
        raise RuntimeError(f"Inputs or code changed since the fit: {mismatches}")
    inputs, _, warm, drive = fit.read_experiment()
    head = replace(fit.read_head(inputs["head"]), cold_exchanger_coupling=summary["cold_exchanger_coupling"])
    parameters = fit.Parameters(
        packing_constant=tuple(summary["packing_constants"]),
        conductance=tuple(summary["conductances_W_per_Pa2"]),
        exchanger_effectiveness=summary["exchanger_effectiveness"],
        frequency=summary["frequency_Hz"],
    )
    templates = [np.array(row["temperature_K"]) for row in summary["points"]]
    indices = [head.stage_1, head.stage_2]
    matrix_indices = fit.column_indices(head)
    maximum_balance = maximum_pressure_change = 0.0
    for row, t in zip(summary["points"], templates, strict=True):
        mean, swing = fit.calc_pressures(head, parameters, drive, t)
        actual = fit.evaluate_head(head, parameters, drive, t, mean, swing)
        maximum_balance = max(maximum_balance, float(max(abs(actual.heat_flow[matrix_indices]))))
        maximum_pressure_change = max(
            maximum_pressure_change,
            abs(mean + swing / 2 - row["supply_pressure_Pa"]),
            abs(mean - swing / 2 - row["return_pressure_Pa"]),
        )
    if maximum_balance > 1e-6 or maximum_pressure_change > 0.01:
        raise RuntimeError("Archived fit failed independent drive and column balance checks.")

    def powers(p):
        return np.concatenate(
            [
                fit.solve_columns(
                    head, p, drive, warm_temperature=warm, stage_temperature=tuple(t[indices]), guess=t
                ).heat_flow[indices]
                for t in templates
            ]
        )

    order = ["log(a_1)", "log(a_2)", "C_1 / (1e-10 W/Pa^2)", "log(C_2)"]
    if head.cold_exchanger_coupling:
        order.append("log(eta)")
    jacobians = []
    for step in (1e-4, 5e-5):
        jacobians.append(
            np.column_stack(
                [
                    (powers(perturbed(parameters, j, np.exp(step))) - powers(perturbed(parameters, j, np.exp(-step))))
                    / (2 * step)
                    for j in range(len(order))
                ]
            )
        )
    singular = np.linalg.svd(jacobians[-1], compute_uv=False)

    def tube_density(warm_end, cold_end, pressure, order):
        nodes, weights = roots_legendre(order)
        temperature = np.exp(np.log(cold_end) + (nodes + 1) / 2 * np.log(warm_end / cold_end))
        return float(weights @ (fit.helium.density(temperature, pressure) * temperature) / (weights @ temperature))

    quadrature = []
    for row, t in zip(summary["points"], templates, strict=True):
        for cold in t[indices]:
            for pressure in (
                row["supply_pressure_Pa"],
                row["return_pressure_Pa"],
                (row["supply_pressure_Pa"] + row["return_pressure_Pa"]) / 2,
            ):
                quadrature.append(
                    abs(tube_density(t[0], cold, pressure, 48) / tube_density(t[0], cold, pressure, 96) - 1)
                )
    power_changes, temperature_changes = [], []
    with patch.object(head_equations, "_mean_tube_density", lambda w, c, p: tube_density(w, c, p, 96)):
        for row, t in zip(summary["points"], templates, strict=True):
            column = fit.solve_columns(
                head, parameters, drive, warm_temperature=warm, stage_temperature=tuple(t[indices]), guess=t
            )
            original = fit.evaluate_head(
                head,
                parameters,
                drive,
                t,
                (row["supply_pressure_Pa"] + row["return_pressure_Pa"]) / 2,
                row["supply_pressure_Pa"] - row["return_pressure_Pa"],
            )
            power_changes.append(float(max(abs(column.heat_flow[indices] - original.heat_flow[indices]))))
            temperature_changes.append(float(max(abs(column.temperature - t))))
    verification = dict(
        input_and_code_fingerprints_match=True,
        max_column_balance_W=maximum_balance,
        max_pressure_change_Pa=maximum_pressure_change,
        log_parameter_order=order,
        log_parameter_jacobian_W=jacobians[-1].tolist(),
        log_parameter_singular_values_W=singular.tolist(),
        log_parameter_condition_number=float(singular[0] / singular[-1]),
        jacobian_relative_change_1e4_to_5e5=float(
            np.linalg.norm(jacobians[-1] - jacobians[0]) / np.linalg.norm(jacobians[-1])
        ),
        max_tube_density_relative_change_48_to_96=max(quadrature),
        max_stage_power_change_W_48_to_96=max(power_changes),
        max_node_temperature_change_K_48_to_96=max(temperature_changes),
        interpretation="Local sensitivity on the selected branch; no statistical error model is assumed.",
    )
    (folder / "verification.json").write_text(json.dumps(verification, indent=2) + "\n")
    print(
        json.dumps({key: value for key, value in verification.items() if key != "log_parameter_jacobian_W"}, indent=2)
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folder", nargs="?", type=fit.Path, default=fit.HERE)
    main(parser.parse_args().folder)
