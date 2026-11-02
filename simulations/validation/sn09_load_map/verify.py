"""Independent numerical checks for the archived SN09 branch fit.

Run after fit.py. This checks pressure recomputation, quadrature convergence,
local parameter sensitivity, and the recorded input/code fingerprints.
"""

import hashlib
import json
from pathlib import Path
from unittest.mock import patch

import fit
import numpy as np
from scipy.special import roots_legendre

from pulse_tube_cooldown.library.links.cooler import head as head_equations


def main() -> None:
    """Writes numerical convergence and local identifiability diagnostics."""
    folder = Path(__file__).parent
    summary = json.loads((folder / "results.json").read_text(encoding="utf-8"))
    root = folder.resolve().parents[2]
    mismatches = [
        name
        for name, digest in summary["sha256"].items()
        if hashlib.sha256((root / name).read_bytes()).hexdigest() != digest
    ]
    if mismatches:
        raise RuntimeError(f"Inputs or code changed since the fit: {mismatches}")
    _, load_map, warm, drive = fit.read_experiment()
    head = fit.read_head("PT1.1")
    parameters = fit.Parameters(
        packing_constant=tuple(summary["packing_constants"]),
        conductance=tuple(summary["conductances_W_per_Pa2"]),
        frequency=summary["frequency_Hz"],
    )
    templates = [np.array(row["temperature_K"]) for row in summary["points"]]
    stage_indices = [head.stage_1, head.stage_2]
    basis_a1 = []
    basis_a2 = []
    basis_c1 = []
    normalized = [
        fit.Parameters(
            packing_constant=(parameters.packing_constant[0], 0.0),
            conductance=(0.0, 0.0),
            frequency=parameters.frequency,
        ),
        fit.Parameters(
            packing_constant=(0.0, parameters.packing_constant[1]),
            conductance=(0.0, parameters.conductance[1]),
            frequency=parameters.frequency,
        ),
        fit.Parameters(
            packing_constant=(0.0, 0.0), conductance=(parameters.conductance[0], 0.0), frequency=parameters.frequency
        ),
    ]
    for t in templates:
        mean, swing = fit.calc_pressures(head, parameters, drive, t)
        values = [fit.evaluate_head(head, p, drive, t, mean, swing).heat_flow[stage_indices] for p in normalized]
        basis_a1.extend(values[0])
        basis_a2.extend(values[1])
        basis_c1.extend(values[2])
    jacobians = []
    for step in [1e-4, 1e-5]:
        responses = []
        for sign in [1, -1]:
            changed = fit.Parameters(
                packing_constant=parameters.packing_constant,
                conductance=(parameters.conductance[0], parameters.conductance[1] * np.exp(sign * step)),
                frequency=parameters.frequency,
            )
            responses.append(
                np.concatenate(
                    [
                        fit.solve_columns(
                            head,
                            changed,
                            drive,
                            warm_temperature=warm,
                            stage_temperature=(t[head.stage_1], t[head.stage_2]),
                            guess=t,
                        ).heat_flow[stage_indices]
                        for t in templates
                    ]
                )
            )
        c2_derivative = (responses[0] - responses[1]) / (2 * step)
        jacobians.append(np.column_stack([basis_a1, np.array(basis_a2) - c2_derivative, basis_c1, c2_derivative]))
    singular_values = np.linalg.svd(jacobians[-1], compute_uv=False)

    def tube_density(warm_end: float, cold_end: float, pressure: float, order: int) -> float:
        nodes, weights = roots_legendre(order)
        temperature = np.exp(np.log(cold_end) + (nodes + 1) / 2 * np.log(warm_end / cold_end))
        return float(weights @ (fit.helium.density(temperature, pressure) * temperature) / (weights @ temperature))

    quadrature = []
    for row, t in zip(summary["points"], templates, strict=True):
        for cold in t[stage_indices]:
            for pressure in [
                row["supply_pressure_Pa"],
                row["return_pressure_Pa"],
                (row["supply_pressure_Pa"] + row["return_pressure_Pa"]) / 2,
            ]:
                density_48 = tube_density(t[0], cold, pressure, 48)
                density_96 = tube_density(t[0], cold, pressure, 96)
                quadrature.append(abs(density_48 / density_96 - 1))

    def density_96(warm_end: float, cold_end: float, pressure: float) -> float:
        return tube_density(warm_end, cold_end, pressure, 96)

    power_changes = []
    temperature_changes = []
    with patch.object(head_equations, "_mean_tube_density", density_96):
        for t in templates:
            column = fit.solve_columns(
                head,
                parameters,
                drive,
                warm_temperature=warm,
                stage_temperature=(t[head.stage_1], t[head.stage_2]),
                guess=t,
            )
            power_changes.append(column)
            temperature_changes.append(float(max(abs(column.temperature - t))))
    stage_power_changes = []
    for row, t, column in zip(summary["points"], templates, power_changes, strict=True):
        original = fit.evaluate_head(
            head,
            parameters,
            drive,
            t,
            (row["supply_pressure_Pa"] + row["return_pressure_Pa"]) / 2,
            row["supply_pressure_Pa"] - row["return_pressure_Pa"],
        )
        stage_power_changes.append(float(max(abs(column.heat_flow[stage_indices] - original.heat_flow[stage_indices]))))

    profile_rows = []
    heater = np.column_stack([load_map.Q1.values, load_map.Q2.values])
    for record in [
        *json.loads((folder / "scan_cold.json").read_text()),
        *json.loads((folder / "scan_cold_down.json").read_text()),
    ]:
        if "temperatures" not in record:
            continue
        stage2_unit = []
        for t in np.array(record["temperatures"]):
            mean, swing = fit.calc_pressures(head, parameters, drive, t)
            unit = fit.Parameters(
                packing_constant=(0.0, 1.0), conductance=(0.0, record["ratio"]), frequency=parameters.frequency
            )
            stage2_unit.append(fit.evaluate_head(head, unit, drive, t, mean, swing).heat_flow[-1])
        wall_power = []
        for t in np.array(record["temperatures"]):
            wall = fit.calc_link_heat_flow(
                fit.build_walls(), {fit.WARM_FLANGE: warm, fit.STAGE_1: t[head.stage_1], fit.STAGE_2: t[head.stage_2]}
            )
            wall_power.append(wall[fit.STAGE_2])
        unit = np.array(stage2_unit)
        target = -heater[:, 1] - np.array(wall_power)
        a2 = max(0.0, float(unit @ target / (unit @ unit)))
        profile_rows.append(
            dict(ratio=record["ratio"], a2=a2, rms_W=float(np.sqrt(np.mean((a2 * unit - target) ** 2))))
        )
    best_stage2 = min(profile_rows, key=lambda row: row["rms_W"])
    verification = dict(
        input_and_code_fingerprints_match=True,
        best_stage2_only_fit_on_scanned_cooling_branch=best_stage2,
        log_parameter_order=["a_1", "a_2", "C_1", "C_2"],
        log_parameter_jacobian_W=jacobians[-1].tolist(),
        log_parameter_singular_values_W=singular_values.tolist(),
        log_parameter_condition_number=float(singular_values[0] / singular_values[-1]),
        jacobian_relative_change_1e4_to_1e5=float(
            np.linalg.norm(jacobians[-1] - jacobians[0]) / np.linalg.norm(jacobians[-1])
        ),
        max_tube_density_relative_change_48_to_96=max(quadrature),
        max_stage_power_change_W_48_to_96=max(stage_power_changes),
        max_node_temperature_change_K_48_to_96=max(temperature_changes),
        interpretation="Local sensitivity on the selected branch; no statistical error model is assumed.",
    )
    (folder / "verification.json").write_text(json.dumps(verification, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps({key: value for key, value in verification.items() if key != "log_parameter_jacobian_W"}, indent=2)
    )


if __name__ == "__main__":
    main()
