"""Compare stage priorities on the current May cooling branch.

Run from the repository root with --stage 1, --stage 2, or --stage joint.
The unselected stage has residual weight 1e-4 to select otherwise unconstrained
coefficients. Output records unweighted residuals and independent balance checks.
"""

import argparse
from dataclasses import asdict, replace
import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np
from scipy.optimize import minimize

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import fit


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=("1", "2", "joint"), required=True)
    parser.add_argument("--weight2", type=float)
    parser.add_argument("--exchanger-reference", action="store_true")
    parser.add_argument("--search", choices=("gradient", "simplex"), default="simplex")
    args = parser.parse_args()
    weights = {"1": (1.0, 1e-4), "2": (1e-4, 1.0), "joint": (1.0, 1.0)}[args.stage]
    if args.weight2 is not None:
        weights = (1.0, args.weight2)
    weights = np.array(weights)
    inputs, data, warm, drive = fit.read_experiment()
    head = fit.read_head(inputs["head"])
    saved = json.loads((fit.HERE / "results.json").read_text())
    initial = fit.Parameters(
        packing_constant=tuple(saved["packing_constants"]),
        conductance=tuple(saved["conductances_W_per_Pa2"]),
        exchanger_effectiveness=saved["exchanger_effectiveness"],
        frequency=saved["frequency_Hz"],
    )
    if args.exchanger_reference:
        import pulse_tube_cooldown.library.links.cooler.head as physics

        original_boundary = physics._stage_2_boundary_temperature

        def boundary(temperature, stage_1):
            result = original_boundary(temperature, stage_1)
            result[-1] = temperature[-1]
            return result

        physics._stage_2_boundary_temperature = boundary

    original_nnls = fit.nnls
    original_fit_at_ratio = fit.fit_at_ratio

    def weighted_nnls(matrix, target):
        scale = np.tile(weights, len(target) // 2)
        return original_nnls(matrix * scale[:, None], target * scale)

    evaluated = []

    def weighted_fit_at_ratio(*positional, **keywords):
        result = original_fit_at_ratio(*positional, **keywords)
        evaluated.append(result)
        return replace(result, residual=result.residual * weights)

    fit.nnls = weighted_nnls
    fit.fit_at_ratio = weighted_fit_at_ratio
    start = time.perf_counter()
    failure = None
    try:
        if args.search == "gradient":
            fitted = fit.fit_load_map(
                head, drive, fit.build_walls(), data,
                warm_temperature=warm, initial=initial,
                ratio_bounds=(0.3e-9, 3.0e-9),
                guesses=[np.array(row["temperature_K"]) for row in saved["points"]],
                verbose=True,
            )
            result = replace(fitted, residual=fitted.residual / weights)
        else:
            cache = {}
            ratio0 = initial.conductance[1] / initial.packing_constant[1]
            end0 = initial.exchanger_effectiveness / initial.packing_constant[1]

            def objective(x):
                if tuple(x) in cache:
                    return float(np.sum((cache[tuple(x)].residual * weights) ** 2))
                ratio, end = np.exp(x) * [ratio0, end0]
                if not (0.3e-9 <= ratio <= 3e-9 and 1e-3 <= end <= 100):
                    return 1e6 + float(x @ x)
                guesses = ([c.temperature for c in cache[min(cache, key=lambda k: np.linalg.norm(x-k))].columns]
                           if cache else [np.array(row["temperature_K"]) for row in saved["points"]])
                try:
                    candidate = original_fit_at_ratio(
                        head, drive, fit.build_walls(), data, warm_temperature=warm,
                        ratio=ratio, effectiveness_ratio=end, frequency=initial.frequency, guesses=guesses,
                    )
                except RuntimeError:
                    print("FAILED", x.tolist(), flush=True)
                    return 1e6 + float(x @ x)
                cache[tuple(x)] = candidate
                evaluated.append(candidate)
                score = float(np.sum((candidate.residual * weights) ** 2))
                print("FIT", x.tolist(), np.sqrt(np.mean(candidate.residual**2, axis=0)), flush=True)
                return score

            optimum = minimize(objective, np.zeros(2), method="Nelder-Mead", options={
                "initial_simplex": np.array([[0., 0.], [0.02, 0.], [0., 0.02]]),
                "maxiter": 70, "xatol": 1e-3, "fatol": 1e-6,
            })
            if not optimum.success:
                failure = str(optimum.message)
            result = min(evaluated, key=lambda r: np.sum((r.residual * weights) ** 2))
    except RuntimeError as error:
        failure = str(error)
        if not evaluated:
            raise
        result = min(evaluated, key=lambda r: np.sum((r.residual * weights) ** 2))
    rows = []
    for index, column in enumerate(result.columns):
        stage = np.array([data.T1[index].item(), data.T2[index].item()])
        heater = np.array([data.Q1[index].item(), data.Q2[index].item()])
        mean, swing = fit.calc_pressures(head, result.parameters, drive, column.temperature)
        verified = fit.evaluate_head(head, result.parameters, drive, column.temperature, mean, swing)
        walls = fit.calc_link_heat_flow(fit.build_walls(), {fit.WARM_FLANGE: warm, fit.STAGE_1: stage[0], fit.STAGE_2: stage[1]})
        residual = verified.heat_flow[[head.stage_1, head.stage_2]] + heater + [walls[fit.STAGE_1], walls[fit.STAGE_2]]
        np.testing.assert_allclose(residual, result.residual[index], atol=1e-8)
        rows.append(dict(
            heater_W=heater.tolist(), plate_K=stage.tolist(),
            residual_W=residual.tolist(), temperature_K=column.temperature.tolist(),
            mean_pressure_Pa=mean, swing_Pa=swing,
            column_balance_W=float(np.max(abs(verified.heat_flow[fit.column_indices(head)]))),
            growth_rate_per_s=fit.calc_growth_rate(head, result.parameters, drive, verified),
        ))
    output = dict(
        stage=args.stage, residual_weights=weights.tolist(),
        exchanger_reference=args.exchanger_reference,
        search=args.search,
        optimizer_failure=failure, elapsed_s=time.perf_counter()-start,
        parameters=asdict(result.parameters),
        rms_W=np.sqrt(np.mean(result.residual**2, axis=0)).tolist(),
        points=rows,
        source_sha256={str(p.relative_to(fit.HERE)): hashlib.sha256(p.read_bytes()).hexdigest()
                       for p in (fit.HERE / "fit.py", fit.HERE / "inputs.toml", fit.HERE / "results.json")},
    )
    suffix = f"weight2_{args.weight2:g}" if args.weight2 is not None else args.stage
    if args.exchanger_reference:
        suffix += "_exchanger"
    (HERE / f"fit_{suffix}.json").write_text(json.dumps(output, indent=2) + "\n")
    print("FINAL", json.dumps({k:v for k,v in output.items() if k not in ("points", "source_sha256")}), flush=True)


if __name__ == "__main__":
    main()
