"""Exploratory branch searches for the SN09 load-map experiment.

Run with warm, cooling-up, cooling-down, low, roots, low-roots, or mixed.
JSON records retain converged roots and the first failed continuation step.
"""

import argparse
import itertools
import json
from pathlib import Path

import fit
import numpy as np
from scipy.optimize import nnls


def find_roots() -> None:
    """Searches several initial profiles at the first and last May points."""
    inputs, lm, warm, drive = fit.read_experiment()
    head = fit.read_head(inputs["head"])
    rows = []
    for point in [0, 5]:
        stages = (float(lm.T1[point]), float(lm.T2[point]))
        for ratio in [5e-10, 1e-9, 3e-9, 1e-8]:
            params = fit.Parameters(
                packing_constant=(0.15, 0.01), conductance=(2.5e-10, 0.01 * ratio), frequency=inputs["frequency"]
            )
            roots = []
            for temps in [
                [stages[0] * 0.9, stages[0] * 0.6, 8, 4.5],
                [stages[0] * 0.6, 10, 5, 4],
                [25, 8, 4, 3],
                [stages[0] - 1, stages[0] - 2, 15, 5.1],
                [10, 7, 4, 2.15],
            ]:
                guess = np.array([warm, 220, 90, stages[0], *temps, stages[1]])
                try:
                    c = fit.solve_columns(
                        head, params, drive, warm_temperature=warm, stage_temperature=stages, guess=guess, max_nfev=120
                    )
                    if any(np.max(abs(c.temperature - r.temperature)) < 1e-3 for r in roots):
                        continue
                    roots.append(c)
                    growth = fit.calc_growth_rate(head, params, drive, c)
                    row = dict(
                        point=point,
                        ratio=ratio,
                        temperature=c.temperature.tolist(),
                        growth=growth,
                        stage_flow=c.heat_flow[[3, 8]].tolist(),
                    )
                    rows.append(row)
                    print("ROOT", row, flush=True)
                except RuntimeError as error:
                    print("FAIL", point, ratio, temps, str(error), flush=True)
    Path(__file__).with_name("branches.json").write_text(json.dumps(rows, indent=2))


def find_low_roots() -> None:
    """Continues distinct low-temperature roots through all May points."""
    inputs, lm, warm, drive = fit.read_experiment()
    head = fit.read_head(inputs["head"])
    params = fit.Parameters(packing_constant=(1.0, 1.0), conductance=(0.0, 1e-9), frequency=inputs["frequency"])
    roots = json.loads(Path(__file__).with_name("branches.json").read_text())
    out = []
    for point in [0, 1, 2, 5, 4, 3]:
        t1, t2 = float(lm.T1[point]), float(lm.T2[point])
        found = []
        seeds = [r["temperature"] for r in roots if r["ratio"] == 1e-9 and r["growth"] < 0]
        seeds.extend([r["temperature"] for r in out])
        for seed in seeds:
            guess = np.array(seed)
            guess[[0, 3, 8]] = [warm, t1, t2]
            try:
                c = fit.solve_columns(
                    head, params, drive, warm_temperature=warm, stage_temperature=(t1, t2), guess=guess
                )
                if any(max(abs(c.temperature - np.array(r["temperature"]))) < 1e-3 for r in found):
                    continue
                p = fit.Parameters(
                    packing_constant=(0.15, 0.01), conductance=(2.5e-10, 1e-11), frequency=inputs["frequency"]
                )
                growth = fit.calc_growth_rate(head, p, drive, c)
                row = dict(point=point, ratio=1e-9, growth=growth, temperature=c.temperature.tolist())
                found.append(row)
                out.append(row)
                print("ROOT", row, flush=True)
            except RuntimeError as e:
                print("FAIL", point, str(e), flush=True)
    Path(__file__).with_name("low_branch.json").write_text(json.dumps(out, indent=2))


def compare_branches() -> None:
    """Compares all positive hot/cold branch assignments on an overlapping ratio interval."""
    inputs, lm, warm, drive = fit.read_experiment()
    head = fit.read_head(inputs["head"])
    hot = json.loads(Path(__file__).with_name("scan.json").read_text())
    cold = json.loads(Path(__file__).with_name("scan_cold_down.json").read_text())
    rows = []
    for ratio in np.linspace(4.1e-10, 5.2e-10, 9):
        matrices = []
        profiles = []
        targets = []
        for point, (t1, t2) in enumerate(zip(lm.T1.values, lm.T2.values, strict=True)):
            variants = []
            states = []
            for branch, data in [("hot", hot), ("cold", cold)]:
                seed = min((r for r in data if "temperatures" in r), key=lambda r: abs(np.log(r["ratio"] / ratio)))
                p = fit.Parameters(packing_constant=(1.0, 1.0), conductance=(0.0, ratio), frequency=inputs["frequency"])
                try:
                    col = fit.solve_columns(
                        head,
                        p,
                        drive,
                        warm_temperature=warm,
                        stage_temperature=(t1, t2),
                        guess=np.array(seed["temperatures"][point]),
                    )
                except RuntimeError:
                    continue
                mean = (col.supply_pressure + col.return_pressure) / 2
                sw = col.supply_pressure - col.return_pressure
                basis = [
                    fit.Parameters(packing_constant=(1.0, 0.0), conductance=(0.0, 0.0), frequency=inputs["frequency"]),
                    fit.Parameters(
                        packing_constant=(0.0, 1.0), conductance=(0.0, ratio), frequency=inputs["frequency"]
                    ),
                    fit.Parameters(
                        packing_constant=(0.0, 0.0), conductance=(1e-10, 0.0), frequency=inputs["frequency"]
                    ),
                ]
                matrix = np.column_stack(
                    [fit.evaluate_head(head, p, drive, col.temperature, mean, sw).heat_flow[[3, 8]] for p in basis]
                )
                variants.append(matrix)
                states.append((branch, col))
            matrices.append(variants)
            profiles.append(states)
            wall = fit.calc_link_heat_flow(fit.build_walls(), {fit.WARM_FLANGE: warm, fit.STAGE_1: t1, fit.STAGE_2: t2})
            targets.extend([-float(lm.Q1[point]) - wall[fit.STAGE_1], -float(lm.Q2[point]) - wall[fit.STAGE_2]])
        best = None
        for assignment in itertools.product(*(range(len(v)) for v in matrices)):
            matrix = np.concatenate([variants[choice] for variants, choice in zip(matrices, assignment, strict=True)])
            coeff, norm = nnls(matrix, np.array(targets))
            if min(coeff) > 0 and (best is None or norm < best[0]):
                best = (norm, coeff, assignment, matrix)
        norm, co, assign, matrix = best
        a1, a2, c1 = co
        params = fit.Parameters(
            packing_constant=(a1, a2), conductance=(c1 * 1e-10, ratio * a2), frequency=inputs["frequency"]
        )
        chosen = [states[choice] for states, choice in zip(profiles, assign, strict=True)]
        growth = [fit.calc_growth_rate(head, params, drive, c) for _, c in chosen]
        row = dict(
            ratio=float(ratio),
            cost=norm**2,
            parameters=[a1, a2, c1 * 1e-10, ratio * a2],
            branches=[b for b, _ in chosen],
            growth=growth,
            residual=(matrix @ co - np.array(targets)).reshape(-1, 2).tolist(),
            temperatures=[c.temperature.tolist() for _, c in chosen],
        )
        rows.append(row)
        print("FIT", ratio, row["cost"], row["parameters"], row["branches"], flush=True)
    Path(__file__).with_name("mixed_branches.json").write_text(json.dumps(rows, indent=2))


def scan_branch(name: str) -> None:
    """Continues one branch in ratio and records the associated nonnegative fits."""
    inputs, load_map, warm, drive = fit.read_experiment()
    head = fit.read_head(inputs["head"])
    if name == "warm":
        ratios = np.geomspace(1e-13, 3e-9, 65)
        previous = None
        output = "scan.json"
    elif name == "low":
        ratios = np.geomspace(1e-9, 3e-10, 30)
        seeds = json.loads(Path(__file__).with_name("low_branch.json").read_text())
        previous = [np.array(next(row["temperature"] for row in seeds if row["point"] == i)) for i in range(6)]
        output = "scan_low.json"
    else:
        seeds = json.loads(Path(__file__).with_name("cooling_seed.json").read_text())
        previous = [np.array(t) for t in seeds["temperatures"]]
        if name == "cooling-up":
            ratios = np.geomspace(5e-10, 3e-9, 30)
            output = "scan_cold.json"
        else:
            ratios = np.geomspace(5e-10, 1e-10, 25)
            output = "scan_cold_down.json"
    rows = []
    for ratio in ratios:
        try:
            result = fit.fit_at_ratio(
                head,
                drive,
                fit.build_walls(),
                load_map,
                warm_temperature=warm,
                ratio=ratio,
                frequency=inputs["frequency"],
                guesses=previous,
            )
            previous = [column.temperature for column in result.columns]
            growth = [fit.calc_growth_rate(head, result.parameters, drive, c) for c in result.columns]
            row = dict(
                ratio=ratio,
                parameters=[*result.parameters.packing_constant, *result.parameters.conductance],
                cost=float(np.sum(result.residual**2)),
                residual=result.residual.tolist(),
                growth=growth,
                temperatures=[c.temperature.tolist() for c in result.columns],
            )
            rows.append(row)
            print("FIT", ratio, row["cost"], row["parameters"], "growth", growth, flush=True)
        except RuntimeError as error:
            print("FAILED", ratio, str(error), flush=True)
            rows.append(dict(ratio=ratio, error=str(error)))
            break
    Path(__file__).with_name(output).write_text(json.dumps(rows, indent=2) + "\n", encoding="utf-8")


def main() -> None:
    """Runs one explicitly selected exploratory search."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("search", choices=["warm", "cooling-up", "cooling-down", "low", "roots", "low-roots", "mixed"])
    arguments = parser.parse_args()
    searches = {"roots": find_roots, "low-roots": find_low_roots, "mixed": compare_branches}
    if arguments.search in searches:
        searches[arguments.search]()
    else:
        scan_branch(arguments.search)


if __name__ == "__main__":
    main()
