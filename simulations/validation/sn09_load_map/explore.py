"""Seeds, scans, and refines the May fit using the current head equations."""

import argparse
import json
from dataclasses import replace

import fit
import numpy as np


def refine_gos(head, count):
    """Divides the complete cold GOS layer into equal-volume entries."""
    warm = next(i for i, material in enumerate(head.material[1]) if material.name == "Gd2O2S")
    return replace(
        head,
        material=(head.material[0], head.material[1][:warm] + (head.material[1][warm],) * count),
        mass=(head.mass[0], np.r_[head.mass[1][:warm], np.full(count, head.mass[1][warm:].sum() / count)]),
        pore_volume=(
            head.pore_volume[0],
            np.r_[head.pore_volume[1][:warm], np.full(count, head.pore_volume[1][warm:].sum() / count)],
        ),
    )


def interpolate_profile(old, new, temperature):
    """Interpolates the GOS entry centers while retaining the plates and other materials."""
    first = next(i for i, material in enumerate(old.material[1]) if material.name == "Gd2O2S") + old.stage_1 + 1
    old_count, new_count = old.stage_2 - first, new.stage_2 - first
    gos = np.interp(
        (np.arange(new_count) + 0.5) / new_count, (np.arange(old_count) + 0.5) / old_count, temperature[first:-1]
    )
    return np.r_[temperature[:first], gos, temperature[-1]]


def parameters(saved):
    """Reads the five fitted parameters from a results record."""
    return fit.Parameters(
        packing_constant=tuple(saved["packing_constants"]),
        conductance=tuple(saved["conductances_W_per_Pa2"]),
        exchanger_effectiveness=saved["exchanger_effectiveness"],
        frequency=saved["frequency_Hz"],
    )


def mesh(folder=fit.HERE, counts=(4, 8, 16), gos_warm_guess=None):
    """Checks fixed-parameter stage powers and profiles across GOS entry counts."""
    saved = json.loads((folder / "results.json").read_text())
    inputs, _, warm, drive = fit.read_experiment()
    original = replace(fit.read_head(inputs["head"]), cold_exchanger_coupling=saved["cold_exchanger_coupling"])
    p = parameters(saved)
    old = original
    templates = [np.array(row["temperature_K"]) for row in saved["points"]]
    rows = []
    for count in counts:
        head = refine_gos(original, count)
        columns = []
        for point, t in zip(saved["points"], templates, strict=True):
            guess = interpolate_profile(old, head, t)
            if gos_warm_guess is not None:
                first = next(i for i, m in enumerate(head.material[1]) if m.name == "Gd2O2S") + head.stage_1 + 1
                guess[first] = gos_warm_guess
            try:
                columns.append(
                    fit.solve_columns(
                        head,
                        p,
                        drive,
                        warm_temperature=warm,
                        stage_temperature=tuple(point["stage_temperature_K"]),
                        guess=guess,
                        max_nfev=200,
                    )
                )
            except RuntimeError as error:
                rows.append(dict(gos_entries=count, failed_point=point["point"], error=str(error)))
                (folder / "mesh_convergence.json").write_text(json.dumps(rows, indent=2) + "\n")
                print(rows[-1], flush=True)
                return
        powers = np.array([c.heat_flow[[head.stage_1, head.stage_2]] for c in columns])
        row = dict(
            gos_entries=count,
            stage_heat_flow_W=powers.tolist(),
            temperatures_K=[c.temperature.tolist() for c in columns],
            clamp_temperature_K=[float(c.temperature[head.clamp]) for c in columns],
            cold_temperature_K=[float(c.temperature[-2]) for c in columns],
        )
        row["growth_rate_per_s"] = [fit.calc_growth_rate(head, p, drive, c) for c in columns]
        row["gos_warm_guess_K"] = gos_warm_guess
        if rows:
            row["max_power_change_from_previous_W"] = np.max(
                abs(powers - np.array(rows[-1]["stage_heat_flow_W"])), axis=0
            ).tolist()
        rows.append(row)
        print(
            json.dumps({k: v for k, v in row.items() if k not in ("temperatures_K", "stage_heat_flow_W")}), flush=True
        )
        (folder / "mesh_convergence.json").write_text(json.dumps(rows, indent=2) + "\n")
        old, templates = head, [c.temperature for c in columns]


def seed():
    """Resolves a cold branch from the baseline's warm layers and a cold GOS guess."""
    inputs, lm, warm, drive = fit.read_experiment()
    head = fit.read_head(inputs["head"])
    baseline = json.loads((fit.HERE / "baseline/results.json").read_text())
    p = fit.Parameters(
        packing_constant=tuple(inputs["initial"]["packing_constant"]),
        conductance=tuple(inputs["initial"]["conductance"]),
        exchanger_effectiveness=inputs["initial"]["exchanger_effectiveness"],
        frequency=inputs["frequency"],
    )
    count = head.stage_2 - 7
    columns = []
    for row in baseline["points"]:
        t = np.array(row["temperature_K"])
        guess = np.r_[t[:7], np.linspace(5.1, t[-1] - 0.3, count), t[-1]]
        columns.append(
            fit.solve_columns(
                head,
                p,
                drive,
                warm_temperature=warm,
                stage_temperature=tuple(row["stage_temperature_K"]),
                guess=guess,
                max_nfev=500,
            )
        )
    output = dict(
        source="explore.py seed; baseline warm layers, linear GOS guess from 5.1 K to T2 - 0.3 K",
        temperatures=[c.temperature.tolist() for c in columns],
    )
    (fit.HERE / "cooling_seed.json").write_text(json.dumps(output, indent=2) + "\n")


def scan(ratios, effectiveness_ratio):
    """Continues a cold branch in C_2/a_2 at an explicit eta/a_2."""
    inputs, lm, warm, drive = fit.read_experiment()
    head = fit.read_head(inputs["head"])
    guesses = [np.array(t) for t in json.loads((fit.HERE / "cooling_seed.json").read_text())["temperatures"]]
    rows = []
    for ratio in ratios:
        result = fit.fit_at_ratio(
            head,
            drive,
            fit.build_walls(),
            lm,
            warm_temperature=warm,
            ratio=ratio,
            effectiveness_ratio=effectiveness_ratio,
            frequency=inputs["frequency"],
            guesses=guesses,
        )
        guesses = [c.temperature for c in result.columns]
        row = dict(
            ratio=ratio,
            effectiveness_ratio=effectiveness_ratio,
            rms_W=np.sqrt(np.mean(result.residual**2, axis=0)).tolist(),
            eta=result.parameters.exchanger_effectiveness,
            temperatures=[t.tolist() for t in guesses],
            growth_rate_per_s=[fit.calc_growth_rate(head, result.parameters, drive, c) for c in result.columns],
        )
        rows.append(row)
        print(json.dumps({k: v for k, v in row.items() if k != "temperatures"}), flush=True)
        (fit.HERE / "coupled_scan.json").write_text(json.dumps(rows, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("seed", "mesh", "scan"))
    parser.add_argument("--counts", nargs="+", type=int, default=[4, 8, 16])
    parser.add_argument("--gos-warm-guess", type=float)
    parser.add_argument("--ratios", nargs="+", type=float)
    parser.add_argument("--effectiveness-ratio", type=float, default=2.5)
    args = parser.parse_args()
    if args.action == "seed":
        seed()
    elif args.action == "mesh":
        mesh(counts=args.counts, gos_warm_guess=args.gos_warm_guess)
    else:
        scan(args.ratios, args.effectiveness_ratio)
