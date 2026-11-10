"""Plot and summarize saved candidate predictions against matched measurements."""

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

from simulations.validation.head_comparison.comparison import cell_derivatives, load_population


def main() -> None:
    """Write per-head diagnostics without rerunning or changing the fit."""
    here = Path(__file__).resolve().parent
    results = here / "results"
    output = results / "diagnostics"
    output.mkdir(exist_ok=True)
    heads = load_population(here.parents[1] / "validation/head_comparison/records.toml", root=here.parents[2])
    plt.switch_backend("Agg")
    summaries = []
    for head in heads:
        serial = head.serial
        cooldown = xr.load_dataset(results / f"{serial}_cooldown.nc")
        predicted = xr.load_dataset(results / f"{serial}_map.nc")
        measured = head.load_map.data
        best = json.loads((results / f"{serial}_best.json").read_text())
        optimizer_path = results / f"{serial}_optimizer.json"
        optimizer = json.loads(optimizer_path.read_text()) if optimizer_path.exists() else None
        evaluations = [json.loads(line) for line in (results / f"{serial}_evaluations.jsonl").read_text().splitlines()]
        summary = {
            "head": serial,
            "optimizer": optimizer,
            "logged_calls": len(evaluations),
            "failed_calls": sum("error" in row for row in evaluations),
            "best_evaluation": best["evaluation"],
            "best_score": best["score"],
            "blocks": best["blocks"],
        }
        figure, axes = plt.subplots(4, 2, figsize=(14, 16))
        figure.subplots_adjust(left=0.08, right=0.97, bottom=0.08, top=0.94, hspace=0.55, wspace=0.28)
        figure.suptitle(
            f"{serial}: saved candidate vs measurements | {head.cooldown.compressors} compressor(s)", fontsize=16
        )
        for column, stage in enumerate(("T1", "T2")):
            ax = axes[0, column]
            observed = head.cooldown.data[stage]
            if stage == "T1":
                observed = observed.where(observed < head.stage_1_max_temperature)
            ax.plot(head.cooldown.data.time / 3600, observed, color="black", label="Measured")
            ax.plot(cooldown.time / 3600, cooldown[stage], color="#c44e28", ls="--", label="Candidate")
            ax.set(title=f"Stage {column + 1} cooldown", xlabel="Time since start (h)", ylabel=f"{stage} (K)")
            if stage == "T2":
                ax.set_yscale("log")
            ax.legend()
            difference = predicted[stage].values - measured[stage].values
            summary[f"map_{stage}_error_K"] = {
                "min": float(difference.min()),
                "mean": float(difference.mean()),
                "max": float(difference.max()),
            }
            summary[f"cooldown_{stage}_final_K"] = float(cooldown[stage][-1])
            summary[f"measured_{stage}_last_valid_K"] = float(observed.dropna("time")[-1])
            ax = axes[1, column]
            for index, power in enumerate(np.unique(measured.Q2)):
                color = plt.get_cmap("tab10")(index % 10)
                for points, linestyle, face in ((measured, "-", color), (predicted, "--", "white")):
                    selection = points.where(points.Q2 == power, drop=True).sortby("Q1")
                    ax.plot(
                        selection.Q1,
                        selection[stage],
                        linestyle=linestyle,
                        marker="o",
                        color=color,
                        markerfacecolor=face,
                        label=f"Stage 2 heater: {power:g} W" if points is measured else None,
                    )
            ax.set(title=f"Stage {column + 1} load map", xlabel="Stage 1 heater power (W)", ylabel=f"{stage} (K)")
            ax.legend(fontsize=8, ncols=2)
        minimum = int(np.argmin(cooldown.T1.values))
        summary["stage_1_trough"] = {
            "time_h": float(cooldown.time[minimum]) / 3600,
            "temperature_K": float(cooldown.T1[minimum]),
            "rise_to_final_K": float(cooldown.T1[-1] - cooldown.T1[minimum]),
        }
        measured_cells, predicted_cells = cell_derivatives(measured), cell_derivatives(predicted)
        labels = [
            f"{a:g}-{b:g}\n{c:g}-{d:g}"
            for a, b, c, d in zip(
                measured_cells.P1_low.values,
                measured_cells.P1_high.values,
                measured_cells.P2_low.values,
                measured_cells.P2_high.values,
                strict=True,
            )
        ]
        positions = np.arange(len(labels))
        for ax, name in zip(axes[2:].flat, ("dQ1dT1", "dQ1dT2", "dQ2dT1", "dQ2dT2"), strict=True):
            observed, estimate = measured_cells[name].values, predicted_cells[name].values
            ax.vlines(positions, observed, estimate, color="0.7", lw=1.5)
            ax.scatter(positions, observed, color="black", label="Measured", zorder=3)
            ax.scatter(positions, estimate, facecolors="white", edgecolors="#c44e28", label="Candidate", zorder=3)
            if name in ("dQ1dT2", "dQ2dT1"):
                ax.axhline(0, color="0.6", lw=0.7)
            ax.set_xticks(positions, labels, fontsize=8)
            ax.set(
                title=f"dQ{name[2]}/dT{name[-1]}: matched heater cells",
                ylabel="Heat-flow slope (W/K)",
                xlabel="Heater intervals (W): stage 1 / stage 2",
            )
            ax.legend(fontsize=9)
            summary[name] = {"measured": observed.tolist(), "candidate": estimate.tolist()}
        for ax in axes.flat:
            ax.grid(alpha=0.2)
            ax.set_axisbelow(True)
        figure.supxlabel(
            "Load maps: solid / filled = measured; dashed / open = candidate. Q in slopes is heat into the stage.",
            fontsize=11,
        )
        figure.savefig(output / f"{serial}.png", dpi=150)
        plt.close(figure)
        summaries.append(summary)
        print(
            json.dumps(
                {
                    key: value
                    for key, value in summary.items()
                    if key not in ("optimizer", "blocks", "dQ1dT1", "dQ1dT2", "dQ2dT1", "dQ2dT2")
                }
            )
        )
    (output / "summary.json").write_text(json.dumps(summaries, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
