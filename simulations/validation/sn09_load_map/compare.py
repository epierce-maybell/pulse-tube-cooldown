"""Compares the archived baseline, finite-capacity control, and coupled May fit."""

import json

import fit
import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def curvature(temperature, power):
    """Returns the change in secant slope over three ascending-temperature points, in W/K."""
    order = np.argsort(temperature)
    slopes = np.diff(np.array(power)[order]) / np.diff(np.array(temperature)[order])
    return float(slopes[1] - slopes[0])


def main():
    """Writes the contribution comparison, clamp table, and measured-point overlays."""
    folders = [fit.HERE / "baseline", fit.HERE / "uncoupled", fit.HERE]
    names = ["Baseline", "Finite capacity only", "Finite capacity + exchanger"]
    saved = [json.loads((folder / "results.json").read_text()) for folder in folders]
    terms = [json.loads((folder / "decomposition.json").read_text()) for folder in folders]
    base, control, coupled = saved
    improvement = 100 * (1 - coupled["rms_residual_W"][1] / base["rms_residual_W"][1])
    lines = [
        "SN09 May map: task 6.3",
        "",
        f"Stage-2 RMS falls by {improvement:.2f}% relative to the preserved baseline.",
        "RMS is sqrt(mean(residual^2)) over the same six heater points, in W.",
        "The joint objective sums unweighted squared residuals from both stages.",
        "Both new refits use eight GOS nodes; their difference isolates the exchanger",
        "coupling on that mesh. The baseline used one GOS node, so comparison to it",
        "also includes the resolution change.",
        "",
        "Model                              RMS1/W     RMS2/W       SSE/W^2",
    ]
    for name, result in zip(names, saved, strict=True):
        lines.append(
            f"{name:32s} {result['rms_residual_W'][0]:10.6f} "
            f"{result['rms_residual_W'][1]:10.6f} {result['sum_squared_residual_W2']:13.6f}"
        )
    lines += [
        "",
        f"Fitted exchanger effectiveness eta = {coupled['exchanger_effectiveness']:.9g}.",
        f"Coupled conductances C1, C2 = {coupled['conductances_W_per_Pa2']} W/Pa^2.",
        "C1 reaches its nonnegative bound of zero. This is a boundary fit, with stage-1",
        "cooling supplied by transport from stage 2 in the fitted balances. The smaller",
        "residual establishes a better May fit; it does not establish a transferable",
        "per-unit calibration or reachability from room temperature.",
        "",
        "All archived refits balance their columns with an independently solved drive",
        "and have negative largest real thermal eigenvalues at fixed plate temperatures.",
        "The global minimum across equilibrium branches has not been established.",
        "",
        "Clamp predictions (K), with no external clamp load:",
        " Q1/W Q2/W    baseline    finite only    coupled",
    ]
    for j, point in enumerate(base["points"]):
        lines.append(
            f"{point['heater_power_W'][0]:5.1f} {point['heater_power_W'][1]:4.1f} "
            + " ".join(f"{r['points'][j]['clamp_temperature_K']:12.5f}" for r in saved)
        )
    lines += [
        "",
        "TPR is a pulse-tube thermometer, not a clamp measurement; it is not scored here.",
        "",
        "Stage-2 power decomposition:",
        "Q2 = corrected gross - packing loss - pressure-span correction - wall load.",
        "Corrected gross = beta(Tm)*Tm*H2/(1 + eps_cold) in the coupled fit, with Tm",
        "the matrix minimum; the uncoupled control uses the plate temperature T2.",
        "The pressure-span correction resolves the last matrix boundary's temperature",
        "relative to the cold endpoint. It is zero for the uncoupled plate boundary.",
        "",
        "Curvature diagnostic: high-T2 secant slope minus low-T2 secant slope (W/K).",
        "These connect three measured operating points; T1 also changes within a group.",
        "They are not derivatives at fixed T1.",
        "Model                           Q1 heater/W    Q2 slope change    gross    -packing    -span    -wall",
    ]
    for name, result, decomposition in zip(names, saved, terms, strict=True):
        for indices in ([2, 1, 0], [5, 4, 3]):
            t = [result["points"][j]["stage_temperature_K"][1] for j in indices]
            q = [result["points"][j]["heater_power_W"][1] - result["points"][j]["residual_W"][1] for j in indices]
            gross = [decomposition[j]["corrected_gross_W"] for j in indices]
            packing = [
                decomposition[j].get("cold_packing_loss_W", decomposition[j]["stage2_passive_boundary_W"][-1])
                for j in indices
            ]
            span = [decomposition[j].get("cold_pressure_span_W", 0.0) for j in indices]
            wall = [decomposition[j]["wall2_W"] for j in indices]
            values = [
                curvature(t, q),
                curvature(t, gross),
                -curvature(t, packing),
                -curvature(t, span),
                -curvature(t, wall),
            ]
            lines.append(
                f"{name:32s} {result['points'][indices[0]]['heater_power_W'][0]:6.1f} "
                + " ".join(f"{v:11.6f}" for v in values)
            )
    lines += ["", "GOS resolution:"]
    mesh = json.loads((fit.HERE / "mesh_convergence.json").read_text())
    for row in mesh:
        if "error" in row:
            lines.append(f"{row['gos_entries']} entries: root solve failed at May point {row['failed_point']}.")
        elif "max_power_change_from_previous_W" in row:
            lines.append(
                f"{row['gos_entries']} entries: max fixed-parameter stage-power change from the previous "
                f"grid = {row['max_power_change_from_previous_W']} W."
            )
    lines += [
        "Eight GOS entries resolve the cold-end plateau. The four-entry joint fit gave",
        "RMS1 = 1.143649 W, RMS2 = 0.297175 W, eta = 0.194946; the eight-entry fit gives",
        f"RMS1 = {coupled['rms_residual_W'][0]:.6f} W, RMS2 = {coupled['rms_residual_W'][1]:.6f} W.",
        "The stage-2 improvement survives refinement. Stage-1 fixed-parameter powers",
        "remain more mesh-sensitive; finer-grid branch selection remains a limitation.",
        "",
        "Reproduce from the repository root:",
        "  .venv/Scripts/python scripts/data/derive_head_model.py",
        "  .venv/Scripts/python simulations/validation/sn09_load_map/fit.py",
        "  .venv/Scripts/python simulations/validation/sn09_load_map/fit.py --uncoupled",
        "  .venv/Scripts/python simulations/validation/sn09_load_map/diagnose.py",
        "  .venv/Scripts/python simulations/validation/sn09_load_map/diagnose.py "
        "simulations/validation/sn09_load_map/uncoupled",
        "  .venv/Scripts/python simulations/validation/sn09_load_map/verify.py",
        "  .venv/Scripts/python simulations/validation/sn09_load_map/verify.py "
        "simulations/validation/sn09_load_map/uncoupled",
        "  .venv/Scripts/python simulations/validation/sn09_load_map/explore.py mesh --gos-warm-guess 5.1",
        "  .venv/Scripts/python simulations/validation/sn09_load_map/compare.py",
        "",
        "The baseline is copied unchanged from the pre-6.2 results at commit 5d6f250.",
        "Historical exploratory scans outside baseline/ also predate tasks 6.2 and 6.4.",
        "The candidate remains in this experiment; the August prediction belongs to 2.5.",
    ]
    report = "\n".join(lines) + "\n"
    (fit.HERE / "comparison.txt").write_text(report)
    print(report)
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), layout="constrained")
    colors = ["0.55", "#cc7a00", "#1565aa"]
    for group, indices in enumerate(([2, 1, 0], [5, 4, 3])):
        t = [base["points"][j]["stage_temperature_K"][1] for j in indices]
        for name, result, color in zip(names, saved, colors, strict=True):
            q = [result["points"][j]["heater_power_W"][1] - result["points"][j]["residual_W"][1] for j in indices]
            clamp = [result["points"][j]["clamp_temperature_K"] for j in indices]
            axes[0, group].plot(t, q, "o-", label=name, color=color)
            axes[1, group].plot(t, clamp, "o-", label=name, color=color)
        axes[0, group].plot(
            t, [base["points"][j]["heater_power_W"][1] for j in indices], "s:", color="black", label="Measured heater"
        )
        axes[0, group].set(
            title=f"Stage-1 heater: {base['points'][indices[0]]['heater_power_W'][0]:g} W",
            ylabel="Available stage-2 cooling (W)",
        )
        axes[1, group].set(ylabel="Clamp prediction (K)")
    for ax in axes.flat:
        ax.set(xlabel="Measured stage-2 temperature (K)")
        ax.grid(alpha=0.2)
        ax.legend(fontsize=8)
    fig.suptitle("SN09 May map: baseline and section 6 refits")
    fig.savefig(fit.HERE / "comparison.png", dpi=170)
    plt.close(fig)


if __name__ == "__main__":
    main()
