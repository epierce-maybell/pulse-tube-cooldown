"""Draws the load maps read as heater power against stage temperatures, Q1(T1, T2) and Q2(T1, T2).

Grids are measured at fixed heater powers. A cut at a fixed temperature of one stage interpolates linearly along each
line of the grid that holds the other stage's power fixed, without extrapolating, so each line contributes one point
to the cut. Slopes come from the cell partials of `characterize.py`; power laws are fitted to them in log-log, pooled
over runs and run by run. Figures go to `figures/` and the fitted numbers to `functional_forms.txt`.
"""

import matplotlib.pyplot as plt
import numpy as np
from characterize import HERE, MUTED, cell_partials, collapse, collapse_groups, grids_of, levels, save, value
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.figure import Figure
from numpy.typing import NDArray
from runs import Run, load_runs

RAMP = LinearSegmentedColormap.from_list("blue", ["#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#104281", "#0d366b"])
GROUP_COLOR = {"one compressor": "#2a78d6", "two compressors": "#eb6834"}
T1_SCALE = Normalize(35, 70)
T2_SCALE = Normalize(2.0, 5.0)


def grid_lines(run: Run, stage: int) -> list[NDArray[np.float64]]:
    """Returns the lines of a grid that hold the other stage's power fixed.

    Each line is an array of rows (own power, T1, T2), ordered by own power, with missing temperatures dropped.
    """
    other = 3 - stage
    lines = []
    for fixed in levels(getattr(run, f"q{other}")):
        powers = levels(getattr(run, f"q{stage}")[getattr(run, f"q{other}") == fixed])
        pairs = [(p, fixed) if stage == 1 else (fixed, p) for p in powers]
        points = [(p, value(run, *pair, "t1"), value(run, *pair, "t2")) for p, pair in zip(powers, pairs, strict=True)]
        points = [point for point in points if not np.isnan(point[1:]).any()]
        if len(points) >= 2:
            lines.append(np.array(points))
    return lines


def cuts(run: Run, stage: int, spacing: float) -> list[tuple[float, NDArray[np.float64]]]:
    """Returns cuts through a grid at fixed temperatures of one stage.

    Args:
        run: A grid.
        stage: The stage whose temperature each cut holds fixed.
        spacing: Spacing of the fixed temperatures in K.

    Returns:
        For each fixed temperature crossed by at least two lines, the rows (other stage's temperature, own power),
        ordered by the other stage's temperature.
    """
    lines = grid_lines(run, stage)
    if not lines:
        return []
    own = np.concatenate([line[:, stage] for line in lines])
    targets = np.arange(np.ceil(own.min() / spacing) * spacing, own.max() + 1e-9, spacing)
    out = []
    for target in targets:
        points = []
        for line in lines:
            temps = line[:, stage]
            if np.all(np.diff(temps) > 0) and temps[0] <= target <= temps[-1]:
                points.append((np.interp(target, temps, line[:, 3 - stage]), np.interp(target, temps, line[:, 0])))
        if len(points) >= 2:
            out.append((float(target), np.array(sorted(points))))
    return out


def power_law(temperature: NDArray[np.float64], slope: NDArray[np.float64]) -> tuple[float, float]:
    """Fits slope = prefactor * temperature^exponent in log-log to positive slopes; returns (exponent, prefactor)."""
    keep = slope > 0
    exponent, intercept = np.polyfit(np.log(temperature[keep]), np.log(slope[keep]), 1)
    return float(exponent), float(np.exp(intercept))


def run_cells(run: Run, x_key: str, y_key: str, span: float) -> tuple[NDArray, NDArray] | None:
    """Returns a run's cell temperatures and partials if the cells span at least `span` K, else None."""
    cells = cell_partials(run)
    x, y = np.array([c[x_key] for c in cells]), np.array([c[y_key] for c in cells])
    return (x, y) if x.size >= 3 and np.ptp(x) >= span else None


def colorbar(fig: Figure, ax: plt.Axes | list, scale: Normalize, label: str) -> None:
    """Adds a colorbar for the blue temperature ramp."""
    fig.colorbar(plt.cm.ScalarMappable(norm=scale, cmap=RAMP), ax=ax, label=label, shrink=0.85, pad=0.02)


def compressor_groups(runs: list[Run]) -> dict[str, list[Run]]:
    """Returns the grids on one compressor and on two."""
    grids = grids_of(runs)
    return {
        "one compressor": [r for r in grids if r.compressors == 1],
        "two compressors": [r for r in grids if r.compressors == 2],
    }


def plot_slope_fits(ax: plt.Axes, runs: list[Run], x_key: str, y_key: str, span: float, notes: list[str]) -> None:
    """Draws cell partials in log-log with the pooled power law per compressor group and each run's own power law."""
    for name, group in compressor_groups(runs).items():
        color = GROUP_COLOR[name]
        cells = [c for run in group for c in cell_partials(run)]
        x, y = np.array([c[x_key] for c in cells]), np.array([c[y_key] for c in cells])
        ax.scatter(x, y, s=12, color=color, alpha=0.35, linewidths=0)
        exponent, prefactor = power_law(x, y)
        line = np.linspace(x.min(), x.max(), 50)
        ax.plot(line, prefactor * line**exponent, color=color, linewidth=2.2, label=f"{name}, pooled: ^{exponent:.2f}")
        own = []
        for run in group:
            fitted = run_cells(run, x_key, y_key, span)
            if fitted:
                run_exponent, run_prefactor = power_law(*fitted)
                line = np.linspace(fitted[0].min(), fitted[0].max(), 20)
                ax.plot(line, run_prefactor * line**run_exponent, color=color, linewidth=0.9, linestyle="--")
                own.append(f"{run.label} {run_exponent:.2f}")
        notes.append(
            f"{y_key} against {x_key.upper()}, {name}: pooled {prefactor:.4g} W/K * ({x_key.upper()} / 1 K)^"
            f"{exponent:.2f} over {y.size} cells; per run: " + ("; ".join(own) or "none span enough")
        )
    ax.plot([], [], color=MUTED, linestyle="--", linewidth=0.9, label="single runs")
    ax.set(xscale="log", yscale="log")
    ax.legend(loc="best")


def plot_stage2(runs: list[Run], notes: list[str]) -> None:
    """Q2 against T2 after each run's constant, colored by T1; Q2 at fixed T2 against T1; dQ2/dT2 in log-log."""
    fig, axes = plt.subplots(1, 4, figsize=(22, 5), gridspec_kw={"wspace": 0.35})
    groups = collapse_groups(runs)
    for ax, name in zip(axes[:2], ["every variant, one compressor", "every variant, two compressors"], strict=True):
        fit = collapse(groups[name], stage=2)
        order = np.argsort(fit["temperature"])
        t1 = np.concatenate([run.t1[~np.isnan(run.t1) & ~np.isnan(run.t2)] for run in groups[name]])
        ax.scatter(fit["temperature"], fit["power"] + fit["constant"][fit["run"]], c=t1, cmap=RAMP, norm=T1_SCALE, s=16)
        x = fit["temperature"][order] - 2.0
        ax.plot(x + 2.0, fit["shared"][0] * x + fit["shared"][1] * x**2, color="#0b0b0b", linewidth=1.0)
        ax.set(
            title=f"Q2 + each run's constant against T2,\n{name}; black: shared curve",
            xlabel="T2 (K)",
            ylabel="Q2 + run constant (W)",
        )
    colorbar(fig, axes[1], T1_SCALE, "T1 (K)")
    for run in grids_of(runs):
        for target, points in cuts(run, stage=2, spacing=0.25):
            axes[2].plot(
                points[:, 0],
                points[:, 1] - points[0, 1],
                color=RAMP(T2_SCALE(target)),
                linewidth=0.9,
                alpha=0.8,
                linestyle="-" if run.compressors == 1 else "--",
            )
    axes[2].set(
        title="Q2 at fixed T2 against T1, relative to the cut's lowest T1\n(solid one compressor, dashed two)",
        xlabel="T1 (K)",
        ylabel="Q2 − Q2 at the cut's lowest T1 (W)",
        ylim=(-0.5, 0.5),
    )
    colorbar(fig, axes[2], T2_SCALE, "T2 of the cut (K)")
    plot_slope_fits(axes[3], runs, "t2", "dq2_dt2", span=1.0, notes=notes)
    axes[3].set(title="dQ2/dT2 at fixed T1, log-log", xlabel="T2 (K)", ylabel="dQ2/dT2 (W/K)")
    axes[3].set_xticks([2.5, 3, 4, 5], ["2.5", "3", "4", "5"])
    save(fig, "12_q2_against_t1_t2.png")


def small_multiples(count: int, columns: int = 5) -> tuple[Figure, NDArray]:
    """Returns a figure with a grid of shared-axis panels and the unused ones hidden."""
    rows = int(np.ceil(count / columns))
    fig, axes = plt.subplots(rows, columns, figsize=(3.6 * columns, 2.9 * rows), sharex=True, sharey=True)
    for ax in axes.flat[count:]:
        ax.set_axis_off()
    return fig, axes


def plot_stage1_against_t1(runs: list[Run]) -> None:
    """Q1 against T1 with Q2 held at each level, colored by T2, for grids with at least three stage-1 powers."""
    chosen = [r for r in grids_of(runs) if len(levels(r.q1)) >= 3]
    fig, axes = small_multiples(len(chosen))
    for ax, run in zip(axes.flat, chosen, strict=False):
        for line in grid_lines(run, stage=1):
            ax.plot(line[:, 1], line[:, 0], color=RAMP(T2_SCALE(line[:, 2].mean())), marker="o", markersize=3)
        ax.set_title(f"{run.label} ({run.compressors or '?'} comp.)", fontsize=8, loc="left")
    for ax in axes[-1]:
        ax.set_xlabel("T1 (K)")
    for ax in axes[:, 0]:
        ax.set_ylabel("Q1 (W)")
    colorbar(fig, list(axes.flat), T2_SCALE, "mean T2 along the line (K)")
    fig.suptitle(
        "Q1 against T1, one line per stage-2 power; T2 varies by at most 0.6 K along a line, so each line "
        "is close to Q1(T1) at fixed T2",
        x=0.01,
        ha="left",
    )
    save(fig, "13_q1_against_t1.png")


def plot_stage1_against_t2(runs: list[Run]) -> None:
    """Q1 at fixed T1 against T2, colored by T1, for grids with at least three stage-2 powers."""
    chosen = [r for r in grids_of(runs) if len(levels(r.q2)) >= 3 and len(levels(r.q1)) >= 2]
    fig, axes = small_multiples(len(chosen), columns=6)
    for ax, run in zip(axes.flat, chosen, strict=False):
        for target, points in cuts(run, stage=1, spacing=2.5):
            ax.plot(points[:, 0], points[:, 1], color=RAMP(T1_SCALE(target)), marker="o", markersize=2.5)
        ax.set_title(f"{run.label} ({run.compressors or '?'} comp.)", fontsize=8, loc="left")
    for ax in axes[-1]:
        ax.set_xlabel("T2 (K)")
    for ax in axes[:, 0]:
        ax.set_ylabel("Q1 (W)")
    colorbar(fig, list(axes.flat), T1_SCALE, "T1 of the cut (K)")
    fig.suptitle(
        "Q1 at fixed T1 against T2: cuts every 2.5 K of T1, interpolated along each line of fixed Q2", x=0.01, ha="left"
    )
    save(fig, "14_q1_against_t2.png")


def plot_stage1_summary(runs: list[Run], notes: list[str]) -> None:
    """dQ1/dT1 against T1 in log-log, Q1 at fixed T1 against T2 for every grid, and dQ1/dT2 against T1."""
    fig, axes = plt.subplots(1, 3, figsize=(20, 5), gridspec_kw={"wspace": 0.3})
    plot_slope_fits(axes[0], runs, "t1", "dq1_dt1", span=8.0, notes=notes)
    axes[0].set(title="dQ1/dT1 at fixed T2, log-log", xlabel="T1 (K)", ylabel="dQ1/dT1 (W/K)")
    axes[0].set_xticks([40, 50, 60], ["40", "50", "60"])
    for run in grids_of(runs):
        for target, points in cuts(run, stage=1, spacing=2.5):
            axes[1].plot(
                points[:, 0],
                points[:, 1] - points[0, 1],
                color=RAMP(T1_SCALE(target)),
                linewidth=0.9,
                alpha=0.8,
                linestyle="-" if run.compressors == 1 else "--",
            )
    axes[1].set(
        title="Q1 at fixed T1 against T2, relative to the cut's lowest T2\n(solid one compressor, dashed two)",
        xlabel="T2 (K)",
        ylabel="Q1 − Q1 at the cut's lowest T2 (W)",
    )
    colorbar(fig, axes[1], T1_SCALE, "T1 of the cut (K)")
    for name, group in compressor_groups(runs).items():
        color = GROUP_COLOR[name]
        cells = [c for run in group for c in cell_partials(run)]
        t1, t2, y = (np.array([c[k] for c in cells]) for k in ("t1", "t2", "dq1_dt2"))
        axes[2].scatter(t1, y, s=12, color=color, alpha=0.35, linewidths=0)
        design = np.column_stack([np.ones_like(t1), t1 - 50, t2 - 3.5])
        coefficients, *_ = np.linalg.lstsq(design, y, rcond=None)
        line = np.linspace(t1.min(), t1.max(), 20)
        axes[2].plot(
            line,
            coefficients[0] + coefficients[1] * (line - 50),
            color=color,
            linewidth=2.2,
            label=f"{name}, pooled at T2 = 3.5 K",
        )
        own = []
        for run in group:
            fitted = run_cells(run, "t1", "dq1_dt2", span=8.0)
            if fitted:
                slope, intercept = np.polyfit(*fitted, 1)
                line = np.linspace(fitted[0].min(), fitted[0].max(), 10)
                axes[2].plot(line, intercept + slope * line, color=color, linewidth=0.9, linestyle="--")
                own.append(f"{run.label} {slope:+.3f}")
        notes.append(
            f"dQ1/dT2, {name}: {coefficients[0]:.2f} W/K {coefficients[1]:+.3f} W/K^2 (T1 - 50 K) "
            f"{coefficients[2]:+.3f} W/K^2 (T2 - 3.5 K), rms {np.std(y - design @ coefficients):.2f} W/K over "
            f"{y.size} cells; per-run slope against T1 in W/K^2: " + ("; ".join(own) or "none span enough")
        )
        tilt = np.array([-c["dq1_dt2"] / c["dq1_dt1"] for c in cells])
        notes.append(
            f"dT1/dT2 at fixed Q1, {name}: median {np.median(tilt):.2f}, quartiles {np.percentile(tilt, 25):.2f} to "
            f"{np.percentile(tilt, 75):.2f}"
        )
    axes[2].plot([], [], color=MUTED, linestyle="--", linewidth=0.9, label="single runs")
    axes[2].axhline(0, color=MUTED, linewidth=0.8)
    axes[2].set(title="dQ1/dT2 at fixed T1 against T1", xlabel="T1 (K)", ylabel="dQ1/dT2 (W/K)", ylim=(-5, 1))
    axes[2].legend(loc="lower right")
    save(fig, "15_q1_summary.png")


def main() -> None:
    """Draws the functional-form figures and writes their fitted numbers."""
    runs = load_runs()
    notes: list[str] = []
    plot_stage2(runs, notes)
    plot_stage1_against_t1(runs)
    plot_stage1_against_t2(runs)
    plot_stage1_summary(runs, notes)
    (HERE / "functional_forms.txt").write_text("\n".join(notes) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
