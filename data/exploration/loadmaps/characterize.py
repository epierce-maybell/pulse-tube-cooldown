"""Draws the load-map characterization figures into `figures/` and writes the numbers they show to `summary.txt`.

Every run in `runs.py` is its own dataset. Temperatures at a repeated pair of heater powers are averaged except in
the repeatability comparison. Each figure's docstring says what it shows.
"""

import csv
import itertools
from collections.abc import Callable, Iterable
from datetime import date
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from numpy.typing import NDArray
from runs import Run, load_runs

HERE = Path(__file__).parent
STAND_COLOR = {"Münster": "#eb6834", "Denver": "#2a78d6"}
VARIANT_MARKER = {"regular 1.5": "o", "regular 1.5 HX1 prototype": "D", "IBM": "s", "unknown": "^"}
MUTED = "#8a8985"
plt.rcParams.update(
    {
        "axes.grid": True,
        "grid.color": "#e4e3df",
        "grid.linewidth": 0.6,
        "axes.edgecolor": "#8a8985",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "font.size": 9,
        "axes.titlesize": 9,
        "legend.fontsize": 8,
        "legend.frameon": False,
        "lines.linewidth": 1.2,
        "savefig.dpi": 150,
    }
)
STANDARD_POINTS = [(0.0, 0.0), (0.0, 0.5), (20.0, 0.0), (20.0, 0.5)]
TIERS = [
    "repeat within a run",
    "same head and stand,\n1 day apart",
    "same head and stand,\n8 to 103 days apart",
    "same head,\nother stand",
    "different heads,\nDenver",
    "different heads,\nMünster",
]


def value(run: Run, q1: float, q2: float, quantity: str) -> float:
    """Returns a run's mean temperature at a pair of heater powers, or NaN if the run lacks that pair.

    Args:
        run: The run.
        q1: Stage-1 heater power in W.
        q2: Stage-2 heater power in W.
        quantity: "t1", "t2", or "t3".
    """
    values = getattr(run, quantity)[(run.q1 == q1) & (run.q2 == q2)]
    return float(np.nanmean(values)) if values.size and not np.isnan(values).all() else np.nan


def levels(values: NDArray[np.float64]) -> list[float]:
    """Returns the distinct heater powers in increasing order."""
    return sorted(set(values.tolist()))


def grids_of(runs: Iterable[Run]) -> list[Run]:
    """Returns the runs that are maps over pairs of heater powers."""
    return [r for r in runs if r.kind == "grid"]


def run_style(run: Run, runs: list[Run]) -> dict:
    """Returns line and marker styles that tell a head's runs on one stand apart, in date order."""
    same = sorted(r.date for r in grids_of(runs) if r.head == run.head and r.stand == run.stand)
    color = STAND_COLOR[run.stand]
    return {
        "color": color,
        "linestyle": ["-", "--", ":", "-.", (0, (5, 1, 1, 1))][same.index(run.date)],
        "marker": VARIANT_MARKER[run.variant],
        "markersize": 4,
        "markeredgecolor": color,
        "markerfacecolor": color if run.compressors == 1 else "white",
    }


def point_style(run: Run, runs: list[Run], size: float = 4.5) -> dict:
    """Returns `run_style` for markers drawn without lines."""
    style = run_style(run, runs)
    return {**style, "linestyle": "", "markersize": size}


def style_legend(target: plt.Axes | Figure, include_runs: bool = False, **kwargs) -> None:
    """Adds the shared legend for stand color, compressor fill, head variant, and run order."""
    handles = [Line2D([], [], color=c, marker="o", label=s) for s, c in STAND_COLOR.items()]
    handles += [
        Line2D([], [], color=MUTED, marker="o", linestyle="", label="one compressor"),
        Line2D([], [], color=MUTED, marker="o", markerfacecolor="white", linestyle="", label="two, or unrecorded"),
    ]
    handles += [Line2D([], [], color=MUTED, marker=m, linestyle="", label=v) for v, m in VARIANT_MARKER.items()]
    if include_runs:
        handles += [
            Line2D([], [], color=MUTED, linestyle=d, label=label)
            for d, label in zip(
                ["-", "--", ":", "-."], ["first run on a stand", "second", "third", "fourth"], strict=True
            )
        ]
    target.legend(handles=handles, **kwargs)


def figure_legend(fig: Figure, include_runs: bool = False) -> None:
    """Places the shared legend to the right of a figure's panels."""
    style_legend(fig, include_runs, loc="upper left", bbox_to_anchor=(1.0, 0.92))


def save(fig: Figure, name: str) -> None:
    """Writes a figure to `figures/` and closes it."""
    fig.savefig(HERE / "figures" / name, bbox_inches="tight")
    plt.close(fig)


def heads_in_order(runs: Iterable[Run]) -> list[str]:
    """Returns the heads by serial, with heads the overview lacks last."""
    return sorted({r.head for r in runs}, key=lambda h: (not h.startswith("S0"), h))


def plot_grids(runs: list[Run]) -> None:
    """Every grid in the T1-T2 plane, one panel per head, drawn as the manufacturer draws a map."""
    grids = grids_of(runs)
    heads = heads_in_order(grids)
    columns = 6
    fig, axes = plt.subplots(int(np.ceil(len(heads) / columns)), columns, figsize=(19, 12.5), sharex=True, sharey=True)
    for ax, head in zip(axes.flat, heads, strict=False):
        for run in (r for r in grids if r.head == head):
            style = run_style(run, runs)
            for q1 in levels(run.q1):
                q2s = [q2 for q2 in levels(run.q2) if not np.isnan(value(run, q1, q2, "t1"))]
                ax.plot([value(run, q1, q2, "t1") for q2 in q2s], [value(run, q1, q2, "t2") for q2 in q2s], **style)
            for q2 in levels(run.q2):
                q1s = [q1 for q1 in levels(run.q1) if not np.isnan(value(run, q1, q2, "t1"))]
                ax.plot(
                    [value(run, q1, q2, "t1") for q1 in q1s],
                    [value(run, q1, q2, "t2") for q1 in q1s],
                    color=style["color"],
                    linestyle=style["linestyle"],
                    linewidth=0.5,
                    alpha=0.7,
                )
        ax.set_title(f"{head}  ({next(r.variant for r in grids if r.head == head)})", loc="left")
    for ax in axes.flat[len(heads) :]:
        ax.set_axis_off()
    for ax in axes[-1]:
        ax.set_xlabel("T1 (K)")
    for ax in axes[:, 0]:
        ax.set_ylabel("T2 (K)")
    style_legend(axes.flat[-1], include_runs=True, loc="center")
    fig.suptitle(
        "Load maps by head: bold lines hold stage-1 power fixed, thin lines hold stage-2 power fixed",
        x=0.01,
        y=0.93,
        ha="left",
    )
    save(fig, "01_grids_by_head.png")


def plot_standard_points(runs: list[Run]) -> None:
    """T1 and T2 at the standard points of every run, one column per head."""
    grids = grids_of(runs)
    heads = heads_in_order(grids)
    fig, axes = plt.subplots(len(STANDARD_POINTS), 2, figsize=(17, 10), sharex=True)
    for row, (q1, q2) in enumerate(STANDARD_POINTS):
        for column, quantity in enumerate(("t1", "t2")):
            ax = axes[row, column]
            for run in grids:
                same = sorted(r.date for r in grids if r.head == run.head and r.stand == run.stand)
                offset = (0.1 + 0.08 * same.index(run.date)) * (1 if run.stand == "Denver" else -1)
                ax.plot(heads.index(run.head) + offset, value(run, q1, q2, quantity), **point_style(run, runs, 5))
            ax.set_ylabel(f"{quantity.upper()} (K) at {q1:g} W, {q2:g} W")
    for ax in axes[-1]:
        ax.set_xticks(range(len(heads)), heads, rotation=60, ha="right")
    figure_legend(fig)
    fig.suptitle(
        "Standard points of every run, labeled (stage-1 power, stage-2 power); Münster runs left of each "
        "head, Denver runs right, in date order",
        x=0.01,
        y=0.93,
        ha="left",
    )
    save(fig, "02_standard_points.png")


def days_between(a: Run, b: Run) -> int:
    """Returns the number of days between two runs."""
    return abs((date.fromisoformat(a.date) - date.fromisoformat(b.date)).days)


def tier(a: Run, b: Run) -> str | None:
    """Classifies a pair of runs, or returns None for pairs outside the comparison.

    Different heads are compared only within one stand, among regular 1.5 heads on one compressor.
    """
    if a.head == b.head and a.head != "unknown":
        if a.stand != b.stand:
            return TIERS[3]
        return TIERS[1] if days_between(a, b) <= 1 else TIERS[2]
    single = all(r.variant == "regular 1.5" and r.compressors == 1 for r in (a, b))
    if single and a.stand == b.stand:
        return TIERS[4] if a.stand == "Denver" else TIERS[5]
    return None


def repeat_differences(run: Run, quantity: str) -> list[float]:
    """Returns the spread of each pair of heater powers a run measured more than once."""
    out = []
    for q1, q2 in {(a, b) for a, b in zip(run.q1, run.q2, strict=True)}:
        values = getattr(run, quantity)[(run.q1 == q1) & (run.q2 == q2)]
        if values.size > 1 and not np.isnan(values).any():
            out.append(float(values.max() - values.min()))
    return out


def tier_differences(runs: list[Run], q1: float, q2: float, quantity: str) -> dict[str, NDArray[np.float64]]:
    """Returns the absolute differences between pairs of datasets at one point, by tier.

    The first tier holds every repeated point within a run, whatever its powers.
    """
    grids = grids_of(runs)
    by_tier = {t: [] for t in TIERS}
    for run in grids:
        by_tier[TIERS[0]] += repeat_differences(run, quantity)
    for a, b in itertools.combinations(grids, 2):
        kind = tier(a, b)
        difference = abs(value(a, q1, q2, quantity) - value(b, q1, q2, quantity))
        if kind and not np.isnan(difference):
            by_tier[kind].append(difference)
    return {t: np.array(d) for t, d in by_tier.items()}


def plot_variability(runs: list[Run]) -> None:
    """Differences between datasets at the standard points, by how closely related the two datasets are."""
    fig, axes = plt.subplots(2, len(STANDARD_POINTS), figsize=(17, 8), sharey=True, sharex=True)
    rng = np.random.default_rng(0)
    for column, (q1, q2) in enumerate(STANDARD_POINTS):
        for row, quantity in enumerate(("t1", "t2")):
            ax = axes[row, column]
            for i, d in enumerate(tier_differences(runs, q1, q2, quantity).values()):
                if d.size:
                    ax.plot(
                        i + rng.uniform(-0.18, 0.18, d.size),
                        d,
                        "o",
                        markersize=3.5,
                        color="#2a78d6",
                        alpha=0.55,
                        markeredgewidth=0,
                    )
                    ax.plot([i - 0.28, i + 0.28], [np.median(d)] * 2, color="#0b0b0b", linewidth=1.5)
            ax.set_yscale("log")
            ax.set_ylim(1e-3, 30)
            ax.set_title(f"|ΔT{quantity[1]}| at {q1:g} W, {q2:g} W")
    for ax in axes[-1]:
        ax.set_xticks(range(len(TIERS)), TIERS, rotation=60, ha="right")
    for ax in axes[:, 0]:
        ax.set_ylabel("absolute difference (K)")
    fig.suptitle(
        "Differences between pairs of datasets, grouped by how related they are; bars are medians. The first "
        "group is every repeated point within a run. Different heads: regular 1.5 on one compressor.",
        x=0.01,
        ha="left",
    )
    save(fig, "03_variability_tiers.png")


def secants(x: list[float], y: list[float]) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Returns the midpoints and slopes of the segments joining successive points."""
    x, y = np.asarray(x), np.asarray(y)
    return (x[1:] + x[:-1]) / 2, np.diff(y) / np.diff(x)


def own_curve(run: Run, stage: int, other_power: float) -> tuple[list[float], list[float]]:
    """Returns a stage's own heater powers and temperatures with the other stage held at one power."""
    own_powers = levels(run.q1 if stage == 1 else run.q2)
    pairs = [(p, other_power) if stage == 1 else (other_power, p) for p in own_powers]
    temps = [value(run, *pair, f"t{stage}") for pair in pairs]
    kept = [(p, t) for p, t in zip(own_powers, temps, strict=True) if not np.isnan(t)]
    return [p for p, _ in kept], [t for _, t in kept]


def plot_own_response(runs: list[Run], stage: int) -> None:
    """A stage's temperature against its own heater power with the other stage at its lowest power.

    Panels: absolute, shifted to the run's own unloaded value, the secant slope against power, and the mean slope with
    the other stage at its lowest power against the mean slope at its highest.
    """
    other = "q2" if stage == 1 else "q1"
    fig, axes = plt.subplots(1, 4, figsize=(19, 4.6))
    sweep_style = {"color": "#1baf7a", "marker": "o", "markersize": 4}
    for run in runs:
        if run.kind == "sweep" and stage == 1:
            continue
        style = run_style(run, runs) if run.kind == "grid" else sweep_style
        powers, temps = own_curve(run, stage, levels(getattr(run, other))[0])
        if len(powers) < 2:
            continue
        axes[0].plot(powers, temps, **style)
        axes[1].plot(powers, np.array(temps) - temps[0], **style)
        axes[2].plot(*secants(powers, temps), **style)
        if run.kind == "grid" and len(levels(getattr(run, other))) > 1:
            top_powers, top_temps = own_curve(run, stage, levels(getattr(run, other))[-1])
            if top_powers == powers:
                span = powers[-1] - powers[0]
                axes[3].plot(
                    (temps[-1] - temps[0]) / span, (top_temps[-1] - top_temps[0]) / span, **point_style(run, runs)
                )
    own_name, other_name = f"Q{stage}", f"Q{3 - stage}"
    axes[0].set(
        xlabel=f"{own_name} (W)", ylabel=f"T{stage} (K)", title=f"T{stage} at the lowest {other_name} of each run"
    )
    axes[1].set(
        xlabel=f"{own_name} (W)",
        ylabel=f"T{stage} − T{stage}({own_name} = 0) (K)",
        title="shifted to the run's unloaded value",
    )
    axes[2].set(
        xlabel=f"{own_name} at segment midpoint (W)",
        ylabel=f"ΔT{stage} / Δ{own_name} (K/W)",
        title="secant slope between successive powers",
    )
    axes[3].set(
        xlabel=f"mean slope at lowest {other_name} (K/W)",
        ylabel=f"mean slope at highest {other_name} (K/W)",
        title=f"does the slope depend on {other_name}?",
    )
    low, high = min(axes[3].get_xlim()[0], axes[3].get_ylim()[0]), max(axes[3].get_xlim()[1], axes[3].get_ylim()[1])
    axes[3].plot([low, high], [low, high], color=MUTED, linewidth=0.8, linestyle="--")
    figure_legend(fig)
    if stage == 2:
        fig.legend(
            handles=[Line2D([], [], **sweep_style, label="S0028 Münster sweep\nat 20 W, two compressors")],
            loc="upper left",
            bbox_to_anchor=(1.0, 0.25),
        )
    save(fig, f"0{3 + stage}_stage{stage}_response.png")


def plot_cross_response(runs: list[Run]) -> None:
    """Each stage's temperature against the other stage's heater power, at the lowest and highest own power."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 9))
    for run in grids_of(runs):
        style = run_style(run, runs)
        q1s, q2s = levels(run.q1), levels(run.q2)
        for column, q1 in enumerate((q1s[0], q1s[-1])):
            temps = np.array([value(run, q1, q2, "t1") for q2 in q2s])
            if len(q2s) > 1 and not np.isnan(temps[0]):
                axes[0, column].plot(q2s, temps - temps[0], **style)
        for column, q2 in enumerate((q2s[0], q2s[-1])):
            temps = np.array([value(run, q1, q2, "t2") for q1 in q1s])
            if len(q1s) > 1:
                axes[1, column].plot(q1s, temps - temps[0], **style)
    axes[0, 0].set(title="T1 against stage-2 power, stage 1 unloaded", xlabel="Q2 (W)", ylabel="T1 − T1(Q2 = 0) (K)")
    axes[0, 1].set(
        title="T1 against stage-2 power, stage 1 at the run's highest power (20 to 35 W)",
        xlabel="Q2 (W)",
        ylabel="T1 − T1(Q2 = 0) (K)",
    )
    axes[1, 0].set(title="T2 against stage-1 power, stage 2 unloaded", xlabel="Q1 (W)", ylabel="T2 − T2(Q1 = 0) (K)")
    axes[1, 1].set(
        title="T2 against stage-1 power, stage 2 at the run's highest power (0.5 to 2.5 W)",
        xlabel="Q1 (W)",
        ylabel="T2 − T2(Q1 = 0) (K)",
    )
    figure_legend(fig, include_runs=True)
    save(fig, "06_cross_response.png")


def cell_partials(run: Run) -> list[dict]:
    """Returns the local partial derivatives of heater power with respect to temperature in each grid cell.

    A cell is the rectangle of four measured points at adjacent levels of both powers. The partials of temperature
    with respect to power are differences along the cell's edges, averaged over its two parallel edges, and the
    partials of power with respect to temperature are their matrix inverse.
    """
    cells = []
    for (a1, b1), (a2, b2) in itertools.product(itertools.pairwise(levels(run.q1)), itertools.pairwise(levels(run.q2))):
        corner = {
            (p1, p2): np.array([value(run, p1, p2, "t1"), value(run, p1, p2, "t2")])
            for p1 in (a1, b1)
            for p2 in (a2, b2)
        }
        if np.isnan(np.array(list(corner.values()))).any():
            continue
        along_q1 = (corner[(b1, a2)] - corner[(a1, a2)] + corner[(b1, b2)] - corner[(a1, b2)]) / (2 * (b1 - a1))
        along_q2 = (corner[(a1, b2)] - corner[(a1, a2)] + corner[(b1, b2)] - corner[(b1, a2)]) / (2 * (b2 - a2))
        inverse = np.linalg.inv(np.column_stack([along_q1, along_q2]))
        center = np.mean(list(corner.values()), axis=0)
        cells.append(
            {
                "run": run,
                "t1": center[0],
                "t2": center[1],
                "dq1_dt1": inverse[0, 0],
                "dq1_dt2": inverse[0, 1],
                "dq2_dt1": inverse[1, 0],
                "dq2_dt2": inverse[1, 1],
                "dq1_dq2": -along_q2[0] / along_q1[0],
            }
        )
    return cells


PARTIALS = [
    ("dq1_dt1", "∂Q1/∂T1 at fixed T2 (W/K)", "t1", (0, 2.6)),
    ("dq1_dt2", "∂Q1/∂T2 at fixed T1 (W/K)", "t2", (-5, 1)),
    ("dq2_dt1", "∂Q2/∂T1 at fixed T2 (W/K)", "t1", (-0.05, 0.05)),
    ("dq2_dt2", "∂Q2/∂T2 at fixed T1 (W/K)", "t2", (0, 1.8)),
]


def plot_capacity_partials(runs: list[Run]) -> None:
    """The local partials of heater power with respect to stage temperature, from every grid cell.

    This inverts each cell's measured response; it is the map read as Q1(T1, T2) and Q2(T1, T2).
    """
    cells = [c for run in grids_of(runs) for c in cell_partials(run)]
    fig, axes = plt.subplots(2, 2, figsize=(14, 9))
    outside = 0
    for ax, (key, label, x_key, limits) in zip(axes.flat, PARTIALS, strict=True):
        for cell in cells:
            ax.plot(cell[x_key], cell[key], alpha=0.8, **point_style(cell["run"], runs, 4))
        outside += sum(not limits[0] <= c[key] <= limits[1] for c in cells)
        ax.set(ylabel=label, xlabel=f"cell mean {x_key.upper()} (K)", ylim=limits)
        ax.axhline(0, color=MUTED, linewidth=0.8)
    figure_legend(fig)
    fig.suptitle(
        f"The map read as heater power against stage temperatures: local partials from {len(cells)} grid "
        f"cells; {outside} of {4 * len(cells)} values fall outside the plotted ranges",
        x=0.01,
        ha="left",
    )
    save(fig, "07_capacity_partials.png")


def plot_capacity_ratio(runs: list[Run]) -> None:
    """Plots (dQ1/dT2) / (dQ2/dT2) at fixed T1 against each cell's mean stage temperatures.

    The dimensionless ratio is dQ1/dQ2 at fixed T1: the change in stage-1 heater power per added watt on stage 2.
    Each point uses the partials from one measured grid cell.
    """
    cells = [cell for run in grids_of(runs) for cell in cell_partials(run)]
    fig, axes = plt.subplots(1, 2, figsize=(14, 4.6), sharey=True)
    for ax, temperature in zip(axes, ("t1", "t2"), strict=True):
        for cell in cells:
            ratio = cell["dq1_dt2"] / cell["dq2_dt2"]
            ax.plot(cell[temperature], ratio, alpha=0.8, **point_style(cell["run"], runs, 4))
        ax.set(xlabel=f"cell mean {temperature.upper()} (K)")
        ax.axhline(0, color=MUTED, linewidth=0.8)
    axes[0].set_ylabel("(dQ1/dT2) / (dQ2/dT2) at fixed T1 (W/W)")
    figure_legend(fig)
    fig.suptitle(f"Stage-1 power change per watt on stage 2 at fixed T1: {len(cells)} grid cells")
    save(fig, "16_capacity_ratio.png")


def collapse(runs: list[Run], stage: int, run_slope: bool = False, t2_terms: bool = False) -> dict:
    """Fits every run's heater power on one stage to a shared shape plus a constant for each run.

    The power on stage s is a shared quadratic in T_s, minus a constant for each run. With `run_slope`, the linear
    term is the run's own; with `t2_terms`, a shared quadratic in T2 is added, for stage 1. Temperatures are taken
    relative to 40 K on stage 1 and 2 K on stage 2.

    Args:
        runs: Grids to fit together.
        stage: 1 or 2.
        run_slope: Whether each run has its own linear term.
        t2_terms: Whether stage-1 power also depends on T2.

    Returns:
        The points' own-stage temperatures, powers, run indices, residuals, each run's constant, and the shared
        coefficients.
    """
    reference = 40.0 if stage == 1 else 2.0
    x, t2, power, index = [], [], [], []
    for i, run in enumerate(runs):
        kept = ~np.isnan(run.t1) & ~np.isnan(run.t2)
        x.append(getattr(run, f"t{stage}")[kept] - reference)
        t2.append(run.t2[kept] - 2.0)
        power.append(getattr(run, f"q{stage}")[kept])
        index.append(np.full(kept.sum(), i))
    x, t2, power, index = map(np.concatenate, (x, t2, power, index))
    shared = [x**2] if run_slope else [x, x**2]
    if t2_terms:
        shared += [t2, t2**2]
    membership = (index[:, None] == np.arange(len(runs))).astype(float)
    design = np.column_stack(shared + list(-membership.T) + (list((membership * x[:, None]).T) if run_slope else []))
    coefficients, *_ = np.linalg.lstsq(design, power, rcond=None)
    return {
        "temperature": x + reference,
        "power": power,
        "run": index,
        "residual": power - design @ coefficients,
        "constant": coefficients[len(shared) : len(shared) + len(runs)],
        "shared": coefficients[: len(shared)],
    }


def collapse_groups(runs: list[Run]) -> dict[str, list[Run]]:
    """Returns the groups of grids the collapse tests fit together."""
    grids = grids_of(runs)
    return {
        "regular 1.5 and HX1, one compressor": [
            r for r in grids if r.compressors == 1 and r.variant in ("regular 1.5", "regular 1.5 HX1 prototype")
        ],
        "every variant, one compressor": [r for r in grids if r.compressors == 1],
        "every variant, two compressors": [r for r in grids if r.compressors == 2],
    }


COLLAPSE_TESTS = [
    ("stage 2: shared quadratic in T2, constant per run", 2, False, False),
    ("stage 2: as above, slope per run", 2, True, False),
    ("stage 1: shared quadratic in T1, constant per run", 1, False, False),
    ("stage 1: as above, slope per run", 1, True, False),
    ("stage 1: shared quadratic in T1, shared quadratic in T2, constant per run", 1, False, True),
    ("stage 1: as above, slope per run", 1, True, True),
]


def plot_collapse(runs: list[Run]) -> None:
    """Stage-2 power against T2 for every run, raw and collapsed onto one shared curve by a constant per run.

    Rows are one compressor (regular 1.5 and HX1) and two compressors (every variant). The last panel is the
    residual of the best stage-1 collapse in the summary, for comparison.
    """
    groups = collapse_groups(runs)
    fig, axes = plt.subplots(2, 4, figsize=(19, 9), gridspec_kw={"hspace": 0.4})
    for row, name in enumerate(["regular 1.5 and HX1, one compressor", "every variant, two compressors"]):
        axes[row, 0].annotate(name, (-0.32, 0.5), xycoords="axes fraction", rotation=90, va="center", fontweight="bold")
        group = groups[name]
        fit = collapse(group, stage=2)
        for i, run in enumerate(group):
            mine = fit["run"] == i
            style = point_style(run, runs, 4)
            axes[row, 0].plot(fit["temperature"][mine], fit["power"][mine], **style)
            axes[row, 1].plot(fit["temperature"][mine], fit["power"][mine] + fit["constant"][i], **style)
            axes[row, 2].plot(fit["temperature"][mine], fit["residual"][mine], **style)
        grid = np.linspace(fit["temperature"].min(), fit["temperature"].max(), 100) - 2.0
        axes[row, 1].plot(
            grid + 2.0, fit["shared"][0] * grid + fit["shared"][1] * grid**2, color="#0b0b0b", linewidth=1.0
        )
        rms = np.sqrt(np.mean(fit["residual"] ** 2))
        axes[row, 0].set(title="Q2 against T2, every point", xlabel="T2 (K)", ylabel="Q2 (W)")
        axes[row, 1].set(
            title="shifted by each run's constant; black: shared curve", xlabel="T2 (K)", ylabel="Q2 + run constant (W)"
        )
        axes[row, 2].set(
            title=f"residual, rms {rms:.3f} W over {fit['power'].size} points", xlabel="T2 (K)", ylabel="residual (W)"
        )
        best = collapse(group, stage=1, run_slope=True, t2_terms=True)
        for i, run in enumerate(group):
            mine = best["run"] == i
            axes[row, 3].plot(best["temperature"][mine], best["residual"][mine], **point_style(run, runs, 4))
        axes[row, 3].set(
            title=f"stage 1, slope and constant per run, shared T2 terms:\nresidual rms "
            f"{np.sqrt(np.mean(best['residual'] ** 2)):.2f} W",
            xlabel="T1 (K)",
            ylabel="Q1 residual (W)",
        )
    figure_legend(fig)
    save(fig, "11_shape_collapse.png")


def plot_pulse_tube(runs: list[Run]) -> None:
    """T3, the Denver stand's pulse-tube thermometer, against stage-1 power and against T1."""
    denver = [r for r in runs if r.stand == "Denver"]
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    labels = []
    for run in denver:
        excursion = np.nanmax(run.t3) > 37 and run.compressors == 1
        color, alpha = ("#e34948", 0.9) if excursion else ("#2a78d6", 0.5)
        q2s = levels(run.q2)
        for q2, dash in ((q2s[0], "-"), (q2s[-1], "--")):
            powers = [p for p in levels(run.q1) if not np.isnan(value(run, p, q2, "t3"))]
            axes[0].plot(
                powers,
                [value(run, p, q2, "t3") for p in powers],
                color=color,
                linestyle=dash,
                marker="o",
                markersize=3,
                alpha=alpha,
            )
        axes[1].plot(
            run.t1,
            run.t3,
            linestyle="",
            marker="o" if run.compressors == 1 else "s",
            markersize=4,
            color=color,
            markerfacecolor=color if run.compressors == 1 else "white",
            alpha=alpha,
        )
        if excursion:
            labels.append((float(np.nanmax(run.t3)), run.label.replace(" Denver", "")))
    position = np.inf
    for t3, text in sorted(labels, reverse=True):
        position = min(t3, position - 1.6)
        axes[0].annotate(text, (31, position), fontsize=7, color="#52514e", va="center")
    axes[0].set(
        xlabel="Q1 (W)",
        ylabel="T3 (K)",
        xlim=(-2, 44),
        xticks=range(0, 35, 5),
        title="T3 against stage-1 power; solid lowest Q2, dashed highest; red runs exceed 37 K",
    )
    axes[1].set(xlabel="T1 (K)", ylabel="T3 (K)", title="T3 against T1, every Denver point")
    axes[1].plot([30, 70], [30, 70], color=MUTED, linewidth=0.8, linestyle="--")
    axes[1].legend(
        handles=[
            Line2D([], [], color="#e34948", marker="o", linestyle="", label="one compressor, T3 above 37 K"),
            Line2D([], [], color="#2a78d6", marker="o", linestyle="", label="one compressor, other runs"),
            Line2D(
                [],
                [],
                color="#2a78d6",
                marker="s",
                markerfacecolor="white",
                linestyle="",
                label="two or unrecorded (S0028, S0009 Aug)",
            ),
            Line2D([], [], color=MUTED, linestyle="--", label="T3 = T1"),
        ],
        loc="upper left",
    )
    save(fig, "08_pulse_tube_t3.png")


def plot_valve_settings(runs: list[Run]) -> None:
    """Münster map features of regular 1.5 heads against each shipped valve setting."""
    with open(HERE / "valve_settings.csv", encoding="utf-8") as file:
        settings = {row["head"]: row for row in csv.DictReader(file)}
    muenster = [r for r in grids_of(runs) if r.stand == "Münster" and r.variant == "regular 1.5"]
    features: list[tuple[str, Callable[[Run], float]]] = [
        ("T1 at 0 W, 0 W (K)", lambda r: value(r, 0, 0, "t1")),
        ("T2 at 0 W, 0 W (K)", lambda r: value(r, 0, 0, "t2")),
        ("T1 rise from 0 to 0.5 W\non stage 2, Q1 = 0 (K)", lambda r: value(r, 0, 0.5, "t1") - value(r, 0, 0, "t1")),
        ("T2 rise from 0 to 0.5 W\non stage 2, Q1 = 0 (K)", lambda r: value(r, 0, 0.5, "t2") - value(r, 0, 0, "t2")),
        ("T1 rise from 0 to 20 W\non stage 1, Q2 = 0 (K)", lambda r: value(r, 20, 0, "t1") - value(r, 0, 0, "t1")),
    ]
    valves = ["ORV1", "DIV1A", "DIV1B", "ORV2", "DIV2A", "DIV2B", "valve_frequency"]
    fig, axes = plt.subplots(len(features), len(valves), figsize=(19, 12), sharex="col", sharey="row")
    for run in muenster:
        row = settings[run.head]
        for (i, (_, feature)), (j, valve) in itertools.product(enumerate(features), enumerate(valves)):
            if row[valve]:
                axes[i, j].plot(float(row[valve]), feature(run), **point_style(run, runs, 5))
                axes[i, j].annotate(
                    run.head[-2:],
                    (float(row[valve]), feature(run)),
                    fontsize=6,
                    xytext=(3, 2),
                    textcoords="offset points",
                    color="#52514e",
                )
    for i, (label, _) in enumerate(features):
        axes[i, 0].set_ylabel(label)
    for j, valve in enumerate(valves):
        axes[-1, j].set_xlabel("valve frequency (Hz)" if valve == "valve_frequency" else f"{valve} (turns open)")
    fig.suptitle(
        "Münster maps of regular 1.5 heads against shipped valve settings; open markers two compressors; "
        "labels are the last two digits of the serial",
        x=0.01,
        ha="left",
    )
    save(fig, "09_valve_settings.png")


def plot_digitizing_check(runs: list[Run]) -> None:
    """Digitized Münster values against S0007's numeric Münster file and against the overview's key points."""
    with open(HERE / "muenster_maps.csv", encoding="utf-8") as file:
        digitized = [r for r in csv.DictReader(file) if r["kind"] == "grid"]
    with open(HERE / "muenster_key_points.csv", encoding="utf-8") as file:
        keys = {(k["head"], float(k["Q1"]), float(k["Q2"])): k for k in csv.DictReader(file)}
    numeric = next(r for r in runs if r.head == "S0007" and r.stand == "Münster")
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.5))
    for row in digitized:
        q1, q2, t1, t2 = (float(row[k]) for k in ("Q1", "Q2", "T1", "T2"))
        if row["head"] == "S0007":
            axes[0].plot(t1 - value(numeric, q1, q2, "t1"), t2 - value(numeric, q1, q2, "t2"), "o", color="#2a78d6")
        key = keys.get((row["head"], q1, q2))
        if key:
            d1, d2 = t1 - float(key["T1"]), t2 - float(key["T2"])
            axes[1].plot(d1, d2, "o", color="#eb6834")
            if abs(d2) > 0.05:
                axes[1].annotate(
                    f"{row['head']} {q1:g}/{q2:g} W", (d1, d2), fontsize=7, xytext=(4, 2), textcoords="offset points"
                )
    axes[0].set(
        xlabel="digitized − numeric T1 (K)",
        ylabel="digitized − numeric T2 (K)",
        title="S0007: digitized figure against the numeric Münster file (15 points)",
    )
    axes[1].set(
        xlabel="digitized − key point T1 (K)",
        ylabel="digitized − key point T2 (K)",
        title="every head: digitized figure against the overview's printed key points",
    )
    save(fig, "10_digitizing_check.png")


def summarize(runs: list[Run]) -> str:
    """Returns the numbers the note cites, as text."""
    lines = ["Differences between datasets: median |dT| in K (count of pairs)"]
    for q1, q2 in STANDARD_POINTS:
        for quantity in ("t1", "t2"):
            medians = [
                f"{np.median(d):.3f} ({d.size})" if d.size else "-"
                for d in tier_differences(runs, q1, q2, quantity).values()
            ]
            lines.append(f"  {quantity.upper()} at {q1:g} W, {q2:g} W: " + " | ".join(medians))
    lines.append("  tiers: " + " | ".join(t.replace("\n", " ") for t in TIERS))

    cells = [c for run in grids_of(runs) for c in cell_partials(run)]
    lines.append(
        "\nLocal partials of heater power by cell, median [25th, 75th percentile] in W/K; dq1_dq2 is the "
        "change in Q1 at fixed T1 per W on stage 2, in W/W"
    )
    for label, chosen in [
        ("one compressor", [c for c in cells if c["run"].compressors == 1]),
        ("two compressors", [c for c in cells if c["run"].compressors == 2]),
    ]:
        for key, _, x_key, _ in [*PARTIALS, ("dq1_dq2", "", "t2", ())]:
            edges = [40, 45, 50, 55, 60, 65] if x_key == "t1" else [2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
            text = []
            for low, high in itertools.pairwise(edges):
                v = np.array([c[key] for c in chosen if low <= c[x_key] < high])
                if v.size >= 3:
                    text.append(
                        f"{low:g}-{high:g} K: {np.median(v):+.3f} [{np.percentile(v, 25):+.3f}, "
                        f"{np.percentile(v, 75):+.3f}] n={v.size}"
                    )
            lines.append(f"  {label}, {key} by {x_key.upper()}: " + "; ".join(text))

    lines.append("\nCollapse tests: rms and largest residual in W (points, runs)")
    for name, group in collapse_groups(runs).items():
        lines.append(f"  {name}")
        for label, stage, run_slope, t2_terms in COLLAPSE_TESTS:
            fit = collapse(group, stage, run_slope, t2_terms)
            lines.append(
                f"    {label}: rms {np.sqrt(np.mean(fit['residual'] ** 2)):.3f}, "
                f"max {np.abs(fit['residual']).max():.3f} ({fit['power'].size}, {len(group)})"
            )
        fit = collapse(group, stage=2)
        lines.append(
            f"    stage-2 shared curve: Q2 = {fit['shared'][0]:.3f} (T2 - 2 K) + {fit['shared'][1]:.3f} "
            f"(T2 - 2 K)^2 - constant; constants {fit['constant'].min():.3f} to {fit['constant'].max():.3f} W"
        )

    lines.append(
        "\nPer-run responses: T1 rise from Q2 = 0 to 0.5 W at Q1 = 0; T1 rise from Q1 = 0 to 20 W at Q2 = 0;"
        " T2 rise from Q2 = 0 to 0.5 W at Q1 = 0 (K)"
    )
    for run in grids_of(runs):
        rises = [
            value(run, 0, 0.5, "t1") - value(run, 0, 0, "t1"),
            value(run, 20, 0, "t1") - value(run, 0, 0, "t1"),
            value(run, 0, 0.5, "t2") - value(run, 0, 0, "t2"),
        ]
        lines.append(
            f"  {run.label:32s} {run.variant:26s} compressors {run.compressors}: "
            + "  ".join(f"{r:6.2f}" for r in rises)
        )
    return "\n".join(lines) + "\n"


def main() -> None:
    """Draws every figure and writes the summary."""
    (HERE / "figures").mkdir(exist_ok=True)
    runs = load_runs()
    plot_grids(runs)
    plot_standard_points(runs)
    plot_variability(runs)
    plot_own_response(runs, stage=1)
    plot_own_response(runs, stage=2)
    plot_cross_response(runs)
    plot_capacity_partials(runs)
    plot_capacity_ratio(runs)
    plot_pulse_tube(runs)
    plot_valve_settings(runs)
    plot_digitizing_check(runs)
    plot_collapse(runs)
    (HERE / "summary.txt").write_text(summarize(runs), encoding="utf-8")


if __name__ == "__main__":
    main()
