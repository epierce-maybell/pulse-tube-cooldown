"""Population cooldowns, load maps, and cell-derived stage responses on one page.

Temperatures are in K, times in s from compressor start, and powers in W.
Q1 and Q2 in input datasets are applied heater powers. Plot derivatives use
net non-heater heat into each stage: dQ/dT = -inverse(dT/dP), where P denotes
heater power. Internal states and unloaded connections equilibrate at each point.
"""

import csv
import itertools
import tomllib
from collections.abc import Mapping, Sequence
from contextlib import AbstractContextManager
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scienceplots  # noqa: F401 — registers the named Matplotlib styles.
import xarray as xr
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator, NullFormatter, ScalarFormatter

from pulse_tube_cooldown.measurements import read_cooldown, read_load_map


@dataclass(frozen=True, slots=True, kw_only=True)
class Run:
    """One measurement run, with its stand, drive settings, and source record."""

    date: str
    stand: str
    compressors: int
    valve_frequency: float
    source: Path
    note: str
    data: xr.Dataset


@dataclass(frozen=True, slots=True, kw_only=True)
class HeadData:
    """A head's selected runs and the measured stage-1 temperature ceiling in K."""

    serial: str
    cooldown: Run
    load_map: Run
    stage_1_max_temperature: float


@dataclass(frozen=True, slots=True, kw_only=True)
class Prediction:
    """Candidate temperatures for a head's two selected run conditions.

    cooldown contains T1(time) and T2(time), in K, with numeric time in s from
    compressor start. load_map contains Q1(point), Q2(point), T1(point), and
    T2(point), with exactly the measured heater pairs in measurement order.
    Each load-map point is an equilibrium including the stand's parasitic loads.
    """

    cooldown: xr.Dataset
    load_map: xr.Dataset


def load_population(selection: Path, *, root: Path) -> tuple[HeadData, ...]:
    """Reads the runs selected in a TOML file, resolving data paths against root."""
    with selection.open("rb") as file:
        config = tomllib.load(file)
    heads = []
    for entry in config["heads"]:
        runs = {}
        for kind in ("cooldown", "load_map"):
            spec = entry[kind]
            source = root / spec["file"]
            if kind == "cooldown":
                data = read_cooldown(source).sel(time=slice(spec["start"], spec["end"]))
                seconds = (data.time.values - np.datetime64(spec["start"])) / np.timedelta64(1, "s")
                data = data[["T1", "T2"]].assign_coords(time=("time", seconds, {"units": "s"}))
            elif spec["format"] == "stand":
                data = read_load_map(source)[["Q1", "Q2", "T1", "T2"]]
            elif spec["format"] == "digitized":
                data = _read_digitized(source, serial=entry["serial"], spec=spec)
            else:
                raise ValueError(f"Unknown load-map format: {spec['format']!r}.")
            runs[kind] = Run(
                date=spec["date"],
                stand=spec["stand"],
                compressors=spec["compressors"],
                valve_frequency=spec["valve_frequency"],
                source=source,
                note=spec["note"],
                data=data,
            )
        heads.append(
            HeadData(
                serial=entry["serial"],
                cooldown=runs["cooldown"],
                load_map=runs["load_map"],
                stage_1_max_temperature=config["stage_1_max_temperature"],
            )
        )
    return tuple(heads)


def _read_digitized(source: Path, *, serial: str, spec: dict) -> xr.Dataset:
    """Reads selected, documented heater levels from a digitized map."""
    with source.open(encoding="utf-8", newline="") as file:
        rows = [
            row
            for row in csv.DictReader(file)
            if row["head"] == serial
            and row["date"] == spec["date"]
            and row["kind"] == "grid"
            and float(row["Q2"]) in spec["stage_2_powers"]
        ]
    return xr.Dataset(
        {
            name: ("point", [float(row[name]) for row in rows], {"units": "W" if name[0] == "Q" else "K"})
            for name in ("Q1", "Q2", "T1", "T2")
        }
    )


def cell_derivatives(points: xr.Dataset) -> xr.Dataset:
    """Estimates the four net heat-flow slopes from complete adjacent heater cells.

    Opposite temperature-response edges are averaged, then the resulting 2 by 2
    dT/dP matrix is inverted and negated. Repeated measurements at a heater pair
    are averaged. Incomplete or non-finite cells have no derivative estimate.
    Returns cell-center T1 and T2 in K, and dQ1dT1, dQ1dT2, dQ2dT1, dQ2dT2
    in W/K, with both heater intervals recorded as coordinates in W.
    """
    powers = [np.unique(points[name].values) for name in ("Q1", "Q2")]
    cells, centers, intervals = [], [], []
    for (a, b), (c, d) in itertools.product(itertools.pairwise(powers[0]), itertools.pairwise(powers[1])):
        corners = {}
        for p1, p2 in itertools.product((a, b), (c, d)):
            selected = points.where((points.Q1 == p1) & (points.Q2 == p2), drop=True)
            corners[p1, p2] = np.array([selected[name].mean().item() for name in ("T1", "T2")])
        if not np.isfinite(list(corners.values())).all():
            continue
        along_1 = (corners[b, c] - corners[a, c] + corners[b, d] - corners[a, d]) / (2 * (b - a))
        along_2 = (corners[a, d] - corners[a, c] + corners[b, d] - corners[b, c]) / (2 * (d - c))
        cells.append(-np.linalg.inv(np.column_stack((along_1, along_2))).ravel())
        centers.append(np.mean(list(corners.values()), axis=0))
        intervals.append((a, b, c, d))
    slopes = np.asarray(cells).reshape(-1, 4)
    temperature = np.asarray(centers).reshape(-1, 2)
    bounds = np.asarray(intervals).reshape(-1, 4)
    return xr.Dataset(
        {
            **{
                name: ("cell", slopes[:, i], {"units": "W/K"})
                for i, name in enumerate(("dQ1dT1", "dQ1dT2", "dQ2dT1", "dQ2dT2"))
            },
            **{name: ("cell", temperature[:, i], {"units": "K"}) for i, name in enumerate(("T1", "T2"))},
        },
        coords={
            name: ("cell", bounds[:, i], {"units": "W"})
            for i, name in enumerate(("P1_low", "P1_high", "P2_low", "P2_high"))
        },
    )


def _check_prediction(head: HeadData, prediction: Prediction) -> None:
    """Checks candidate output at the plotting boundary."""
    for data, dimension in ((prediction.cooldown, "time"), (prediction.load_map, "point")):
        for stage in ("T1", "T2"):
            if data[stage].dims != (dimension,) or not np.all(np.isfinite(data[stage]) & (data[stage] > 0)):
                raise ValueError(f"{head.serial}: candidate {stage} must contain finite positive K along {dimension}.")
    time = prediction.cooldown.time.values
    if not np.issubdtype(time.dtype, np.number) or len(time) < 2 or not np.isfinite(time).all():
        raise ValueError(f"{head.serial}: candidate cooldown time must be numeric seconds from compressor start.")
    if time[0] != 0 or not (np.diff(time) > 0).all():
        raise ValueError(f"{head.serial}: candidate cooldown time must increase strictly from zero.")
    for name in ("Q1", "Q2"):
        if not np.array_equal(prediction.load_map[name].values, head.load_map.data[name].values):
            raise ValueError(f"{head.serial}: candidate map must use the measured heater pairs in measurement order.")


def _draw_map(ax: Axes, points: xr.Dataset, *, color: str, candidate: bool, marker: str) -> None:
    """Draws constant-power curves and filled data or open model points."""
    for fixed, varying in (("Q1", "Q2"), ("Q2", "Q1")):
        for power in np.unique(points[fixed]):
            line = points.where(points[fixed] == power, drop=True).sortby(varying)
            ax.plot(
                line.T1,
                line.T2,
                color=color,
                linestyle="--" if candidate else "-",
                linewidth=0.85 if candidate else 0.5,
                alpha=0.8 if candidate else (0.55 if fixed == "Q1" else 0.25),
                zorder=1,
            )
    ax.scatter(
        points.T1,
        points.T2,
        facecolor="none" if candidate else color,
        edgecolor=color,
        marker=marker,
        s=19 if candidate else 13,
        linewidths=0.75 if candidate else 0.4,
        zorder=3 if candidate else 2,
    )


def comparison_style() -> AbstractContextManager[None]:
    """Applies the Nature preset and the comparison page's stylesheet."""
    return plt.style.context(
        [
            "science",
            "no-latex",
            "nature",
            Path(__file__).with_name("page.mplstyle"),
        ]
    )


def plot_comparison(
    heads: Sequence[HeadData],
    *,
    predictions: Mapping[str, Prediction] | None = None,
    candidate_name: str = "Candidate",
) -> Figure:
    """Compares transient and equilibrium behavior across the selected population.

    The diagnostic question is whether a candidate represents the observed
    cooldowns, fixed points, and local stage responses with shared parameters.
    Predictions use each run's recorded stand and drive conditions.

    records.toml identifies sources, heater levels, run dates, settings, trimming,
    and setup differences. S0035 spans a hose change. S0038 combines a Denver
    cooldown with a Munster map, selecting only explicitly labeled heater levels.
    No sensor uncertainty estimates are available. Missing sensor readings remain
    gaps; measured T1 is shown only below its recorded reliability ceiling.
    Scatter axes follow the data; cross-partial panels include the zero reference.
    Colors use the first six entries of Matplotlib's tab10 cycle. Marker shapes
    match the head identity across panels.
    """
    if predictions is not None:
        if set(predictions) != {head.serial for head in heads}:
            raise ValueError("Candidate predictions must cover exactly the selected heads.")
        for head in heads:
            _check_prediction(head, predictions[head.serial])

    markers = ("o", "s", "^", "D", "v", "P")
    with comparison_style():
        colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
        figure = plt.figure()
        grid = figure.add_gridspec(
            3,
            4,
            left=0.07,
            right=0.98,
            bottom=0.075,
            top=0.885,
            hspace=0.65,
            wspace=0.65,
        )
        cooldown_axes = [figure.add_subplot(grid[0, :2]), figure.add_subplot(grid[0, 2:])]
        map_ax = figure.add_subplot(grid[1:, :2])
        derivative_axes = [figure.add_subplot(grid[row, col]) for row in (1, 2) for col in (2, 3)]
        names = ("dQ1dT1", "dQ1dT2", "dQ2dT1", "dQ2dT2")
        horizontal = ("T1", "T2", "T1", "T2")
        handles = []

        for index, head in enumerate(heads):
            color, marker = colors[index % len(colors)], markers[index % len(markers)]
            label = head.serial
            if head.cooldown.compressors != 1:
                label += f" ({head.cooldown.compressors} compressors)"
            handles.append(Line2D([], [], color=color, marker=marker, markersize=3.5, linewidth=1, label=label))
            for stage, ax in zip(("T1", "T2"), cooldown_axes, strict=True):
                measured = head.cooldown.data[stage]
                if stage == "T1":
                    measured = measured.where(measured < head.stage_1_max_temperature)
                ax.plot(head.cooldown.data.time / 3600, measured, color=color, linewidth=1)
                if predictions is not None:
                    predicted = predictions[head.serial].cooldown
                    ax.plot(predicted.time / 3600, predicted[stage], color=color, linestyle="--", linewidth=1.2)

            _draw_map(map_ax, head.load_map.data, color=color, candidate=False, marker=marker)
            measured_cells = cell_derivatives(head.load_map.data)
            for name, coordinate, ax in zip(names, horizontal, derivative_axes, strict=True):
                ax.scatter(
                    measured_cells[coordinate],
                    measured_cells[name],
                    color=color,
                    s=20,
                    marker=marker,
                    linewidths=0.4,
                )

            if predictions is not None:
                predicted = predictions[head.serial].load_map
                _draw_map(map_ax, predicted, color=color, candidate=True, marker=marker)
                predicted_cells = cell_derivatives(predicted)
                for name, coordinate, ax in zip(names, horizontal, derivative_axes, strict=True):
                    ax.scatter(
                        predicted_cells[coordinate],
                        predicted_cells[name],
                        edgecolor=color,
                        facecolor="none",
                        s=26,
                        marker=marker,
                        linewidths=0.8,
                    )

        duration = max(float(head.cooldown.data.time[-1]) for head in heads) / 3600
        for stage, ax in enumerate(cooldown_axes, start=1):
            letter = "a" if stage == 1 else "b"
            ax.set(
                title=f"{letter}   Stage {stage} cooldown",
                xlabel="Time since start (h)",
                ylabel=f"T{stage} (K)",
            )
            ax.set_yscale("log")
            ax.set_yticks((40, 100, 200) if stage == 1 else (3, 10, 30, 100, 300))
            ax.yaxis.set_major_formatter(ScalarFormatter())
            ax.yaxis.set_minor_formatter(NullFormatter())
            ax.tick_params(which="minor", left=False)
            ax.set_xlim(0, duration)
            ax.xaxis.set_major_locator(MaxNLocator(nbins=4, integer=True))
            ax.grid(axis="y", color="0.9", linewidth=0.4)
            ax.set_axisbelow(True)

        map_ax.set(title="c   Load map", xlabel="T1 (K)", ylabel="T2 (K)")
        tick_steps = [1, 2, 2.5, 5, 10]
        map_ax.xaxis.set_major_locator(MaxNLocator(nbins=5, steps=tick_steps))
        map_ax.yaxis.set_major_locator(MaxNLocator(nbins=5, steps=tick_steps))
        for letter, name, coordinate, ax in zip("defg", names, horizontal, derivative_axes, strict=True):
            ax.set(
                title=f"{letter}   dQ{name[2]}/dT{name[-1]} (W/K)",
                xlabel=f"Cell mean {coordinate} (K)",
            )
            ax.xaxis.set_major_locator(MaxNLocator(nbins=4, steps=tick_steps))
            ax.yaxis.set_major_locator(MaxNLocator(nbins=4, steps=tick_steps))
            if name in ("dQ1dT2", "dQ2dT1"):
                ax.axhline(0, color="0.65", linewidth=0.6, zorder=0)
            ax.margins(x=0.12, y=0.12)

        figure.legend(
            handles=handles,
            loc="upper center",
            bbox_to_anchor=(0.52, 0.995),
            ncols=len(heads),
            frameon=False,
            handlelength=1.5,
            columnspacing=1.6,
        )
        if predictions is not None:
            comparison_handles = [
                Line2D([], [], color="0.25", marker="o", markersize=3.5, linewidth=0.8, label="Data"),
                Line2D(
                    [],
                    [],
                    color="0.25",
                    marker="o",
                    markerfacecolor="white",
                    markersize=3.5,
                    linestyle="--",
                    linewidth=0.8,
                    label=candidate_name,
                ),
            ]
            figure.legend(
                handles=comparison_handles,
                loc="upper right",
                bbox_to_anchor=(0.985, 0.955),
                ncols=2,
                frameon=False,
                handlelength=2,
            )
    return figure
