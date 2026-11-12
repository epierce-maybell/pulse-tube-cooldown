"""Simulates a transient validation run and compares it with its measured cooldown and legacy baseline.

A run is a directory holding a `config.toml`. Its `[validation]` table names, relative to the repository root, the
`measured` cooldown conventions and an optional legacy `baseline` trace. Its optional `comparison_nodes` maps a
measured or baseline node name to the model node it is compared with. Its optional `[layout]` table maps each node
name to its `[column, row]` in the graph drawing. Run from the repository root:
    uv run python -m simulations.validation.transient.run simulations/validation/transient/fp01/calibration
"""

import argparse
import dataclasses
import math
import subprocess
import sys
import tomllib
from collections.abc import Callable, Mapping
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import tomli_w
import xarray as xr
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from numpy.typing import NDArray

from pulse_tube_cooldown import configure_model, save_output, simulate_cooldown
from pulse_tube_cooldown.core.types import Model
from pulse_tube_cooldown.measurements import log_temperature_residual, read_refrigerator_cooldown
from pulse_tube_cooldown.plot.graph import plot_graph
from pulse_tube_cooldown.plot.result import plot_heat_flow, plot_measured_temperature


def condensing_time(temperature: xr.DataArray, condition: Callable[[float, NDArray[np.float64]], float]) -> float:
    """Returns the first time `condition` reaches zero, interpolated between samples, or NaN if it never does.

    Args:
        temperature: Temperatures with nodes in the order `condition` reads.
        condition: A procedure's stop condition.
    """
    time = temperature.time.values
    values = np.array(
        [condition(t, row) for t, row in zip(time, temperature.transpose("time", "node").values, strict=True)]
    )
    [crossed] = np.nonzero(values <= 0.0)
    if crossed.size == 0:
        return math.nan
    i = crossed[0]
    if i == 0:
        return float(time[0])
    return float(np.interp(0.0, [values[i], values[i - 1]], [time[i], time[i - 1]]))


def rename_nodes[T: (xr.Dataset, xr.DataArray)](data: T, names: Mapping[str, str]) -> T:
    """Returns `data` with each node in `names` renamed to its value."""
    return data.assign_coords(node=[names.get(str(node), str(node)) for node in data.node.values])


def compare(
    model: Model,
    data: xr.Dataset,
    measured: xr.Dataset,
    measured_condensing: float,
    legacy: xr.DataArray | None = None,
) -> dict[str, dict]:
    """Returns the condensing times and scores of the model and, if given, the legacy baseline.

    Each trace is scored up to the earlier of its condensing time and the measured one.

    Args:
        model: The configured model, whose stop condition defines condensing.
        data: The model's Result data.
        measured: The measured cooldown.
        measured_condensing: The measured condensing time.
        legacy: The legacy baseline's temperatures.
    """
    condition = model.procedure.stop_condition
    traces = {"model": data.temperature}
    if legacy is not None:
        traces["legacy"] = legacy.reindex(node=[node.name for node in model.nodes])

    metrics: dict[str, dict] = {"condensing_time": {"measured": measured_condensing}}
    for name, temperature in traces.items():
        condensing = condensing_time(temperature, condition)
        end = float(np.fmin(condensing, measured_condensing))
        rms, rms_by_node = log_temperature_residual(temperature, measured.isel(time=measured.time.values <= end))
        metrics["condensing_time"][name] = condensing
        metrics[name] = {
            "condensing_error": condensing - measured_condensing,
            "scoring_end": end,
            "rms": rms,
            "rms_by_node": rms_by_node,
        }
    return metrics


def provenance(root: Path) -> dict[str, str | bool]:
    """Returns the repository's commit and whether its working tree has uncommitted changes."""

    def git(*args: str) -> str:
        return subprocess.run(["git", *args], cwd=root, capture_output=True, text=True, check=True).stdout.strip()

    return {"commit": git("rev-parse", "HEAD"), "uncommitted_changes": bool(git("status", "--porcelain"))}


def write_metrics(metrics: dict[str, dict], path: Path) -> None:
    """Writes the metrics as TOML."""
    path.write_text("# Times in s since compressor start.\n" + tomli_w.dumps(metrics), encoding="utf-8")


def format_metrics(metrics: dict[str, dict]) -> str:
    """Returns the metrics as a table with times in h."""
    traces = [name for name in metrics["condensing_time"] if name != "measured"]
    nodes = list(metrics["model"]["rms_by_node"])
    columns = ["condensing (h)", "error (h)", "rms", *nodes]
    rows = [["measured", f"{metrics['condensing_time']['measured'] / 3600.0:.2f}", "", "", *[""] * len(nodes)]]
    for name in traces:
        trace = metrics[name]
        rows.append(
            [
                name,
                f"{metrics['condensing_time'][name] / 3600.0:.2f}",
                f"{trace['condensing_error'] / 3600.0:+.2f}",
                f"{trace['rms']:.3f}",
                *(f"{trace['rms_by_node'][node]:.3f}" for node in nodes),
            ]
        )
    table = [["", *columns], *rows]
    widths = [max(len(row[i]) for row in table) for i in range(len(columns) + 1)]
    return "\n".join("  ".join(text.rjust(width) for text, width in zip(row, widths, strict=True)) for row in table)


def plot_temperature(
    data: xr.Dataset, measured: xr.Dataset, legacy: xr.DataArray | None, condensing: dict[str, float]
) -> Figure:
    """Plots the model, measured, and any legacy temperatures of each thermometer, marking each condensing time."""
    line_styles = {"model": "-", "measured": ":", "legacy": "-."}
    if legacy is None:
        del line_styles["legacy"]
    figure = plot_measured_temperature(data, measured)
    axes = figure.axes
    axes[0].get_legend().remove()
    for ax, node in zip(axes, measured.node.values, strict=True):
        if legacy is not None and node in legacy.node.values:
            color = ax.get_lines()[0].get_color()
            hours = legacy.time.values / 3600.0
            ax.plot(hours, legacy.sel(node=node).values, color=color, linestyle=line_styles["legacy"])
        for name, time in condensing.items():
            if not math.isnan(time):
                ax.axvline(time / 3600.0, color="black", linestyle=line_styles[name], linewidth=1.0)

    handles = [Line2D([], [], color="0.3", linestyle=style) for style in line_styles.values()]
    handles += [Line2D([], [], color="black", linewidth=1.0, linestyle=style) for style in line_styles.values()]
    labels = [name.capitalize() for name in line_styles] + [f"{name.capitalize()} condensing" for name in line_styles]
    figure.legend(handles, labels, loc="outside right upper")
    return figure


def plot_heat_flows(data: xr.Dataset, nodes: list[str]) -> Figure:
    """Plots each link's heat flow into each of `nodes`, one panel per node."""
    figure, axes = plt.subplots(
        len(nodes), 1, figsize=(12.0, 3.0 * len(nodes)), sharex=True, squeeze=False, layout="constrained"
    )
    for ax, node in zip(axes[:, 0], nodes, strict=True):
        plot_heat_flow(data.sel(node=[node]), ax=ax)
        ax.legend(loc="center left", bbox_to_anchor=(1.0, 0.5), fontsize="small", title=node)
        ax.set_xlabel("")
    axes[-1, 0].set_xlabel("Time (h)")
    return figure


def simulate(model: Model, path: Path) -> None:
    """Simulates `model` to its time limit, ignoring its stop condition, and saves the Result to `path`."""
    procedure = dataclasses.replace(model.procedure, stop_condition=None)
    result = simulate_cooldown(dataclasses.replace(model, procedure=procedure))
    save_output(result, path)
    if result.status == "solver_failure":
        print(f"The solver failed: {result.data.attrs['message']}", file=sys.stderr)


def run(directory: Path, *, plot_only: bool = False) -> dict[str, dict]:
    """Simulates and compares the run in `directory`.

    Writes `result.nc`, `metrics.toml`, `graph.png`, `temperature.png`, and `heat_flow.png` into `directory`.

    Args:
        directory: The run directory.
        plot_only: Reuse the directory's `result.nc` instead of simulating.

    Returns:
        The metrics.
    """
    root = Path(__file__).parents[3]
    config_path = directory / "config.toml"
    with config_path.open("rb") as file:
        config = tomllib.load(file)
    validation = config["validation"]
    model = configure_model(config_path)

    result_path = directory / "result.nc"
    if not plot_only:
        simulate(model, result_path)
    data = xr.load_dataset(result_path)

    comparison_nodes = validation.get("comparison_nodes", {})
    conventions_path = root / validation["measured"]
    with conventions_path.open("rb") as file:
        conventions = tomllib.load(file)
    measured = rename_nodes(read_refrigerator_cooldown(conventions_path), comparison_nodes)
    measured_condensing = (conventions["condensing_begins"] - conventions["compressor_start"]).total_seconds()
    legacy = None
    if "baseline" in validation:
        legacy = rename_nodes(xr.load_dataset(root / validation["baseline"]).temperature, comparison_nodes)

    metrics = compare(model, data, measured, measured_condensing, legacy)
    write_metrics({"provenance": provenance(root), **metrics}, directory / "metrics.toml")
    print(format_metrics(metrics))

    figures = {
        "graph": plot_graph(model, config.get("layout"), sensors=measured.node.values),
        "temperature": plot_temperature(data, measured, legacy, metrics["condensing_time"]),
        "heat_flow": plot_heat_flows(data, list(measured.node.values)),
    }
    for name, figure in figures.items():
        figure.savefig(directory / f"{name}.png", dpi=150)
        plt.close(figure)
    return metrics


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("directory", type=Path, help="run directory holding config.toml")
    parser.add_argument("--plot-only", action="store_true", help="reuse result.nc instead of simulating")
    args = parser.parse_args()
    plt.switch_backend("Agg")
    run(args.directory, plot_only=args.plot_only)
