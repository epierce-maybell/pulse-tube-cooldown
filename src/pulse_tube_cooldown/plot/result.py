"""Plots of a Result's data."""

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.lines import Line2D

__all__ = ["plot_heat_flow", "plot_measured_temperature", "plot_temperature"]


def plot_temperature(data: xr.Dataset, ax: Axes | None = None, log: bool = True) -> Axes:
    """Plots each node's temperature over time, marking phase switches.

    Args:
        data: A Result's data, either `result.data` or a Dataset reopened from `save_output`.
        ax: Axes to draw on. If omitted, draws on a new figure.
        log: Whether the temperature axis is logarithmic.

    Returns:
        The Axes drawn on, with time in hours and one labelled line per node.
    """
    if ax is None:
        _, ax = plt.subplots()
    hours = data.time.values / 3600.0
    for node in data.node.values:
        ax.plot(hours, data.temperature.sel(node=node).values, label=node)

    _mark_phase_switches(data, ax)
    ax.set_yscale("log" if log else "linear")
    ax.set_xlabel("Time (h)")
    ax.set_ylabel(f"Temperature ({data.temperature.attrs['units']})")
    ax.legend(title="Node")
    return ax


def plot_measured_temperature(data: xr.Dataset, measured: xr.Dataset, *, log: bool = True) -> Figure:
    """Plots each measured thermometer over its simulated node, one panel per thermometer.

    Args:
        data: A Result's data, either `result.data` or a Dataset reopened from `save_output`.
        measured: `temperature(time, node)` on the same time origin as `data`, with each node named in `data`.
        log: Whether the temperature axes are logarithmic.

    Returns:
        The Figure, with panels in `measured`'s order sharing time in hours.
    """
    nodes = measured.node.values
    figure, axes = plt.subplots(
        len(nodes), 1, figsize=(10.0, 2.0 * len(nodes)), sharex=True, squeeze=False, layout="constrained"
    )
    simulated_nodes = list(data.node.values)
    for ax, node in zip(axes[:, 0], nodes, strict=True):
        color = f"C{simulated_nodes.index(node)}"
        ax.plot(data.time.values / 3600.0, data.temperature.sel(node=node).values, color=color)
        ax.plot(measured.time.values / 3600.0, measured.temperature.sel(node=node).values, color=color, linestyle=":")
        _mark_phase_switches(data, ax, label=ax is axes[0, 0])
        ax.set_yscale("log" if log else "linear")
        ax.set_ylabel(f"{node} ({data.temperature.attrs['units']})")

    styles = [Line2D([], [], color="0.3"), Line2D([], [], color="0.3", linestyle=":")]
    axes[0, 0].legend(styles, ["Simulated", "Measured"], loc="upper right")
    axes[-1, 0].set_xlabel("Time (h)")
    return figure


def plot_heat_flow(data: xr.Dataset, ax: Axes | None = None) -> Axes:
    """Plots signed heat flow for each connected link and node, marking phase switches.

    Positive flow enters the labelled node; negative flow leaves it. Inactive links
    appear as gaps in their traces.

    Args:
        data: A Result's data, either `result.data` or a Dataset reopened from `save_output`.
            Select links or nodes with `data.sel(link=[...], node=[...])`.
        ax: Axes to draw on. If omitted, draws on a new figure.

    Returns:
        The Axes drawn on, with time in hours and one labelled line per connected pair.
    """
    if ax is None:
        _, ax = plt.subplots()
    ax.axhline(0.0, color="0.5", linestyle="--", linewidth=1, zorder=0)
    hours = data.time.values / 3600.0
    for link in data.link.values:
        for node in data.node.values:
            if data.connected.sel(link=link, node=node).item():
                ax.plot(hours, data.heat_flow.sel(link=link, node=node).values, label=f"{link} → {node}")

    _mark_phase_switches(data, ax)
    ax.set_yscale("linear")
    ax.set_xlabel("Time (h)")
    ax.set_ylabel(f"Heat flow ({data.heat_flow.attrs['units']}; positive into node)")
    ax.legend(title="Link → node")
    return ax


def _mark_phase_switches(data: xr.Dataset, ax: Axes, label: bool = True) -> None:
    """Marks each phase transition at its first sample in the new phase, optionally naming both phases."""
    hours = data.time.values / 3600.0
    phases = data.phase.values
    for i in np.flatnonzero(phases[1:] != phases[:-1]):
        ax.axvline(hours[i + 1], color="0.4", linestyle="--", linewidth=1)
        if not label:
            continue
        ax.annotate(
            f"{phases[i]} → {phases[i + 1]}",
            (hours[i + 1], 1.0),
            xycoords=("data", "axes fraction"),
            xytext=(3, -3),
            textcoords="offset points",
            va="top",
            fontsize="small",
            color="0.3",
        )
