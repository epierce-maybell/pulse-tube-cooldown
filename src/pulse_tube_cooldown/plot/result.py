"""Plots of a Result's data."""

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from matplotlib.axes import Axes

__all__ = ["plot_heat_flow", "plot_temperature"]


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


def _mark_phase_switches(data: xr.Dataset, ax: Axes) -> None:
    """Marks each phase transition at its first sample in the new phase."""
    hours = data.time.values / 3600.0
    phases = data.phase.values
    for i in np.flatnonzero(phases[1:] != phases[:-1]):
        ax.axvline(hours[i + 1], color="0.4", linestyle="--", linewidth=1)
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
