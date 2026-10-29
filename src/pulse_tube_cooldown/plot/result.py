"""Plots of a Result's data."""

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from matplotlib.axes import Axes

__all__ = ["plot_temperature"]


def plot_temperature(data: xr.Dataset, /, *, ax: Axes | None = None, log: bool = True) -> Axes:
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

    ax.set_yscale("log" if log else "linear")
    ax.set_xlabel("Time (h)")
    ax.set_ylabel(f"Temperature ({data.temperature.attrs['units']})")
    ax.legend(title="Node")
    return ax
