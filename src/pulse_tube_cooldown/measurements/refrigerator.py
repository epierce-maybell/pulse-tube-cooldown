"""Refrigerator cooldown records and their configured comparison windows."""

import os
import tomllib
from pathlib import Path

import numpy as np
import xarray as xr
from numpy.typing import NDArray

__all__ = ["read_refrigerator_cooldown"]


def read_refrigerator_cooldown(conventions_path: str | os.PathLike[str], /) -> xr.Dataset:
    """Reads a cooldown using its comparison conventions.

    Args:
        conventions_path: TOML file naming the record relative to itself, channel-to-node mappings, scored nodes,
            compressor start, condensing start, and the number of consecutive equal samples treated as stuck.

    Returns:
        Temperature in K over time in s since compressor start, within the configured window, and a scored-node
        mask. Nonpositive, NaN, and stuck readings are NaN. Stuck runs are detected over the complete source record.
    """
    path = Path(conventions_path)
    with path.open("rb") as file:
        conventions = tomllib.load(file)
    record = path.parent / conventions["record"]
    with record.open(encoding="utf-8") as file:
        columns = file.readline().split()
    table = np.loadtxt(record, skiprows=1, ndmin=2)
    start = conventions["compressor_start"].timestamp()
    time = table[:, 0] - start
    window = (time >= 0.0) & (time < conventions["condensing_begins"].timestamp() - start)

    channels = conventions["channels"]
    temperature = np.column_stack(
        [
            _drop_invalid(table[:, columns.index(channel)], stuck_samples=int(conventions["stuck_samples"]))
            for channel in channels
        ]
    )
    nodes = list(channels.values())
    return xr.Dataset(
        {
            "temperature": (("time", "node"), temperature[window], {"units": "K"}),
            "scored": ("node", [node in conventions["scored"] for node in nodes]),
        },
        coords={"time": ("time", time[window], {"units": "s"}), "node": nodes},
    )


def _drop_invalid(readings: NDArray[np.float64], *, stuck_samples: int) -> NDArray[np.float64]:
    starts_run = np.r_[True, readings[1:] != readings[:-1]]
    run = np.cumsum(starts_run) - 1
    run_length = np.bincount(run)[run]
    return np.where((readings > 0.0) & (run_length < stuck_samples), readings, np.nan)
