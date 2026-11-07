"""Readers for the pulse-tube test stand's cooldown traces and load maps, in the format of `data/test_stand.md`.

T1 is stage 1 and T2 stage 2. TPR is the raw column name for T3, the thermometer mounted on the pulse tube itself.
Q1 and Q2 are the stage-1 and stage-2 heater powers. Readings the stand logs as `nan` stay NaN.
"""

import csv
import os
from datetime import datetime
from pathlib import Path

import numpy as np
import xarray as xr

__all__ = ["read_cooldown", "read_load_map", "read_load_map_set"]


def read_cooldown(path: str | os.PathLike[str]) -> xr.Dataset:
    """Reads a cooldown trace.

    Args:
        path: A `cooldown_<date>.dat` file. Its header, where present, ends at the column line, which starts with
            `T1` or `Date`. Each row starts with its timestamp, as `Thu May 14 10:16:45 2026` or as ISO 8601.

    Returns:
        T1, TPR, and T2 in K over the `time` of each row, from its timestamp.
    """
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    start = next((number + 1 for number, line in enumerate(lines) if line.startswith(("T1 ", "Date "))), 0)
    time, readings = zip(*(_split_timestamp(line.split()) for line in lines[start:]), strict=True)
    temperature = np.array([reading[1:4] for reading in readings], dtype=np.float64)
    return xr.Dataset(
        {
            name: ("time", column, {"units": "K"})
            for name, column in zip(("T1", "TPR", "T2"), temperature.T, strict=True)
        },
        coords={"time": np.array(time, dtype="datetime64[s]")},
    )


def _split_timestamp(fields: list[str]) -> tuple[datetime, list[str]]:
    """Splits a cooldown row into its timestamp and the fields after it: elapsed minutes, T1, TPR, T2, and so on."""
    if fields[0][0].isdigit():
        return datetime.fromisoformat(fields[0]), fields[1:]
    return datetime.strptime(" ".join(fields[:5]), "%a %b %d %H:%M:%S %Y"), fields[5:]


def read_load_map(path: str | os.PathLike[str]) -> xr.Dataset:
    """Reads a load map.

    Args:
        path: A `load_map_<date>.dat` file.

    Returns:
        Q1 and Q2 in W, and T1, TPR, and T2 in K, at each steady-state `point`.
    """
    rows = np.loadtxt(path, skiprows=3, usecols=(0, 1, 3, 4, 5), ndmin=2)
    units = {"Q1": "W", "Q2": "W", "T1": "K", "TPR": "K", "T2": "K"}
    return xr.Dataset(
        {name: ("point", column, {"units": units[name]}) for name, column in zip(units, rows.T, strict=True)}
    )


def read_load_map_set(path: str | os.PathLike[str]) -> xr.Dataset:
    """Reads the load maps a set table lists, as one Dataset of points.

    Args:
        path: A CSV with one row per map and columns serial, date, file, source, compressors, valve_frequency, and mli.
            `file` is the map's path relative to the CSV's directory, and `source` its name in the stand's archive.

    Returns:
        Q1, Q2, T1, TPR, and T2 along `point`, in table and measurement order. Each point carries its map's serial,
        date, compressor count, valve frequency in Hz, and MLI state ("full", "reduced", "none", or "unrecorded")
        as coordinates.
    """
    table = Path(path)
    with open(table, encoding="utf-8", newline="") as file:
        rows = list(csv.DictReader(file))

    maps = []
    for row in rows:
        points = read_load_map(table.parent / row["file"])
        count = points.sizes["point"]
        maps.append(
            points.assign_coords(
                serial=("point", [row["serial"]] * count),
                date=("point", [row["date"]] * count),
                compressors=("point", [int(row["compressors"])] * count),
                valve_frequency=("point", [float(row["valve_frequency"])] * count, {"units": "Hz"}),
                mli=("point", [row["mli"]] * count),
            )
        )
    return xr.concat(maps, dim="point")
