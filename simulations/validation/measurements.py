"""SN09's and FP1's measured cooldowns and the log-temperature residual of a run against them.

Each system's record and comparison conventions are in `data/<system>/cooldown.toml`.
"""

import tomllib
from pathlib import Path

import numpy as np
import xarray as xr

DATA = Path(__file__).parents[2] / "data"


def read_measurements(system: str) -> xr.Dataset:
    """Reads `temperature(time, node)` in K over the window, time in s since compressor start, and `scored(node)`.

    Readings that are NaN, not positive, or in a run of 10 or more identical consecutive samples are NaN.
    """
    with (DATA / system / "cooldown.toml").open("rb") as file:
        conventions = tomllib.load(file)
    record = DATA / system / conventions["record"]
    with record.open(encoding="utf-8") as file:
        columns = file.readline().split()
    table = np.loadtxt(record, skiprows=1, ndmin=2)

    start = conventions["compressor_start"].timestamp()
    time = table[:, 0] - start
    window = (time >= 0.0) & (time < conventions["condensing_begins"].timestamp() - start)

    channels = conventions["channels"]
    temperature = np.column_stack([drop_invalid(table[:, columns.index(channel)]) for channel in channels])
    nodes = list(channels.values())
    return xr.Dataset(
        {
            "temperature": (("time", "node"), temperature[window], {"units": "K"}),
            "scored": ("node", [node in conventions["scored"] for node in nodes]),
        },
        coords={"time": ("time", time[window], {"units": "s"}), "node": nodes},
    )


def drop_invalid(readings: np.ndarray) -> np.ndarray:
    """Sets NaN, non-positive, and stuck readings to NaN."""
    starts_run = np.r_[True, readings[1:] != readings[:-1]]
    run = np.cumsum(starts_run) - 1
    run_length = np.bincount(run)[run]
    return np.where((readings > 0.0) & (run_length < 10), readings, np.nan)


def residual(simulated: xr.DataArray, measured: xr.Dataset) -> tuple[float, dict[str, float]]:
    """Returns the RMS of `ln(T_simulated / T_measured)` pooled over the scored nodes, and for each scored node.

    `simulated` is `temperature(time, node)` with time in s since compressor start, interpolated at each valid
    measured sample linearly in time in log temperature. A run that ends before the window does is not scored.
    """
    if simulated.time[-1] < measured.time[-1]:
        raise ValueError(
            f"Run ends at {float(simulated.time[-1])} s, before the window's {float(measured.time[-1])} s."
        )

    log_ratios = {}
    for node in measured.node.values[measured.scored.values]:
        temperature = measured.temperature.sel(node=node).values
        valid = ~np.isnan(temperature)
        log_simulated = np.interp(measured.time.values[valid], simulated.time.values, np.log(simulated.sel(node=node)))
        log_ratios[str(node)] = log_simulated - np.log(temperature[valid])

    def rms(values):
        return float(np.sqrt(np.mean(values**2)))

    return rms(np.concatenate(list(log_ratios.values()))), {node: rms(ratio) for node, ratio in log_ratios.items()}
