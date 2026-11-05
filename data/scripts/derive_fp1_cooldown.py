"""Writes a cooldown record from an InfluxDB annotated-CSV export of FP1's thermometer log.

The export holds the `thermometer` measurement's `temperature` field, one row per device and second, over the
query window. The record is in the form of `data/sn09/cooldown_2026-06-04.dat`: a header line, then tab-separated
rows of `_time` in POSIX seconds on each whole minute of the window, followed by one `<device>_temperature`
column per logged device holding that minute's mean reading in K. A minute without samples for a device is NaN.

Usage: `uv run --with pandas python data/scripts/derive_fp1_cooldown.py <export.csv> [<output.dat>]`.
The output defaults to `data/fp1/cooldown_<date>.dat`, dated by the first sample's UTC day.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

BIN = "60s"


def read_export(path):
    """Returns each device's readings in K indexed by sample time, from an annotated CSV export."""
    frame = pd.read_csv(path, comment="#", usecols=["_time", "_value", "device"])
    frame["_time"] = pd.to_datetime(frame["_time"], utc=True)
    return frame


def bin_means(frame, bin_width=BIN):
    """Returns each device's mean reading over every bin of the export's window, devices as columns."""
    frame = frame.assign(_time=frame["_time"].dt.floor(bin_width))
    means = frame.pivot_table(index="_time", columns="device", values="_value", aggfunc="mean")
    minutes = pd.date_range(means.index.min(), means.index.max(), freq=bin_width)
    return means.reindex(minutes)


def write_record(means, path):
    """Writes the binned means as a tab-separated record with POSIX-second times."""
    seconds = (means.index - pd.Timestamp(0, tz="UTC")) // pd.Timedelta(seconds=1)
    columns = [f"{device}_temperature" for device in means.columns]
    table = np.column_stack([seconds, means.to_numpy()])
    np.savetxt(
        path,
        table,
        fmt=["%d"] + ["%.6g"] * len(columns),
        delimiter="\t",
        header="\t".join(["_time", *columns]),
        comments="",
    )


def main(argv):
    """Derives the record from the export named on the command line."""
    export = Path(argv[1])
    means = bin_means(read_export(export))
    date = means.index[0].strftime("%Y-%m-%d")
    output = (
        Path(argv[2]) if len(argv) > 2 else Path(__file__).resolve().parents[2] / "data/fp1" / f"cooldown_{date}.dat"
    )
    write_record(means, output)
    print(f"{output}: {len(means)} minutes, {means.index[0]} to {means.index[-1]}, columns {list(means.columns)}")


if __name__ == "__main__":
    main(sys.argv)
