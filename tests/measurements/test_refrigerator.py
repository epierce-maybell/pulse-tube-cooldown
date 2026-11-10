"""Sanity checks of the measured cooldowns and the residual."""

import math
from pathlib import Path

import numpy as np
import pytest
import xarray as xr

from pulse_tube_cooldown.measurements import log_temperature_residual, read_refrigerator_cooldown

DATA = Path(__file__).parents[2] / "data"


def test_sn09_window_and_mapping():
    """00:21 UTC on 2026-06-04 to 14:14 UTC on 2026-06-06, the minute before condensing."""
    sn09 = read_refrigerator_cooldown(DATA / "sn09" / "cooldown.toml")

    assert float(sn09.time[0]) == 0.0
    assert float(sn09.time[-1]) == 222780.0
    assert float(sn09.temperature.isel(time=0).sel(node="PRP")) == 295.2309482758622607
    assert set(sn09.node.values[sn09.scored.values]) == {"PRP", "RGP", "CFP", "STP", "ICP", "DHX"}
    # RGP's stuck run from 00:01 to 00:29 UTC crosses compressor start.
    assert sn09.temperature.sel(node="RGP", time=slice(0.0, 480.0)).isnull().all()
    assert float(sn09.temperature.sel(node="RGP", time=540.0)) == 294.9520015791164838


def test_fp1_window_and_mapping():
    """16:31 UTC on 2026-10-29 to 10:41 UTC on 2026-10-31."""
    fp1 = read_refrigerator_cooldown(DATA / "fp1" / "cooldown.toml")

    assert float(fp1.time[-1]) == 151800.0
    assert set(fp1.node.values) == {"PRP", "RGP", "CFP", "STP", "DHX"}
    assert fp1.scored.all()
    assert float(fp1.temperature.isel(time=-1).sel(node="DHX")) == 14.7043


def test_reader_drops_nonpositive_nan_and_stuck_readings(tmp_path):
    """The configured threshold drops repeated samples alongside NaN and nonpositive readings."""
    conventions = tmp_path / "cooldown.toml"
    conventions.write_text(
        'record = "record.dat"\n'
        "compressor_start = 2026-01-01T00:00:00Z\n"
        "condensing_begins = 2026-01-01T00:00:20Z\n"
        "stuck_samples = 3\n"
        'scored = ["A"]\n'
        "[channels]\n"
        'A = "A"\n'
        'B = "B"\n',
        encoding="utf-8",
    )
    rows = [
        "seconds A B",
        "1767225600 290 200",
        "1767225601 0 200",
        "1767225602 -2 200",
        "1767225603 nan 200",
        "1767225604 250 200",
        "1767225605 200 300",
        "1767225606 200 300",
        "1767225607 150 300",
        "1767225608 100 400",
        "1767225609 100 400",
        "1767225610 100 400",
    ]
    (tmp_path / "record.dat").write_text("\n".join(rows) + "\n", encoding="utf-8")

    measured = read_refrigerator_cooldown(conventions)

    np.testing.assert_array_equal(
        measured.temperature.sel(node="A"),
        [290.0, np.nan, np.nan, np.nan, 250.0, 200.0, 200.0, 150.0, np.nan, np.nan, np.nan],
    )


def _temperature(time, columns):
    return xr.DataArray(
        np.column_stack(list(columns.values())),
        dims=("time", "node"),
        coords={"time": np.asarray(time, dtype=np.float64), "node": list(columns)},
    )


def test_residual():
    """Pooled: sqrt((3 * 0.1^2 + 0.2^2) / 4).

    A is high by e^0.1 at three samples, the middle one reached by log interpolation; B by e^0.2 at its one
    valid sample; C is unscored.
    """
    measured = xr.Dataset(
        {
            "temperature": _temperature(
                [0.0, 60.0, 120.0], {"A": [400.0, 200.0, 100.0], "B": [300.0, np.nan, np.nan], "C": [1.0] * 3}
            ),
            "scored": ("node", [True, True, False]),
        }
    )
    # The repeated time is a phase switch.
    simulated = _temperature(
        [0.0, 120.0, 120.0, 180.0],
        {
            "A": np.array([400.0, 100.0, 100.0, 100.0]) * math.exp(0.1),
            "B": [300.0 * math.exp(0.2), 1.0, 1.0, 1.0],
            "C": [300.0] * 4,
        },
    )

    rms, rms_by_node = log_temperature_residual(simulated, measured)

    assert rms == pytest.approx(math.sqrt(0.0175))
    assert rms_by_node == pytest.approx({"A": 0.1, "B": 0.2})
    with pytest.raises(ValueError):
        log_temperature_residual(simulated.isel(time=slice(0, 1)), measured)
