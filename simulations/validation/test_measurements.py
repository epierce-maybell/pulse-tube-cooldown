"""Sanity checks of the measured cooldowns and the residual."""

import math

import numpy as np
import pytest
import xarray as xr
from measurements import drop_invalid, read_measurements, residual


def test_sn09_window_and_mapping():
    """00:21 UTC on 2026-06-04 to 14:14 UTC on 2026-06-06, the minute before condensing."""
    sn09 = read_measurements("sn09")

    assert float(sn09.time[0]) == 0.0
    assert float(sn09.time[-1]) == 222780.0
    assert float(sn09.temperature.isel(time=0).sel(node="PRP")) == 295.2309482758622607
    assert set(sn09.node.values[sn09.scored.values]) == {"PRP", "RGP", "CFP", "STP", "ICP", "DHX"}
    # RGP's stuck run from 00:01 to 00:29 UTC crosses compressor start.
    assert sn09.temperature.sel(node="RGP", time=slice(0.0, 480.0)).isnull().all()
    assert float(sn09.temperature.sel(node="RGP", time=540.0)) == 294.9520015791164838


def test_fp1_window_and_mapping():
    """16:31 UTC on 2026-10-29 to 10:41 UTC on 2026-10-31."""
    fp1 = read_measurements("fp1")

    assert float(fp1.time[-1]) == 151800.0
    assert set(fp1.node.values) == {"PRP", "RGP", "CFP", "STP", "DHX"}
    assert fp1.scored.all()
    assert float(fp1.temperature.isel(time=-1).sel(node="DHX")) == 14.7043


def test_drop_invalid():
    """A run of 9 identical readings is kept; a run of 10 is not."""
    readings = np.array([290.0, 0.0, -2.0, np.nan, 250.0] + [200.0] * 9 + [150.0] + [100.0] * 10)

    expected = [290.0, np.nan, np.nan, np.nan, 250.0] + [200.0] * 9 + [150.0] + [np.nan] * 10
    np.testing.assert_array_equal(drop_invalid(readings), expected)


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

    rms, rms_by_node = residual(simulated, measured)

    assert rms == pytest.approx(math.sqrt(0.0175))
    assert rms_by_node == pytest.approx({"A": 0.1, "B": 0.2})
    with pytest.raises(ValueError):
        residual(simulated.isel(time=slice(0, 1)), measured)
