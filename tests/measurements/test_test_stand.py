"""Tests for the readers of pulse-tube test-stand records, against rows of SN09's files."""

from pathlib import Path

import numpy as np
import pytest

from pulse_tube_cooldown.measurements.test_stand import read_cooldown, read_load_map

_SN09 = Path(__file__).parents[2] / "data" / "sn09" / "pulse_tube"


@pytest.fixture(scope="module")
def cooldown():
    return read_cooldown(_SN09 / "cooldown_2026-05-14.dat")


@pytest.fixture(scope="module")
def load_map():
    return read_load_map(_SN09 / "load_map_2026-05-14.dat")


def test_cooldown_holds_every_row(cooldown):
    # 7349 lines, of which the first four are the calibration curves and the column line.
    assert cooldown.sizes["time"] == 7345


def test_cooldown_first_row(cooldown):
    first = cooldown.isel(time=0)

    assert first.time.values == np.datetime64("2026-05-14T10:16:45")
    assert float(first.T1) == 262.05218395483047
    assert float(first.TPR) == 294.2471368091581
    assert float(first.T2) == 308.6607942011764


def test_cooldown_keeps_dropouts_as_nan(cooldown):
    row = cooldown.sel(time="2026-05-14T10:39:09")

    assert np.isnan(row.T1.item())
    assert float(row.TPR) == 280.8115729003965
    assert float(row.T2) == 194.36424320694445


def test_cooldown_clock_runs_through_appended_run(cooldown):
    # The elapsed column restarts at 14:19:17; the timestamps do not.
    assert np.all(np.diff(cooldown.time.values) > np.timedelta64(0, "s"))
    assert float(cooldown.sel(time="2026-05-14T14:19:17").T1) == 39.837141689921935
    assert cooldown.time.values[-1] == np.datetime64("2026-05-14T14:21:51")


def test_cooldown_temperatures_are_in_kelvin(cooldown):
    for channel in ("T1", "TPR", "T2"):
        assert cooldown[channel].attrs["units"] == "K"


def test_load_map_holds_each_point(load_map):
    np.testing.assert_array_equal(load_map.Q1, [20.0, 20.0, 20.0, 0.0, 0.0, 0.0])
    np.testing.assert_array_equal(load_map.Q2, [1.0, 0.5, 0.0, 1.0, 0.5, 0.0])
    np.testing.assert_array_equal(
        load_map.T1,
        [
            55.0703336630144,
            54.08703004814028,
            53.02827312066305,
            40.46756386544182,
            39.20308102292346,
            38.04663390781175,
        ],
    )
    np.testing.assert_array_equal(
        load_map.TPR,
        [
            35.840206090877665,
            34.86643978603753,
            34.17218776222264,
            33.19724558550911,
            33.16450290354262,
            33.146922380326195,
        ],
    )
    np.testing.assert_array_equal(
        load_map.T2,
        [
            3.983430685634459,
            3.3115453115640037,
            2.4056824622684445,
            3.9470913090330555,
            3.2144582049364914,
            2.312348322774008,
        ],
    )


def test_load_map_units(load_map):
    assert load_map.Q1.attrs["units"] == "W"
    assert load_map.Q2.attrs["units"] == "W"
    for channel in ("T1", "TPR", "T2"):
        assert load_map[channel].attrs["units"] == "K"
