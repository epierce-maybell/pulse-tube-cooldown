"""Tests for the readers of pulse-tube test-stand records, against rows of S0009's files."""

from pathlib import Path

import numpy as np
import pytest

from pulse_tube_cooldown.measurements.test_stand import read_cooldown, read_load_map, read_load_map_set

_S0009 = Path(__file__).parents[2] / "data/test_stand/s0009"


@pytest.fixture(scope="module")
def cooldown():
    return read_cooldown(_S0009 / "cooldown_2026-05-14.dat")


@pytest.fixture(scope="module")
def load_map():
    return read_load_map(_S0009 / "load_map_2026-05-14.dat")


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


def test_cooldown_reads_iso_timestamps(tmp_path):
    # The header and first two rows of S0007's 2026-04-13 log.
    path = tmp_path / "cooldown.dat"
    path.write_text(
        "T2: C:\\Experiments\\S0007\\Z1512_CernoxFormat.cof\n"
        "T1/Taux: C:\\Experiments\\S0007\\Z1518_CernoxFormat.cof\n"
        "T3/Taux: C:\\Experiments\\S0007\\Z1518_CernoxFormat.cof\n"
        "Date TimeElapsed(m) T1 TPR T2 Taux Taux2 Preturn ORV1 DIV1A DIV1B ORV2 DIV2A DIV2B Frequency/Hz\n"
        "2026-04-13T09:12:23   0.04   292.3613868610032   291.5836146927791   307.8949991490293   "
        "292.3613868610032    292.3613868610032    -1 - - - - - - 1.1 0 0\n"
        "2026-04-13T09:12:25   0.07   292.3613868610032   291.5836146927791   307.88376721037844   "
        "292.3613868610032    292.3613868610032    -1 - - - - - - 1.1 0 0\n"
    )

    cooldown = read_cooldown(path)

    np.testing.assert_array_equal(
        cooldown.time, np.array(["2026-04-13T09:12:23", "2026-04-13T09:12:25"], dtype="datetime64[s]")
    )
    np.testing.assert_array_equal(cooldown.T1, [292.3613868610032, 292.3613868610032])
    np.testing.assert_array_equal(cooldown.TPR, [291.5836146927791, 291.5836146927791])
    np.testing.assert_array_equal(cooldown.T2, [307.8949991490293, 307.88376721037844])


def test_cooldown_without_a_header_keeps_its_first_row(tmp_path):
    # The archive's day-split logs, such as S0019's 2026-08-06 log, start with data.
    path = tmp_path / "cooldown.dat"
    path.write_text(
        "Thu Aug  6 16:15:29 2026   0.03   276.7936670406888   292.92669955830337   309.32917943606117   "
        "276.7936670406888   276.7936670406888\n"
        "Thu Aug  6 16:15:31 2026   0.07   276.7936670406888   292.92669955830337   309.32917943606117   "
        "276.7936670406888   276.7936670406888\n"
    )

    cooldown = read_cooldown(path)

    assert cooldown.sizes["time"] == 2
    assert cooldown.time.values[0] == np.datetime64("2026-08-06T16:15:29")
    assert float(cooldown.T2[0]) == 309.32917943606117


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


_HEADER = "ORV1:- DIV1A:- DIV1B:- ORV2:-DIV2A:- DIV2B:- f:1.1\nPower Stage 1: [20.0, 0.0]\nPower Stage 2: [0.5, 0.0]\n"


@pytest.fixture
def load_map_set(tmp_path):
    (tmp_path / "s0001").mkdir()
    (tmp_path / "s0001/load_map_2026-01-02.dat").write_text(
        _HEADER
        + "20.0  0.5   120.0   50.0   34.0   3.5   50.0   50.0\n"
        + "0.0  0.0   240.0   38.0   33.0   2.2   38.0   38.0\n"
    )
    (tmp_path / "s0002").mkdir()
    (tmp_path / "s0002/load_map_2026-03-04.dat").write_text(
        _HEADER + "0.0  0.5   120.0   41.0   33.5   3.3   41.0   41.0\n"
    )
    table = tmp_path / "set.csv"
    table.write_text(
        "serial,date,file,source,compressors,valve_frequency,mli\n"
        "S0001,2026-01-02,s0001/load_map_2026-01-02.dat,a.dat,1,1.0,full\n"
        "S0002,2026-03-04,s0002/load_map_2026-03-04.dat,b.dat,2,1.1,unrecorded\n"
    )
    return read_load_map_set(table)


def test_load_map_set_joins_every_point_in_table_order(load_map_set):
    np.testing.assert_array_equal(load_map_set.Q1, [20.0, 0.0, 0.0])
    np.testing.assert_array_equal(load_map_set.Q2, [0.5, 0.0, 0.5])
    np.testing.assert_array_equal(load_map_set.T1, [50.0, 38.0, 41.0])
    np.testing.assert_array_equal(load_map_set.TPR, [34.0, 33.0, 33.5])
    np.testing.assert_array_equal(load_map_set.T2, [3.5, 2.2, 3.3])


def test_load_map_set_points_carry_their_run(load_map_set):
    np.testing.assert_array_equal(load_map_set.serial, ["S0001", "S0001", "S0002"])
    np.testing.assert_array_equal(load_map_set.date, ["2026-01-02", "2026-01-02", "2026-03-04"])
    np.testing.assert_array_equal(load_map_set.compressors, [1, 1, 2])
    np.testing.assert_array_equal(load_map_set.valve_frequency, [1.0, 1.0, 1.1])
    np.testing.assert_array_equal(load_map_set.mli, ["full", "full", "unrecorded"])


def test_load_map_set_units(load_map_set):
    assert load_map_set.valve_frequency.attrs["units"] == "Hz"
    assert load_map_set.Q1.attrs["units"] == "W"
    assert load_map_set.T2.attrs["units"] == "K"
