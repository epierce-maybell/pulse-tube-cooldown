"""Tests for the canonical head's load-map fit set, `data/test_stand/fit_set.csv`, against the supplied index."""

import csv
from datetime import datetime
from pathlib import Path

import numpy as np
import pytest

from pulse_tube_cooldown.measurements.test_stand import read_load_map, read_load_map_set

_DATA = Path(__file__).parents[2] / "data/test_stand"


@pytest.fixture(scope="module")
def fit_set():
    with open(_DATA / "fit_set.csv", encoding="utf-8") as file:
        return list(csv.DictReader(file))


def test_fit_set_selects_one_map_per_head(fit_set):
    selected = {
        row["serial"]: (row["date"], int(row["compressors"]), float(row["valve_frequency"]), row["mli"])
        for row in fit_set
    }

    assert len(fit_set) == len(selected)
    assert selected == {
        "S0007": ("2026-04-14", 1, 1.0, "full"),
        "S0008": ("2026-05-08", 1, 1.0, "full"),
        "S0009": ("2026-05-14", 1, 1.0, "full"),
        "S0019": ("2026-08-07", 1, 1.1, "unrecorded"),
        "S0023": ("2026-06-19", 1, 1.1, "full"),
        "S0028": ("2026-07-02", 2, 1.1, "full"),
        "S0034": ("2026-08-12", 1, 1.0, "none"),
        "S0035": ("2026-09-21", 1, 1.1, "unrecorded"),
    }


def test_each_map_traces_to_one_denver_pt1_5_index_row(fit_set):
    with open(_DATA / "load_map_index_2026-11-05.csv", encoding="latin-1") as file:
        index = list(csv.DictReader(file))

    for row in fit_set:
        test_date = datetime.strptime(row["date"], "%Y-%m-%d").strftime("%m-%d-%Y")
        matches = [
            entry
            for entry in index
            if entry["Serial Number"] == row["serial"]
            and entry["Test Date (start)"] == test_date
            and entry["Test Stand"].startswith("Denver")
        ]
        assert len(matches) == 1, row["serial"]
        assert matches[0]["PT Type"] == "PT1.5"


def test_each_map_file_is_on_file(fit_set):
    for row in fit_set:
        assert read_load_map(_DATA / row["file"]).sizes["point"] > 0, row["serial"]


def test_renamed_s0023_map_holds_its_archive_points():
    points = read_load_map_set(_DATA / "fit_set.csv")
    s0023 = points.where(points.serial == "S0023", drop=True)

    np.testing.assert_array_equal(s0023.Q1, [20.0, 20.0, 20.0, 0.0, 0.0, 0.0])
    assert float(s0023.T1[0]) == 58.85265388116228
    assert float(s0023.T2[-1]) == 2.5566764625725638
