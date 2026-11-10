"""Tests for the test-stand cooldown set, `data/test_stand/cooldown_set.csv`, against rows of its archive logs."""

import csv
from pathlib import Path

import numpy as np
import pytest

from pulse_tube_cooldown.measurements.test_stand import read_cooldown

_DATA = Path(__file__).parents[2] / "data/test_stand"


@pytest.fixture(scope="module")
def cooldown_set():
    with open(_DATA / "cooldown_set.csv", encoding="utf-8") as file:
        return list(csv.DictReader(file))


# Each copy's first and last rows in its archive log: time, T1, TPR, T2.
_ENDS = {
    "S0007": (
        ("2026-04-13T09:12:23", 292.3613868610032, 291.5836146927791, 307.8949991490293),
        ("2026-04-13T11:56:33", 40.12280344040372, 33.61229140055599, 2.2587616588069355),
    ),
    "S0009": (
        ("2026-05-14T10:16:45", 262.05218395483047, 294.2471368091581, 308.6607942011764),
        ("2026-05-14T14:21:51", 42.323451782030965, 33.18984619445345, 2.8817696138976325),
    ),
    "S0019": (
        ("2026-08-06T16:15:29", 276.7936670406888, 292.92669955830337, 309.32917943606117),
        ("2026-08-06T20:17:15", 40.708054245064275, 33.91660007950796, 2.352880921835454),
    ),
    "S0028": (
        ("2026-06-30T16:37:46", 275.9377951582203, 296.25578677253327, 309.2356363297636),
        ("2026-06-30T21:12:02", 35.14034634983693, 35.4805722722362, 2.2180318041344536),
    ),
    "S0035": (
        ("2026-09-25T11:10:14", 260.63672275202276, 292.14188347324927, 306.43302449057984),
        ("2026-09-25T14:00:02", 42.60355397556917, 32.956498745291185, 2.4146062443740357),
    ),
    "S0038": (
        ("2026-10-12T14:34:34", 262.4052679988977, 294.75948910707444, 309.18534584974617),
        ("2026-10-12T17:00:02", 45.15543096416311, 33.6587384305478, 2.6931854471000016),
    ),
}


@pytest.mark.parametrize("serial", sorted(_ENDS))
def test_each_copy_runs_between_its_archive_rows(cooldown_set, serial):
    (row,) = [row for row in cooldown_set if row["serial"] == serial]
    cooldown = read_cooldown(_DATA / row["file"])

    for index, (time, t1, tpr, t2) in zip((0, -1), _ENDS[serial], strict=True):
        reading = cooldown.isel(time=index)
        assert reading.time.values == np.datetime64(time)
        np.testing.assert_array_equal([reading.T1, reading.TPR, reading.T2], [t1, tpr, t2])
