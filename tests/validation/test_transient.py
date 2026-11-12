"""The transient validation harness: condensing times, scoring windows, and the run's outputs."""

import math
import shutil
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pytest
import xarray as xr

import simulations.validation.transient.run as transient
from pulse_tube_cooldown import configure_model
from pulse_tube_cooldown.__main__ import main as ptc_main
from pulse_tube_cooldown.library.conditions import condensing_begins

FIXTURES = Path(__file__).parent / "fixtures"
TWO_NODE = FIXTURES / "two_node.toml"
CONDUCTION = FIXTURES / "conduction_sn09.toml"
LN2 = math.log(2.0)
TIME = [0.0, 10.0, 20.0, 30.0, 40.0]


def _trace(time, columns):
    """A temperature(time, node) DataArray from a mapping of node name to its temperatures."""
    values = np.array(list(columns.values()), dtype=float).T
    return xr.DataArray(values, coords={"time": time, "node": list(columns)}, dims=("time", "node"))


@pytest.fixture(scope="module")
def two_node_model():
    return configure_model(TWO_NODE)


@pytest.mark.parametrize(
    ("time", "b", "want"),
    [
        ([0.0, 100.0, 200.0], [40.0, 16.0, 4.0], 150.0),  # B - 10 goes 30, 6, -6: zero midway through 100..200 s.
        ([0.0, 70.0, 140.0, 210.0], [40.0, 5.0, 40.0, 5.0], 60.0),  # Later re-crossings do not move the first.
        ([0.0, 100.0, 200.0], [40.0, 30.0, 10.0], 200.0),  # Landing exactly on zero at a sample.
        ([5.0, 10.0], [8.0, 5.0], 5.0),  # Already past the event at the first sample.
        ([0.0, 100.0, 200.0], [40.0, 30.0, 20.0], math.nan),  # Never reaches the threshold.
        ([0.0, 100.0, 200.0], [math.nan] * 3, math.nan),  # A trace without the node it reads.
    ],
)
def test_condensing_time(two_node_model, time, b, want):
    """The first zero crossing is interpolated linearly; no crossing gives NaN.

    Catches interpolating from the wrong bracket, taking a later crossing, or treating NaN as a crossing.
    """
    condition = condensing_begins(two_node_model.graphs[0], "B", 10.0)
    temperature = _trace(time, {"A": [295.0] * len(time), "B": b})

    assert transient.condensing_time(temperature, condition) == pytest.approx(want, nan_ok=True)


@pytest.fixture(scope="module")
def measured():
    """Measured A at 100 K and B at 40, 30, 20, 10, 10 K at 0, 10, 20, 30, 40 s, both scored."""
    temperature = _trace(TIME, {"A": [100.0] * 5, "B": [40.0, 30.0, 20.0, 10.0, 10.0]})
    return xr.Dataset({"temperature": temperature, "scored": ("node", [True, True])})


@pytest.fixture(scope="module")
def model_data():
    """Model: A is 200 K to 20 s then 400 K; B equals the measured B, reaching 10 K at 30 s."""
    return xr.Dataset(
        {"temperature": _trace(TIME, {"A": [200.0, 200.0, 200.0, 400.0, 400.0], "B": [40.0, 30.0, 20.0, 10.0, 10.0]})}
    )


@pytest.fixture(scope="module")
def legacy():
    """Legacy with an extra node Z listed first: A is 100 K to 20 s then 300 K; B reaches 10 K at 20 s."""
    return _trace(
        TIME,
        {
            "Z": [1.0] * 5,
            "B": [40.0, 30.0, 10.0, 10.0, 10.0],
            "A": [100.0, 100.0, 100.0, 300.0, 300.0],
        },
    )


def test_compare_window_ends_at_measured_condensing(two_node_model, measured, model_data, legacy):
    """A model that condenses after the measurement is scored only up to the measured condensing time.

    Catches scoring samples after the window end (A's ratio is ln 2 to 20 s, ln 4 at 30 s) and a sign flip in the
    condensing error.
    """
    metrics = transient.compare(two_node_model, model_data, measured, 20.0, legacy)

    model = metrics["model"]
    assert metrics["condensing_time"]["model"] == pytest.approx(30.0)
    assert model["condensing_error"] == pytest.approx(10.0)
    assert model["scoring_end"] == pytest.approx(20.0)
    assert model["rms_by_node"] == pytest.approx({"A": LN2, "B": 0.0})
    assert model["rms"] == pytest.approx(LN2 / math.sqrt(2.0))


def test_compare_window_ends_at_model_condensing(two_node_model, measured, model_data, legacy):
    """A model that condenses before the measurement is scored only up to its own condensing time.

    Catches using the measured time alone as the window end.
    """
    metrics = transient.compare(two_node_model, model_data, measured, 35.0, legacy)

    model = metrics["model"]
    assert model["condensing_error"] == pytest.approx(-5.0)
    assert model["scoring_end"] == pytest.approx(30.0)
    assert model["rms_by_node"] == pytest.approx({"A": LN2 * math.sqrt(7.0 / 4.0), "B": 0.0})


def test_compare_legacy_is_reindexed_to_model_nodes(two_node_model, measured, model_data, legacy):
    """Legacy nodes the model lacks are dropped and the rest follow the model's node order.

    Catches reading the condition at the legacy trace's own node order (Z, B, A) and scoring legacy to its own
    window rather than the earlier of its condensing time (20 s) and the measured one (25 s).
    """
    metrics = transient.compare(two_node_model, model_data, measured, 25.0, legacy)

    result = metrics["legacy"]
    assert metrics["condensing_time"]["legacy"] == pytest.approx(20.0)
    assert result["condensing_error"] == pytest.approx(-5.0)
    assert result["scoring_end"] == pytest.approx(20.0)
    assert result["rms_by_node"] == pytest.approx({"A": 0.0, "B": LN2 / math.sqrt(3.0)})
    assert result["rms"] == pytest.approx(LN2 / math.sqrt(6.0))


def test_compare_legacy_missing_a_model_node(two_node_model, measured, model_data, legacy):
    """A legacy trace without the condition's node has NaN condensing and is scored to the measured time.

    Catches an exception or a zero condensing time when legacy lacks a model node.
    """
    metrics = transient.compare(two_node_model, model_data, measured, 25.0, legacy.sel(node=["A", "Z"]))

    assert math.isnan(metrics["condensing_time"]["legacy"])
    assert math.isnan(metrics["legacy"]["condensing_error"])
    assert metrics["legacy"]["scoring_end"] == pytest.approx(25.0)
    assert metrics["legacy"]["rms_by_node"]["A"] == pytest.approx(0.0)


@pytest.fixture(scope="module")
def harness_output(tmp_path_factory):
    """`run` on a directory holding the conduction-only SN09 configuration and the result `run-ptc` made from it."""
    plt.switch_backend("Agg")
    directory = tmp_path_factory.mktemp("transient")
    shutil.copy(CONDUCTION, directory / "config.toml")
    code = ptc_main(
        [
            "--config",
            str(directory / "config.toml"),
            "--output",
            str(directory / "result.nc"),
            "--ignore-stop-condition",
        ]
    )
    assert code == 0
    return directory, transient.run("sn09", directory)


def test_run_writes_outputs(harness_output):
    """`run` writes its metrics and three figures into the directory.

    Catches a missing output file.
    """
    directory, _ = harness_output

    for name in ["metrics.toml", "graph.png", "temperature.png", "heat_flow.png"]:
        assert (directory / name).stat().st_size > 0


def test_run_scores_against_sn09_record(harness_output):
    """The measured condensing time is 2 d 13 h 54 min (222840 s) after compressor start, and the scored nodes are
    the record's.

    Catches reading the wrong record times or scoring the wrong thermometers.
    """
    _, metrics = harness_output

    assert metrics["condensing_time"]["measured"] == 222840.0
    assert 0.0 < metrics["condensing_time"]["model"] < 240000.0
    assert set(metrics["model"]["rms_by_node"]) == {"PRP", "RGP", "CFP", "STP", "ICP", "DHX"}
    assert set(metrics["legacy"]["rms_by_node"]) == set(metrics["model"]["rms_by_node"])


def test_run_depends_only_on_the_directory(harness_output, tmp_path):
    """A copy of the directory's config and result gives the same metrics.

    Catches `run` reading state from anywhere but `config.toml` and `result.nc`.
    """
    directory, _ = harness_output
    again = tmp_path / "again"
    again.mkdir()
    for name in ["config.toml", "result.nc"]:
        shutil.copy(directory / name, again / name)

    transient.run("sn09", again)

    assert (again / "metrics.toml").read_text() == (directory / "metrics.toml").read_text()
