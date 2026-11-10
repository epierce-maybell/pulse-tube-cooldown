"""Tests for plotting heat flow from a Result's data."""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest
import xarray as xr

from pulse_tube_cooldown.plot import result as result_plot


@pytest.fixture
def data():
    # The switch is inactive in the middle phase. The manifold touches all three nodes.
    switch = np.array([4.0, 2.0, np.nan, np.nan, 0.0, -1.0])
    manifold = np.array([6.0, 3.0, 3.0, 1.5, 1.5, 0.0])
    flow = np.full((6, 2, 3), np.nan)
    flow[:, 0, 0] = switch
    flow[:, 0, 1] = -switch
    flow[:, 1, 0] = -manifold
    flow[:, 1, 1] = manifold / 3.0
    flow[:, 1, 2] = 2.0 * manifold / 3.0
    return xr.Dataset(
        {
            "phase": ("time", ["precool"] * 2 + ["isolated"] * 2 + ["pushing"] * 2),
            "heat_flow": (("time", "link", "node"), flow, {"units": "W"}),
            "connected": (("link", "node"), [[True, True, False], [True, True, True]]),
        },
        coords={
            "time": ("time", [0.0, 3600.0, 3600.0, 7200.0, 7200.0, 10800.0], {"units": "s"}),
            "link": ["switch", "manifold"],
            "node": ["stage", "shield", "bath"],
        },
    )


@pytest.fixture
def ax():
    fig, ax = plt.subplots()
    yield ax
    plt.close(fig)


def _flow_lines(ax):
    return {line.get_label(): line for line in ax.get_lines() if " → " in line.get_label()}


def _phase_markers(ax):
    return [line.get_xdata()[0] for line in ax.get_lines() if np.ptp(line.get_xdata()) == 0]


def test_draws_signed_flow_for_each_connected_pair_with_inactive_gaps(data, ax):
    result_plot.plot_heat_flow(data, ax=ax)

    expected = {
        "switch → stage": [4.0, 2.0, np.nan, np.nan, 0.0, -1.0],
        "switch → shield": [-4.0, -2.0, np.nan, np.nan, 0.0, 1.0],
        "manifold → stage": [-6.0, -3.0, -3.0, -1.5, -1.5, 0.0],
        "manifold → shield": [2.0, 1.0, 1.0, 0.5, 0.5, 0.0],
        "manifold → bath": [4.0, 2.0, 2.0, 1.0, 1.0, 0.0],
    }
    lines = _flow_lines(ax)
    assert set(lines) == set(expected)
    for label, values in expected.items():
        np.testing.assert_allclose(lines[label].get_xdata(), [0.0, 1.0, 1.0, 2.0, 2.0, 3.0])
        np.testing.assert_allclose(lines[label].get_ydata(), values)


def test_labels_units_sign_convention_and_link_node_pairs(data, ax):
    result_plot.plot_heat_flow(data, ax=ax)

    assert ax.get_xlabel() == "Time (h)"
    assert ax.get_ylabel() == "Heat flow (W; positive into node)"
    assert {text.get_text() for text in ax.get_legend().get_texts()} == set(_flow_lines(ax))


def test_uses_linear_axis_to_show_positive_negative_and_zero_flows(data, ax):
    ax.set_yscale("log")

    result_plot.plot_heat_flow(data, ax=ax)

    assert ax.get_yscale() == "linear"


def test_marks_phase_switches_with_phase_names(data, ax):
    result_plot.plot_heat_flow(data, ax=ax)

    np.testing.assert_allclose(_phase_markers(ax), [1.0, 2.0])
    assert [text.get_text() for text in ax.texts] == ["precool → isolated", "isolated → pushing"]


def test_single_phase_has_no_switch_marker(data, ax):
    result_plot.plot_heat_flow(data.isel(time=slice(0, 2)), ax=ax)

    assert len(_flow_lines(ax)) == 5
    assert not _phase_markers(ax)
    assert not ax.texts


def test_plots_reopened_netcdf_data(data, ax, tmp_path):
    path = tmp_path / "run.nc"
    data.to_netcdf(path, engine="netcdf4")

    with xr.open_dataset(path) as saved:
        result_plot.plot_heat_flow(saved, ax=ax)

    np.testing.assert_allclose(_flow_lines(ax)["switch → stage"].get_ydata(), [4.0, 2.0, np.nan, np.nan, 0.0, -1.0])
    assert ax.get_ylabel() == "Heat flow (W; positive into node)"
    assert len(ax.texts) == 2


def test_creates_axes_when_omitted(data):
    ax = result_plot.plot_heat_flow(data)
    try:
        assert len(_flow_lines(ax)) == 5
    finally:
        plt.close(ax.figure)


def test_plots_selected_links_and_nodes(data, ax):
    result_plot.plot_heat_flow(data.sel(link=["manifold"], node=["shield"]), ax=ax)

    lines = _flow_lines(ax)
    assert set(lines) == {"manifold → shield"}
    np.testing.assert_allclose(lines["manifold → shield"].get_ydata(), [2.0, 1.0, 1.0, 0.5, 0.5, 0.0])
    np.testing.assert_allclose(_phase_markers(ax), [1.0, 2.0])
