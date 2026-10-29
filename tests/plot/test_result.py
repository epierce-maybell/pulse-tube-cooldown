"""Tests for plotting a Result's data."""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest
import xarray as xr

from pulse_tube_cooldown.plot.result import plot_temperature


def _three_phase_data():
    # Switches at 1 h and 2 h, each with a duplicate sample for both sides.
    time = np.array([0.0, 1800.0, 3600.0, 3600.0, 5400.0, 7200.0, 7200.0, 9000.0])
    stage = np.array([300.0, 150.0, 70.0, 70.0, 20.0, 10.0, 10.0, 5.0])
    bath = np.full(time.size, 4.0)
    return xr.Dataset(
        {
            "phase": ("time", ["precool"] * 3 + ["pushing"] * 3 + ["condensing"] * 2),
            "temperature": (("time", "node"), np.column_stack([stage, bath]), {"units": "K"}),
        },
        coords={"time": ("time", time, {"units": "s"}), "node": ["stage", "bath"]},
    )


@pytest.fixture
def ax():
    fig, ax = plt.subplots()
    yield ax
    plt.close(fig)


def _node_lines(ax):
    return {line.get_label(): line for line in ax.get_lines() if line.get_label() in {"stage", "bath"}}


def test_draws_each_node_temperature_against_hours(ax):
    plot_temperature(_three_phase_data(), ax=ax)

    lines = _node_lines(ax)
    assert set(lines) == {"stage", "bath"}
    np.testing.assert_allclose(lines["stage"].get_xdata(), [0.0, 0.5, 1.0, 1.0, 1.5, 2.0, 2.0, 2.5])
    np.testing.assert_allclose(lines["stage"].get_ydata(), [300.0, 150.0, 70.0, 70.0, 20.0, 10.0, 10.0, 5.0])
    np.testing.assert_allclose(lines["bath"].get_ydata(), 4.0)


def test_labels_axes_with_units_and_legend_with_node_names(ax):
    plot_temperature(_three_phase_data(), ax=ax)

    assert ax.get_xlabel() == "Time (h)"
    assert ax.get_ylabel() == "Temperature (K)"
    assert {text.get_text() for text in ax.get_legend().get_texts()} == {"stage", "bath"}


def test_temperature_axis_is_logarithmic_by_default(ax):
    plot_temperature(_three_phase_data(), ax=ax)

    assert ax.get_yscale() == "log"


def test_temperature_axis_is_linear_on_request(ax):
    plot_temperature(_three_phase_data(), ax=ax, log=False)

    assert ax.get_yscale() == "linear"


def test_marks_each_phase_switch_with_both_phase_names(ax):
    plot_temperature(_three_phase_data(), ax=ax)

    nodes = _node_lines(ax).values()
    markers = sorted(line.get_xdata()[0] for line in ax.get_lines() if line not in nodes)
    np.testing.assert_allclose(markers, [1.0, 2.0])
    texts = [text.get_text() for text in ax.texts]
    assert any("precool" in text and "pushing" in text for text in texts)
    assert any("pushing" in text and "condensing" in text for text in texts)


def test_single_phase_has_no_switch_marker(ax):
    data = _three_phase_data().isel(time=slice(0, 3))

    plot_temperature(data, ax=ax)

    assert len(ax.get_lines()) == 2
    assert not ax.texts


def test_marks_phase_switches_in_a_reopened_netcdf_file(ax, tmp_path):
    path = tmp_path / "run.nc"
    _three_phase_data().to_netcdf(path, engine="netcdf4")

    with xr.open_dataset(path) as saved:
        plot_temperature(saved, ax=ax)

    nodes = _node_lines(ax).values()
    markers = sorted(line.get_xdata()[0] for line in ax.get_lines() if line not in nodes)
    np.testing.assert_allclose(markers, [1.0, 2.0])
    assert ax.get_ylabel() == "Temperature (K)"


def test_returns_the_axes_it_drew_on(ax):
    assert plot_temperature(_three_phase_data(), ax=ax) is ax


def test_draws_on_a_new_figure_without_axes():
    ax = plot_temperature(_three_phase_data())
    try:
        assert set(_node_lines(ax)) == {"stage", "bath"}
    finally:
        plt.close(ax.figure)
