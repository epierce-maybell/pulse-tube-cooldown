"""Tests for plotting measured thermometer traces over their simulated nodes."""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest
import xarray as xr
from matplotlib.colors import to_hex

from pulse_tube_cooldown.plot.result import plot_measured_temperature, plot_temperature


def _simulated():
    # A phase switch at 1 h, with a duplicate sample for both sides.
    time = np.array([0.0, 1800.0, 3600.0, 3600.0, 5400.0, 7200.0])
    return xr.Dataset(
        {
            "phase": ("time", ["precool"] * 3 + ["pushing"] * 3),
            "temperature": (
                ("time", "node"),
                np.column_stack(
                    [
                        [300.0, 150.0, 70.0, 70.0, 20.0, 10.0],
                        [300.0, 200.0, 120.0, 120.0, 60.0, 40.0],
                        [300.0, 280.0, 250.0, 250.0, 200.0, 150.0],
                    ]
                ),
                {"units": "K"},
            ),
        },
        coords={"time": ("time", time, {"units": "s"}), "node": ["stage", "shield", "bath"]},
    )


def _measured():
    # Listed in a different order from the simulation; shield is not measured. Bath has an invalid reading.
    return xr.Dataset(
        {
            "temperature": (
                ("time", "node"),
                np.column_stack([[290.0, 240.0, np.nan, 160.0], [295.0, 140.0, 25.0, 12.0]]),
                {"units": "K"},
            ),
            "scored": ("node", [False, True]),
        },
        coords={"time": ("time", [0.0, 2400.0, 4800.0, 7200.0], {"units": "s"}), "node": ["bath", "stage"]},
    )


@pytest.fixture
def figure():
    figure = plot_measured_temperature(_simulated(), _measured())
    yield figure
    plt.close(figure)


def _traces(ax, linestyle):
    return [line for line in ax.get_lines() if line.get_linestyle() == linestyle]


def test_draws_one_panel_per_measured_node_in_measured_order(figure):
    assert [ax.get_ylabel() for ax in figure.axes] == ["bath (K)", "stage (K)"]


def test_draws_simulated_node_solid_and_measured_trace_dotted_against_hours(figure):
    bath, stage = figure.axes

    (simulated,) = _traces(stage, "-")
    np.testing.assert_allclose(simulated.get_xdata(), [0.0, 0.5, 1.0, 1.0, 1.5, 2.0])
    np.testing.assert_allclose(simulated.get_ydata(), [300.0, 150.0, 70.0, 70.0, 20.0, 10.0])
    (measured,) = _traces(stage, ":")
    np.testing.assert_allclose(measured.get_xdata(), [0.0, 0.6667, 1.3333, 2.0], atol=1e-4)
    np.testing.assert_allclose(measured.get_ydata(), [295.0, 140.0, 25.0, 12.0])

    (simulated,) = _traces(bath, "-")
    np.testing.assert_allclose(simulated.get_ydata(), [300.0, 280.0, 250.0, 250.0, 200.0, 150.0])


def test_leaves_a_gap_at_invalid_readings(figure):
    bath = figure.axes[0]

    (measured,) = _traces(bath, ":")
    np.testing.assert_array_equal(np.isnan(measured.get_ydata()), [False, False, True, False])


def test_measured_trace_shares_its_node_colour_and_nodes_differ(figure):
    bath, stage = figure.axes

    for ax in (bath, stage):
        (simulated,) = _traces(ax, "-")
        (measured,) = _traces(ax, ":")
        assert simulated.get_color() == measured.get_color()
    assert _traces(bath, "-")[0].get_color() != _traces(stage, "-")[0].get_color()


def test_panels_share_time_in_hours_labelled_on_the_bottom(figure):
    bath, stage = figure.axes

    assert bath.get_shared_x_axes().joined(bath, stage)
    assert stage.get_xlabel() == "Time (h)"
    assert bath.get_xlabel() == ""


def test_keys_line_styles_on_the_top_panel_only(figure):
    bath, stage = figure.axes

    assert [text.get_text() for text in bath.get_legend().get_texts()] == ["Simulated", "Measured"]
    assert stage.get_legend() is None


def test_marks_phase_switch_on_every_panel_and_names_it_on_the_top(figure):
    bath, stage = figure.axes

    for ax in (bath, stage):
        (marker,) = _traces(ax, "--")
        np.testing.assert_allclose(marker.get_xdata(), 1.0)
    assert any("precool" in text.get_text() and "pushing" in text.get_text() for text in bath.texts)
    assert not stage.texts


def test_temperature_axes_are_logarithmic_by_default(figure):
    assert {ax.get_yscale() for ax in figure.axes} == {"log"}


def test_temperature_axes_are_linear_on_request():
    figure = plot_measured_temperature(_simulated(), _measured(), log=False)
    try:
        assert {ax.get_yscale() for ax in figure.axes} == {"linear"}
    finally:
        plt.close(figure)


def test_single_measured_node_draws_one_panel():
    figure = plot_measured_temperature(_simulated(), _measured().sel(node=["stage"]))
    try:
        (ax,) = figure.axes
        assert ax.get_ylabel() == "stage (K)"
        assert ax.get_xlabel() == "Time (h)"
        assert [text.get_text() for text in ax.get_legend().get_texts()] == ["Simulated", "Measured"]
    finally:
        plt.close(figure)


def test_node_colours_match_the_temperature_plot(figure):
    fig, ax = plt.subplots()
    try:
        plot_temperature(_simulated(), ax=ax)
        expected = {line.get_label(): to_hex(line.get_color()) for line in ax.get_lines()}
    finally:
        plt.close(fig)

    for panel, node in zip(figure.axes, ["bath", "stage"], strict=True):
        assert to_hex(_traces(panel, ":")[0].get_color()) == expected[node]
