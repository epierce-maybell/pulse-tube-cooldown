"""Plotly figures of transient validation runs against their measured cooldowns.

Times are plotted in hours since compressor start. Each run keeps one color across figures. Series colors are
mid-lightness hues that read on light and dark backgrounds; gridlines, axes, and text follow the app theme.
"""

import math

import numpy as np
import plotly.graph_objects as go
import xarray as xr
from plotly.subplots import make_subplots

from simulations.validation.transient.dashboard.runs import Measured, Run, heat_sources, residual

SECONDS_PER_HOUR = 3600.0
COLUMNS = 3
HEAT_DECADES = 5
HEAT_FLOW_HEIGHT = 760
HEAT_FLOW_LEGEND_WIDTH = 280
LOG_TICKS = [scale * decade for decade in (1.0, 10.0, 100.0) for scale in (1, 2, 5)] + [1000.0]
TIME_TITLE = "Time since compressor start (h)"
PANEL_HEIGHT = 240
ROW_GAP = 100
PANELS_TOP = 76
SINGLE_TOP = 48
BOTTOM = 8

CATEGORICAL = ("#3987e5", "#d95926", "#199e70", "#c98500", "#d55181", "#008300", "#9085e9", "#e66767")
RUN_COLORS = ("#3987e5", "#d95926")
MEASURED_COLOR = "#199e70"
LEGACY_COLOR = "#9085e9"
NET_COLOR = "#e66767"
NEUTRAL = "#898781"


def temperature_by_node(runs: list[Run], measured: Measured, nodes: list[str]) -> go.Figure:
    """Plots each node's temperature in its own panel, with its measurement and legacy baseline where they exist."""
    figure = node_panels(nodes, shared_y=False)
    for i, node in enumerate(nodes):
        row, col = panel(i)
        if node in measured.nodes:
            band = measured_trace(measured.cooldown.temperature.sel(node=node), MEASURED_COLOR, "Measured")
            figure.add_trace(band, row=row, col=col)
        if measured.legacy is not None and node in measured.legacy.node:
            figure.add_trace(legacy_trace(measured.legacy.sel(node=node)), row=row, col=col)
        for run, color in zip(runs, RUN_COLORS, strict=False):
            figure.add_trace(line(run.data.temperature.sel(node=node), color, run.label), row=row, col=col)
    one_legend_entry_per_series(figure)
    figure.update_yaxes(type="log", tickvals=LOG_TICKS)
    label_panels(figure, nodes, "Temperature (K)")
    return style(figure, height=panels_height(len(nodes)), top=PANELS_TOP)


def temperature_overlay(run: Run, measured: Measured, nodes: list[str]) -> go.Figure:
    """Plots every node's simulated and any measured temperatures in one panel, one color per node."""
    figure = go.Figure()
    for node, color in zip(nodes, CATEGORICAL, strict=False):
        if node in measured.nodes:
            band = measured_trace(measured.cooldown.temperature.sel(node=node), color, f"{node} measured", node)
            figure.add_trace(band.update(showlegend=False))
        figure.add_trace(line(run.data.temperature.sel(node=node), color, node))
    figure.update_yaxes(type="log", tickvals=LOG_TICKS, title_text="Temperature (K)")
    figure.update_xaxes(title_text=TIME_TITLE)
    return style(figure, height=560, top=SINGLE_TOP)


def residual_by_node(runs: list[Run], measured: Measured) -> go.Figure:
    """Plots each run's and the legacy baseline's ln(T / T_measured) per measured node, on a shared symmetric scale."""
    traces = {run.label: residual(run.data.temperature, measured) for run in runs}
    if measured.legacy is not None:
        legacy = residual(measured.legacy, measured)
    figure = node_panels(measured.nodes, shared_y=True)
    for i, node in enumerate(measured.nodes):
        row, col = panel(i)
        if measured.legacy is not None:
            figure.add_trace(legacy_trace(legacy.sel(node=node), unit=""), row=row, col=col)
        for (label, ratio), color in zip(traces.items(), RUN_COLORS, strict=False):
            figure.add_trace(line(ratio.sel(node=node), color, label, unit=""), row=row, col=col)
        figure.add_hline(y=0.0, line={"color": NEUTRAL, "width": 1}, row=row, col=col)
    one_legend_entry_per_series(figure)
    ratios = [*traces.values(), *([legacy] if measured.legacy is not None else [])]
    extent = max(float(np.nanmax(np.abs(ratio.values))) for ratio in ratios)
    figure.update_yaxes(range=[-1.05 * extent, 1.05 * extent])
    label_panels(figure, measured.nodes, "ln(T / T measured)")
    return style(figure, height=panels_height(len(measured.nodes)), top=PANELS_TOP)


def heat_flow(run: Run, measured: Measured, node: str) -> go.Figure:
    """Plots a node's temperature and net heat inflow above the heat flowing into it from each source.

    Both heat flow axes are signed log scales spanning `HEAT_DECADES` decades below the largest flow.
    """
    sources = heat_sources(run, node)
    net = sources.sum("source")
    scale, ticks, labels = signed_log_axis(float(max(np.abs(sources).max(), np.abs(net).max())))
    figure = make_subplots(
        rows=2,
        cols=1,
        shared_xaxes=True,
        row_heights=[0.4, 0.6],
        vertical_spacing=0.05,
        specs=[[{"secondary_y": True}], [{}]],
    )
    if node in measured.nodes:
        band = measured_trace(measured.cooldown.temperature.sel(node=node), MEASURED_COLOR, "Measured")
        figure.add_trace(band, row=1, col=1)
    figure.add_trace(line(run.data.temperature.sel(node=node), RUN_COLORS[0], run.label), row=1, col=1)
    figure.add_trace(signed_log_line(net, scale, NET_COLOR, "Net heat in"), row=1, col=1, secondary_y=True)
    colors = CATEGORICAL[: sources.sizes["source"]]
    for source, color in zip(map(str, sources.source.values), colors, strict=True):
        trace = signed_log_line(sources.sel(source=source), scale, color, source).update(legend="legend2")
        figure.add_trace(trace, row=2, col=1)
    figure.add_hline(y=0.0, line={"color": NEUTRAL, "width": 1}, row=2, col=1)
    figure.update_yaxes(type="log", tickvals=LOG_TICKS, title_text="Temperature (K)", row=1, col=1, secondary_y=False)
    heat_axis = {"tickvals": ticks, "ticktext": labels}
    figure.update_yaxes(heat_axis, showgrid=False, title_text=f"Net heat into {node}", row=1, col=1, secondary_y=True)
    figure.update_yaxes(heat_axis, title_text=f"Heat into {node}", row=2, col=1)
    figure.update_xaxes(title_text=TIME_TITLE, row=2, col=1)
    style(figure, height=HEAT_FLOW_HEIGHT, top=SINGLE_TOP)
    legend = {"orientation": "v", "xref": "container", "x": 1.0, "xanchor": "right", "yref": "paper", "yanchor": "top"}
    sources_top = figure.get_subplot(2, 1).yaxis.domain[1]
    figure.update_layout(
        margin={"r": HEAT_FLOW_LEGEND_WIDTH},
        legend=legend | {"y": 1.0},
        legend2=legend | {"y": sources_top, "title": f"Heat into {node} from"},
    )
    return figure


def signed_log_axis(largest: float) -> tuple[float, list[float], list[str]]:
    """Returns the linear scale in W, tick positions, and tick labels of a signed log axis reaching `largest` W.

    The axis plots asinh(Q / scale), which is linear within `scale` of zero and logarithmic beyond it. Ticks fall on
    zero and on each power of ten from `HEAT_DECADES` decades below `largest` up to it, with either sign.
    """
    top = math.ceil(math.log10(largest))
    scale = 10.0 ** (top - HEAT_DECADES)
    powers = [10.0**exponent for exponent in range(top - HEAT_DECADES + 1, top + 1)]
    values = [-power for power in reversed(powers)] + [0.0] + powers
    return scale, [float(np.arcsinh(value / scale)) for value in values], [watts(value) for value in values]


def watts(power: float) -> str:
    """Formats zero or a signed power of ten in W with an SI prefix, such as "−10 mW"."""
    if power == 0.0:
        return "0"
    exponent = round(math.log10(abs(power)))
    thousands = 3 * math.floor(exponent / 3)
    prefix = {-9: "n", -6: "µ", -3: "m", 0: "", 3: "k", 6: "M"}[thousands]
    return f"{'−' if power < 0 else '+'}{10 ** (exponent - thousands)} {prefix}W"


def signed_log_line(values: xr.DataArray, scale: float, color: str, name: str) -> go.Scatter:
    """Returns a 2 px line of heat flow over time on a signed log axis of linear scale `scale` W."""
    return go.Scatter(
        x=values.time.values / SECONDS_PER_HOUR,
        y=np.arcsinh(values.values / scale),
        name=name,
        line={"width": 2, "color": color},
        customdata=values.values,
        hovertemplate="%{customdata:.3g} W",
    )


def node_panels(nodes: list[str], *, shared_y: bool) -> go.Figure:
    """Returns a grid with one titled panel per node, with zoom linked across the time axes."""
    rows, cols = math.ceil(len(nodes) / COLUMNS), min(COLUMNS, len(nodes))
    specs = [[{} if row * cols + col < len(nodes) else None for col in range(cols)] for row in range(rows)]
    figure = make_subplots(
        rows=rows,
        cols=cols,
        specs=specs,
        subplot_titles=nodes,
        horizontal_spacing=0.07,
        vertical_spacing=ROW_GAP / plotting_height(rows) if rows > 1 else 0.0,
    )
    figure.update_xaxes(matches="x")
    if shared_y:
        figure.update_yaxes(matches="y")
    return figure


def label_panels(figure: go.Figure, nodes: list[str], y_title: str) -> None:
    """Titles the y axis of each row's first panel and the time axis of each column's lowest panel."""
    for i in range(len(nodes)):
        row, col = panel(i)
        if col == 1:
            figure.update_yaxes(title_text=y_title, row=row, col=col)
        if i + COLUMNS >= len(nodes):
            figure.update_xaxes(title_text=TIME_TITLE, row=row, col=col)


def panel(index: int) -> tuple[int, int]:
    """Returns the 1-based row and column of the panel at `index`."""
    row, col = divmod(index, COLUMNS)
    return row + 1, col + 1


def plotting_height(rows: int) -> int:
    """Returns the height in px of `rows` panel rows and the gaps between them."""
    return rows * PANEL_HEIGHT + (rows - 1) * ROW_GAP


def panels_height(count: int) -> int:
    """Returns the figure height in px for `count` node panels."""
    return plotting_height(math.ceil(count / COLUMNS)) + PANELS_TOP + BOTTOM


def one_legend_entry_per_series(figure: go.Figure) -> None:
    """Shows each series name once in the legend, at its first trace, keeping hidden traces hidden."""
    shown: set[str] = set()
    for trace in figure.data:
        if trace.showlegend is not False:
            trace.showlegend = trace.name not in shown
            shown.add(trace.name)


def line(values: xr.DataArray, color: str, name: str, *, unit: str = " K") -> go.Scatter:
    """Returns a 2 px line over time, grouped in the legend by `name`."""
    return go.Scatter(
        x=values.time.values / SECONDS_PER_HOUR,
        y=values.values,
        name=name,
        legendgroup=name,
        line={"width": 2, "color": color},
        hovertemplate=f"%{{y:.3g}}{unit}",
    )


def measured_trace(values: xr.DataArray, color: str, name: str, legendgroup: str | None = None) -> go.Scatter:
    """Returns a wide translucent band for a measured trace, grouped in the legend by `legendgroup` or `name`."""
    return go.Scatter(
        x=values.time.values / SECONDS_PER_HOUR,
        y=values.values,
        name=name,
        legendgroup=legendgroup or name,
        line={"width": 6, "color": color},
        opacity=0.45,
        hovertemplate="%{y:.3g} K",
    )


def legacy_trace(values: xr.DataArray, *, unit: str = " K") -> go.Scatter:
    """Returns a dotted line for a legacy baseline trace."""
    return go.Scatter(
        x=values.time.values / SECONDS_PER_HOUR,
        y=values.values,
        name="Legacy",
        legendgroup="Legacy",
        line={"width": 2, "color": LEGACY_COLOR, "dash": "dot"},
        hovertemplate=f"%{{y:.3g}}{unit}",
    )


def style(figure: go.Figure, *, height: int, top: int) -> go.Figure:
    """Applies the shared layout: transparent background, legend above, and a 12 h time grid.

    Args:
        figure: The figure to style.
        height: The figure height in px.
        top: The space above the plotting area in px, holding the legend and any panel titles.
    """
    figure.update_layout(
        height=height,
        margin={"l": 8, "r": 8, "t": top, "b": BOTTOM},
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        hovermode="x unified",
        legend={"orientation": "h", "x": 0.0, "y": 1.0, "yanchor": "top", "yref": "container", "title": None},
    )
    figure.update_xaxes(dtick=12, zeroline=False)
    figure.update_yaxes(zeroline=False)
    return figure
