"""Plotly figures of transient validation runs against their measured cooldowns.

Times are plotted in hours since compressor start. Each run keeps one color across figures. Series colors are
mid-lightness hues that read on light and dark backgrounds; gridlines, axes, and text follow the app theme.
"""

import math

import numpy as np
import plotly.graph_objects as go
import xarray as xr
from plotly.subplots import make_subplots

from simulations.validation.transient.dashboard.runs import Measured, Run, residual

SECONDS_PER_HOUR = 3600.0
COLUMNS = 3
OTHER = "Other"
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


def temperature_by_node(runs: list[Run], measured: Measured) -> go.Figure:
    """Plots each measured node's temperature in its own panel, against its measurement and legacy baseline."""
    figure = node_panels(measured.nodes, shared_y=False)
    for i, node in enumerate(measured.nodes):
        row, col = panel(i)
        first = i == 0
        figure.add_trace(
            measured_trace(measured.cooldown.temperature.sel(node=node), MEASURED_COLOR, "Measured", first),
            row=row,
            col=col,
        )
        if measured.legacy is not None:
            figure.add_trace(legacy_trace(measured.legacy.sel(node=node), first), row=row, col=col)
        for run, color in zip(runs, RUN_COLORS, strict=False):
            figure.add_trace(
                line(run.data.temperature.sel(node=node), color, run.label, legendgroup=run.label, showlegend=first),
                row=row,
                col=col,
            )
    figure.update_yaxes(type="log", tickvals=LOG_TICKS)
    label_panels(figure, measured.nodes, "Temperature (K)")
    return style(figure, height=panels_height(len(measured.nodes)), top=PANELS_TOP)


def temperature_overlay(run: Run, measured: Measured) -> go.Figure:
    """Plots every measured node's simulated and measured temperatures in one panel, one color per node."""
    figure = go.Figure()
    for node, color in zip(measured.nodes, CATEGORICAL, strict=False):
        figure.add_trace(measured_trace(measured.cooldown.temperature.sel(node=node), color, node, False))
        figure.add_trace(line(run.data.temperature.sel(node=node), color, node, legendgroup=node, showlegend=True))
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
        first = i == 0
        if measured.legacy is not None:
            figure.add_trace(legacy_trace(legacy.sel(node=node), first), row=row, col=col)
        for (label, ratio), color in zip(traces.items(), RUN_COLORS, strict=False):
            figure.add_trace(
                line(ratio.sel(node=node), color, label, legendgroup=label, showlegend=first, unit=""),
                row=row,
                col=col,
            )
        figure.add_hline(y=0.0, line={"color": NEUTRAL, "width": 1}, row=row, col=col)
    ratios = [*traces.values(), *([legacy] if measured.legacy is not None else [])]
    extent = max(float(np.nanmax(np.abs(ratio.values))) for ratio in ratios)
    figure.update_yaxes(range=[-1.05 * extent, 1.05 * extent])
    label_panels(figure, measured.nodes, "ln(T / T measured)")
    return style(figure, height=panels_height(len(measured.nodes)), top=PANELS_TOP)


def budget(run: Run, node: str, shown: list[str]) -> xr.DataArray:
    """Returns the heat flow(time, link) into `node` from each link in `shown` and from all its other links as one.

    Inactive links carry zero.
    """
    flows = run.data.heat_flow.sel(node=node, link=run.data.connected.sel(node=node)).fillna(0.0)
    links = [str(link) for link in flows.link.values]
    rest = [link for link in links if link not in shown]
    series = flows.sel(link=[link for link in shown if link in links])
    if rest:
        other = flows.sel(link=rest).sum("link").expand_dims(link=[OTHER])
        series = xr.concat([series, other], dim="link")
    return series.transpose("time", "link")


def largest_links(run: Run, node: str, count: int) -> list[str]:
    """Returns the `count` links with the largest mean absolute heat flow into `node`, in descending order."""
    flows = run.data.heat_flow.sel(node=node, link=run.data.connected.sel(node=node)).fillna(0.0)
    size = np.abs(flows).mean("time")
    return [str(link) for link in size.sortby(size, ascending=False).link.values[:count]]


def budget_extent(flows: xr.DataArray) -> tuple[float, float]:
    """Returns the 1st percentile of stacked cooling and 99th percentile of stacked heating in W.

    Percentiles keep a brief start-up transient from setting the scale.
    """
    cooling = flows.clip(max=0.0).sum("link").values
    heating = flows.clip(min=0.0).sum("link").values
    return float(np.percentile(cooling, 1.0)), float(np.percentile(heating, 99.0))


def heat_flow_budget(flows: xr.DataArray, extent: tuple[float, float]) -> go.Figure:
    """Plots a node's heat flows stacked above zero for heating and below for cooling, with their net.

    Args:
        flows: The node's heat flow(time, link) in W, from `budget`.
        extent: The y range in W.
    """
    hours = flows.time.values / SECONDS_PER_HOUR
    colors = dict(zip([str(link) for link in flows.link.values if link != OTHER], CATEGORICAL[:-1], strict=False))
    colors[OTHER] = NEUTRAL
    figure = go.Figure()
    for link in flows.link.values:
        values = flows.sel(link=link).values
        for stack, part in (("heating", np.clip(values, 0.0, None)), ("cooling", np.clip(values, None, 0.0))):
            figure.add_trace(
                go.Scatter(
                    x=hours,
                    y=part,
                    name=str(link),
                    legendgroup=str(link),
                    showlegend=stack == "heating",
                    stackgroup=stack,
                    line={"width": 0.5, "color": colors[str(link)]},
                    fillcolor=colors[str(link)],
                    customdata=values,
                    hovertemplate="%{customdata:.3g} W" if stack == "heating" else None,
                    hoverinfo=None if stack == "heating" else "skip",
                )
            )
    figure.add_trace(
        go.Scatter(
            x=hours,
            y=flows.sum("link").values,
            name="Net",
            line={"width": 2.5, "color": NET_COLOR},
            hovertemplate="%{y:.3g} W",
        )
    )
    low, high = extent
    margin = 0.05 * (high - low)
    figure.update_yaxes(title_text="Heat flow (W)", range=[low - margin, high + margin])
    figure.update_xaxes(title_text=TIME_TITLE)
    figure.add_hline(y=0.0, line={"color": NEUTRAL, "width": 1})
    return style(figure, height=500, top=SINGLE_TOP)


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


def line(
    values: xr.DataArray, color: str, name: str, *, legendgroup: str, showlegend: bool, unit: str = " K"
) -> go.Scatter:
    """Returns a 2 px line over time."""
    return go.Scatter(
        x=values.time.values / SECONDS_PER_HOUR,
        y=values.values,
        name=name,
        legendgroup=legendgroup,
        showlegend=showlegend,
        line={"width": 2, "color": color},
        hovertemplate=f"%{{y:.3g}}{unit}",
    )


def measured_trace(values: xr.DataArray, color: str, name: str, showlegend: bool) -> go.Scatter:
    """Returns a wide translucent band for a measured trace."""
    return go.Scatter(
        x=values.time.values / SECONDS_PER_HOUR,
        y=values.values,
        name=name,
        legendgroup=name if showlegend else None,
        showlegend=showlegend,
        line={"width": 6, "color": color},
        opacity=0.45,
        hoverinfo="skip",
    )


def legacy_trace(values: xr.DataArray, showlegend: bool) -> go.Scatter:
    """Returns a dotted line for a legacy baseline trace."""
    return go.Scatter(
        x=values.time.values / SECONDS_PER_HOUR,
        y=values.values,
        name="Legacy",
        legendgroup="Legacy",
        showlegend=showlegend,
        line={"width": 2, "color": LEGACY_COLOR, "dash": "dot"},
        hovertemplate="%{y:.3g}",
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
