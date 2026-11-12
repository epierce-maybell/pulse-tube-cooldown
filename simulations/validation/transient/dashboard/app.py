"""Streamlit app comparing transient validation runs with their measured cooldowns."""

import math
from pathlib import Path
from typing import Any

import pandas as pd
import streamlit as st

from simulations.validation.transient.dashboard import figures
from simulations.validation.transient.dashboard.runs import (
    Measured,
    Run,
    config_differences,
    find_runs,
    load_measured,
    load_run,
    named_tables,
    plotted_nodes,
)

RUNS = Path(__file__).parents[1]
ROOT = Path(__file__).parents[4]


@st.cache_data
def cached_run(directory: Path, modified: float) -> Run:
    """Loads a run, reloading when `modified` changes."""
    return load_run(directory)


@st.cache_data
def cached_measured(validation: dict) -> Measured:
    """Loads the measured cooldown and legacy baseline named by a run's `[validation]` table."""
    return load_measured({"validation": validation}, ROOT)


def open_run(system: str, label: str) -> Run:
    """Loads `system`'s run `label`, reloading it after the runner rewrites it."""
    directory = RUNS / system / label
    modified = max((directory / name).stat().st_mtime for name in ("result.nc", "metrics.toml"))
    return cached_run(directory, modified)


def hours(seconds: float, *, signed: bool = False) -> str:
    """Formats a time in hours, or "Not reached" for NaN."""
    if math.isnan(seconds):
        return "Not reached"
    return f"{seconds / 3600.0:{'+' if signed else ''}.1f} h"


def scorecard(runs: list[Run], measured: Measured) -> None:
    """Shows the first run's cooldown time error and log-temperature RMS beside the second run's or the legacy's.

    The cooldown time runs from compressor start to the start of condensing.
    """
    run = runs[0]
    other = runs[1] if len(runs) > 1 else None
    other_label = other.label if other else "Legacy"
    other_metrics = other.metrics["model"] if other else run.metrics.get("legacy")

    with st.grid(3, gap="medium"):
        st.metric(
            "Cooldown time error",
            hours(run.metrics["model"]["condensing_error"], signed=True),
            delta=other_metrics and f"{other_label} {hours(other_metrics['condensing_error'], signed=True)}",
            delta_color="off",
            delta_arrow="off",
            border=True,
            height="stretch",
        )
        st.metric(
            "Log-temperature RMS",
            f"{run.metrics['model']['rms']:.2f}",
            delta=other_metrics and f"{other_label} {other_metrics['rms']:.2f}",
            delta_color="off",
            delta_arrow="off",
            border=True,
            height="stretch",
        )
        st.metric("Measured cooldown time", hours(measured.condensing_time), border=True, height="stretch")


def differences(a: Run, b: Run) -> None:
    """Shows the configuration values that differ between two runs."""
    rows = list(config_differences(a.config, b.config))
    with st.expander(f"{len(rows)} configuration differences", expanded=bool(rows)):
        table = pd.DataFrame(
            [(path, describe(value_a), describe(value_b)) for path, value_a, value_b in rows],
            columns=["Parameter", a.label, b.label],
        )
        st.dataframe(
            table, hide_index=True, width="stretch", alt=f"Configuration differences between {a.label} and {b.label}"
        )


def describe(value: Any) -> str:
    """Returns a configuration value as text, naming a whole table only by whether it is present."""
    if value is None:
        return "absent"
    if isinstance(value, dict) or named_tables(value):
        return "present"
    return str(value)


def heat_flow(runs: list[Run], measured: Measured) -> None:
    """Shows one node's temperature and heat inflow by source for each run, one run below the other."""
    plotted = plotted_nodes(runs[0], measured)
    ordered = plotted + [str(node) for node in runs[0].data.node.values if node not in plotted]
    node = st.selectbox("Node", ordered, key="node", bind="query-params")
    for run in runs:
        if len(runs) > 1:
            st.subheader(run.label, divider=False)
        st.plotly_chart(
            figures.heat_flow(run, measured, node),
            config={"displaylogo": False},
            key=f"heat flow {run.label}",
            alt=f"Temperature of and heat flows into {node} in {run.label}",
        )


def network(system: str, runs: list[Run]) -> None:
    """Shows each run's thermal network side by side, drawn for the app's light or dark theme."""
    theme = "dark" if st.context.theme.type == "dark" else "light"
    for column, run in zip(st.columns(len(runs)), runs, strict=True):
        with column:
            if len(runs) > 1:
                st.subheader(run.label, divider=False)
            st.image(RUNS / system / run.label / f"graph_{theme}.svg", alt=f"Thermal network of {run.label}")


def main() -> None:
    """Lays out the selected runs' scorecard, configuration differences, and comparison figures."""
    st.set_page_config(page_title="Cal/val", page_icon=":material/ac_unit:", layout="wide")
    runs_by_system = find_runs(RUNS)
    if not runs_by_system:
        st.info(f"No compared runs under {RUNS}.")
        return

    with st.sidebar:
        system = st.segmented_control(
            "System",
            list(runs_by_system),
            default=next(iter(runs_by_system)),
            format_func=str.upper,
            required=True,
            key="system",
            bind="query-params",
        )
        labels = runs_by_system[system]
        label = st.selectbox("Run", labels, key="run", bind="query-params")
        compared = st.selectbox(
            "Compare with",
            [None, *[other for other in labels if other != label]],
            format_func=lambda other: "Nothing" if other is None else other,
            key="compare",
            bind="query-params",
        )

    runs = [open_run(system, label)]
    if compared:
        runs.append(open_run(system, compared))
    measured = cached_measured(runs[0].config["validation"])

    st.title(f"{system.upper()} · {' vs '.join(run.label for run in runs)}")
    scorecard(runs, measured)
    if compared:
        differences(*runs)

    with st.container(horizontal=True):
        view = st.segmented_control(
            "View",
            ["Temperature", "Residuals", "Heat flow", "Network"],
            default="Temperature",
            required=True,
            label_visibility="collapsed",
            key="view",
            bind="query-params",
        )
        if view == "Temperature":
            layout = st.segmented_control(
                "Nodes",
                ["By node", "All nodes"],
                default="By node",
                required=True,
                label_visibility="collapsed",
                key="nodes",
                bind="query-params",
            )
    if view == "Temperature" and layout == "All nodes":
        st.plotly_chart(
            figures.temperature_overlay(runs[0], measured, plotted_nodes(runs[0], measured)),
            config={"displaylogo": False},
            alt=f"Simulated and measured temperatures of every plotted node in {runs[0].label}",
        )
    elif view == "Temperature":
        st.plotly_chart(
            figures.temperature_by_node(runs, measured, plotted_nodes(runs[0], measured)),
            config={"displaylogo": False},
            alt="Simulated, measured, and legacy temperatures for each plotted node",
        )
    elif view == "Residuals":
        st.plotly_chart(
            figures.residual_by_node(runs, measured),
            config={"displaylogo": False},
            alt="Log-temperature residuals against the measurement for each measured node",
        )
    elif view == "Heat flow":
        heat_flow(runs, measured)
    else:
        network(system, runs)


main()
