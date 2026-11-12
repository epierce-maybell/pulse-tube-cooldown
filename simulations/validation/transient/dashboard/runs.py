"""Transient validation runs and the measured and legacy traces they are compared with."""

import dataclasses
import tomllib
from collections.abc import Iterator
from pathlib import Path
from typing import Any

import numpy as np
import xarray as xr

from pulse_tube_cooldown.measurements import log_temperature_ratio, read_refrigerator_cooldown
from simulations.validation.transient.run import rename_nodes


@dataclasses.dataclass(frozen=True)
class Run:
    """A simulated run with its configuration and the runner's metrics.

    Attributes:
        label: The run directory's name.
        config: The run's configuration.
        data: The run's Result data.
        metrics: The runner's metrics, with times in s since compressor start.
    """

    label: str
    config: dict[str, Any]
    data: xr.Dataset
    metrics: dict[str, Any]

    @property
    def cooldown_node(self) -> str:
        """The node whose temperature ends the procedure, and so defines the cooldown time."""
        return self.config["procedure"]["until"]["node"]


@dataclasses.dataclass(frozen=True)
class Measured:
    """A run's measured cooldown and legacy baseline, with nodes named as in the model.

    Attributes:
        cooldown: The measured temperature(time, node) dataset with its scored(node) mask.
        condensing_time: The measured condensing time in s since compressor start.
        legacy: The legacy baseline's temperature(time, node), if the run names one.
    """

    cooldown: xr.Dataset
    condensing_time: float
    legacy: xr.DataArray | None

    @property
    def nodes(self) -> list[str]:
        """The measured nodes."""
        return [str(node) for node in self.cooldown.node.values]


def find_runs(directory: Path) -> dict[str, list[str]]:
    """Returns the labels of the compared runs under `directory`, by system.

    A compared run is a `<system>/<label>/` directory holding `result.nc` and `metrics.toml`.
    """
    runs: dict[str, list[str]] = {}
    for metrics in sorted(directory.glob("*/*/metrics.toml")):
        run = metrics.parent
        if (run / "result.nc").exists():
            runs.setdefault(run.parent.name, []).append(run.name)
    return runs


def load_run(directory: Path) -> Run:
    """Loads the run in `directory`."""
    with (directory / "config.toml").open("rb") as file:
        config = tomllib.load(file)
    with (directory / "metrics.toml").open("rb") as file:
        metrics = tomllib.load(file)
    return Run(directory.name, config, xr.load_dataset(directory / "result.nc"), metrics)


def load_measured(config: dict[str, Any], root: Path) -> Measured:
    """Loads the measured cooldown and legacy baseline named by a run's `[validation]` table.

    Args:
        config: The run's configuration.
        root: The repository root, which the `[validation]` paths are relative to.
    """
    validation = config["validation"]
    names = validation.get("comparison_nodes", {})
    conventions_path = root / validation["measured"]
    with conventions_path.open("rb") as file:
        conventions = tomllib.load(file)
    cooldown = rename_nodes(read_refrigerator_cooldown(conventions_path), names)
    condensing = (conventions["condensing_begins"] - conventions["compressor_start"]).total_seconds()
    legacy = None
    if "baseline" in validation:
        legacy = rename_nodes(xr.load_dataset(root / validation["baseline"]).temperature, names)
    return Measured(cooldown, condensing, legacy)


def plotted_nodes(run: Run, measured: Measured) -> list[str]:
    """Returns the measured nodes followed by the run's cooldown node if it is unmeasured."""
    nodes = measured.nodes
    if run.cooldown_node not in nodes:
        nodes.append(run.cooldown_node)
    return nodes


def residual(temperature: xr.DataArray, measured: Measured) -> xr.DataArray:
    """Returns ln(T / T_measured)(time, node) at the measured samples."""
    return log_temperature_ratio(temperature.reindex(node=measured.cooldown.node.values), measured.cooldown)


def heat_sources(run: Run, node: str) -> xr.DataArray:
    """Returns the heat flow(time, source) into `node` in W, largest mean magnitude first.

    A link between `node` and one other node belongs to the source named for that node. A link to several nodes is
    its own source, and such links that carry the same flow into `node` are one source named "A + B". Inactive links
    carry zero.
    """
    flows = run.data.heat_flow.sel(node=node, link=run.data.connected.sel(node=node)).fillna(0.0)
    sources: dict[str, list[str]] = {}
    shared: list[list[str]] = []
    for link in map(str, flows.link.values):
        connected = run.data.connected.sel(link=link)
        neighbors = [str(other) for other in connected.node.values[connected.values] if other != node]
        if len(neighbors) == 1:
            sources.setdefault(neighbors[0], []).append(link)
            continue
        match = next((group for group in shared if np.allclose(flows.sel(link=group[0]), flows.sel(link=link))), None)
        if match:
            match.append(link)
        else:
            shared.append([link])
    sources |= {" + ".join(group): group for group in shared}
    data = xr.concat([flows.sel(link=links).sum("link") for links in sources.values()], dim="source")
    data = data.assign_coords(source=list(sources))
    size = np.abs(data).mean("time")
    return data.sortby(size, ascending=False)


def config_differences(a: Any, b: Any, path: str = "") -> Iterator[tuple[str, Any, Any]]:
    """Yields each configuration value that differs between `a` and `b`, as its dotted path and both values.

    Tables present in both are compared value by value. Named tables in arrays, such as nodes and links, are matched
    by name. A value or table present in only one configuration is one difference, with None on the other side.
    """
    if isinstance(a, dict) and isinstance(b, dict):
        for key in [*a, *(key for key in b if key not in a)]:
            yield from config_differences(a.get(key), b.get(key), f"{path}.{key}" if path else key)
    elif named_tables(a) and named_tables(b):
        by_name_a, by_name_b = {item["name"]: item for item in a}, {item["name"]: item for item in b}
        for name in [*by_name_a, *(name for name in by_name_b if name not in by_name_a)]:
            yield from config_differences(by_name_a.get(name), by_name_b.get(name), f"{path}[{name}]")
    elif a != b:
        yield path, a, b


def named_tables(value: Any) -> bool:
    """Returns whether `value` is a non-empty array of tables that each have a name."""
    return isinstance(value, list) and bool(value) and all(isinstance(item, dict) and "name" in item for item in value)
