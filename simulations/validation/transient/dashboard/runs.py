"""Transient validation runs and the measured and legacy traces they are compared with."""

import dataclasses
import tomllib
from collections.abc import Iterator
from pathlib import Path
from typing import Any

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


@dataclasses.dataclass(frozen=True)
class Measured:
    """A run's measured cooldown and legacy baseline, with nodes named as in the model.

    Attributes:
        cooldown: The measured temperature(time, node) dataset with its scored(node) mask.
        condensing_time: The measured condensing time in s since compressor start.
        legacy: The legacy baseline's temperature(time, node) over the measured nodes, if the run names one.
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
        baseline = rename_nodes(xr.load_dataset(root / validation["baseline"]).temperature, names)
        legacy = baseline.reindex(node=cooldown.node.values)
    return Measured(cooldown, condensing, legacy)


def residual(temperature: xr.DataArray, measured: Measured) -> xr.DataArray:
    """Returns ln(T / T_measured)(time, node) at the measured samples."""
    return log_temperature_ratio(temperature.reindex(node=measured.cooldown.node.values), measured.cooldown)


def config_differences(a: dict[str, Any], b: dict[str, Any]) -> list[tuple[str, Any, Any]]:
    """Returns each configuration value that differs between `a` and `b`, as its path and both values.

    Named tables in arrays, such as nodes and links, are matched by name. A value missing from one configuration is
    None there.
    """
    flat_a, flat_b = dict(flatten(a)), dict(flatten(b))
    paths = list(flat_a) + [path for path in flat_b if path not in flat_a]
    return [(path, flat_a.get(path), flat_b.get(path)) for path in paths if flat_a.get(path) != flat_b.get(path)]


def flatten(value: Any, path: str = "") -> Iterator[tuple[str, Any]]:
    """Yields each leaf of a TOML value with its dotted path. Named tables in arrays are keyed by name."""
    if isinstance(value, dict):
        for key, item in value.items():
            yield from flatten(item, f"{path}.{key}" if path else key)
    elif isinstance(value, list) and value and all(isinstance(item, dict) and "name" in item for item in value):
        for item in value:
            yield from flatten({key: entry for key, entry in item.items() if key != "name"}, f"{path}[{item['name']}]")
    else:
        yield path, value
