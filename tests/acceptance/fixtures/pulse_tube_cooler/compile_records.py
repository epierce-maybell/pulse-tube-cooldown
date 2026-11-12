"""Compile test-stand conditions and objectives with absolute tolerances. Units are s, K, W, and W/K."""

import os
import tomllib
from collections import defaultdict
from collections.abc import Sequence
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

from pulse_tube_cooldown.measurements import read_cooldown, read_load_map

__all__ = ["CooldownObjective", "CouplingObjective", "LoadMapObjective", "Objectives", "compile_records"]


@dataclass(frozen=True, slots=True, kw_only=True)
class CooldownObjective:
    """Crossing time in s at a stage's temperature gate in K."""

    stage: int
    temperature: float
    target: float
    tolerance: float


@dataclass(frozen=True, slots=True, kw_only=True)
class LoadMapObjective:
    """Stage temperature in K at stage-1 and stage-2 heater powers in W."""

    stage: int
    heater_power: tuple[float, float]
    target: float
    tolerance: float


@dataclass(frozen=True, slots=True, kw_only=True)
class CouplingObjective:
    """Net non-heater heat-flow response in W/K, positive into the receiving stage."""

    stage: int
    temperature_stage: int
    target: float
    tolerance: float


@dataclass(frozen=True, slots=True, kw_only=True)
class Objectives:
    """Objectives and mean initial TPR room temperature in K for one compressor count."""

    room_temperature: float
    cooldown: tuple[CooldownObjective, ...]
    load_map: tuple[LoadMapObjective, ...]
    coupling: tuple[CouplingObjective, ...]


def compile_records(
    manifest: str | os.PathLike[str] | None = None,
) -> dict[int, Objectives]:
    """Return objectives keyed by compressor count, using the manifest's records and tolerance rules."""
    path = Path(manifest) if manifest is not None else Path(__file__).with_name("manifest.toml")
    with path.open("rb") as file:
        records = tomllib.load(file)
    settings = records["objectives"]
    powers = tuple(sorted((float(p[0]), float(p[1])) for p in settings["heater_powers"]))
    gates = tuple(tuple(float(t) for t in settings[f"stage_{stage}_gates"]) for stage in (1, 2))
    deviations = float(settings["standard_deviations"])
    relative_tolerance = float(settings["single_record_relative_tolerance"])

    maps: dict[int, list[NDArray[np.float64]]] = defaultdict(list)
    traces: dict[int, list[NDArray[np.float64]]] = defaultdict(list)
    room_temperatures: dict[int, list[float]] = defaultdict(list)
    for entry in records["load_maps"]:
        maps[int(entry["compressors"])].append(_load_temperatures(path.parent / entry["file"], powers))
    for entry in records["cooldowns"]:
        room_temperature, crossing_times = _load_cooldown(
            path.parent / entry["file"],
            gates,
            start=entry["compressor_start"],
            end=entry["end"],
            heater_start=entry.get("heater_start"),
            stage_1_max_temperature=float(settings["stage_1_max_temperature"]),
        )
        count = int(entry["compressors"])
        room_temperatures[count].append(room_temperature)
        traces[count].append(crossing_times)

    coupling = {count: np.array([_coupling(t, powers) for t in temperatures]) for count, temperatures in maps.items()}
    coupling_tolerance = deviations * coupling[1].std(axis=0, ddof=1)
    objectives = {}
    for count, temperatures in maps.items():
        map_target, map_tolerance = _summarize(
            temperatures, deviations=deviations, relative_tolerance=relative_tolerance
        )
        time_target, time_tolerance = _summarize(
            traces[count],
            deviations=float(settings["timing_standard_deviations"]),
            relative_tolerance=relative_tolerance,
        )
        coupling_target = coupling[count].mean(axis=0)
        objectives[count] = Objectives(
            room_temperature=float(np.mean(room_temperatures[count])),
            cooldown=tuple(
                CooldownObjective(
                    stage=stage + 1,
                    temperature=gate,
                    target=float(time_target[stage, index]),
                    tolerance=float(time_tolerance[stage, index]),
                )
                for stage, thresholds in enumerate(gates)
                for index, gate in enumerate(thresholds)
            ),
            load_map=tuple(
                LoadMapObjective(
                    stage=stage + 1,
                    heater_power=power,
                    target=float(map_target[index, stage]),
                    tolerance=float(map_tolerance[index, stage]),
                )
                for index, power in enumerate(powers)
                for stage in range(2)
            ),
            coupling=tuple(
                CouplingObjective(
                    stage=stage + 1,
                    temperature_stage=other + 1,
                    target=float(coupling_target[stage, other]),
                    tolerance=float(coupling_tolerance[stage, other]),
                )
                for stage in range(2)
                for other in range(2)
            ),
        )
    return objectives


def _load_temperatures(path: Path, powers: Sequence[tuple[float, float]]) -> NDArray[np.float64]:
    points = read_load_map(path)
    readings = np.column_stack((points.T1.values, points.T2.values))
    temperatures = []
    for first, second in powers:
        selected = readings[(points.Q1.values == first) & (points.Q2.values == second)]
        if not len(selected) or not np.isfinite(selected).all():
            raise ValueError(f"{path.name}: missing or nonfinite temperatures at ({first}, {second}) W.")
        temperatures.append(selected.mean(axis=0))
    return np.array(temperatures)


def _load_cooldown(
    path: Path,
    gates: Sequence[Sequence[float]],
    *,
    start: datetime,
    end: datetime,
    heater_start: datetime | None,
    stage_1_max_temperature: float,
) -> tuple[float, NDArray[np.float64]]:
    """Return the last finite prestart TPR reading in K and interpolated downward crossing times in s."""
    trace = read_cooldown(path)
    prestart = trace.TPR.values[(trace.time.values < np.datetime64(start)) & np.isfinite(trace.TPR.values)]
    if not prestart.size:
        raise ValueError(f"{path.name}: no finite prestart TPR reading.")
    room_temperature = float(prestart[-1])
    trace = trace.sel(time=slice(np.datetime64(start), np.datetime64(end)))
    if heater_start is not None:
        trace = trace.where(trace.time < np.datetime64(heater_start), drop=True)
    seconds = (trace.time.values - np.datetime64(start)) / np.timedelta64(1, "s")
    times = []
    for stage, thresholds in enumerate(gates, start=1):
        readings = trace[f"T{stage}"].values
        valid = np.isfinite(readings)
        valid &= readings < stage_1_max_temperature if stage == 1 else readings <= room_temperature
        retained = readings[valid]
        retained_seconds = seconds[valid]
        crossings = []
        for gate in thresholds:
            downward = np.flatnonzero((retained[:-1] > gate) & (retained[1:] <= gate))
            if not downward.size:
                raise ValueError(f"{path.name}: stage {stage} never crosses {gate} K before cutoff.")
            index = downward[0]
            crossings.append(
                float(
                    np.interp(
                        gate,
                        retained[index : index + 2][::-1],
                        retained_seconds[index : index + 2][::-1],
                    )
                )
            )
        times.append(crossings)
    return room_temperature, np.array(times)


def _coupling(temperatures: NDArray[np.float64], powers: Sequence[tuple[float, float]]) -> NDArray[np.float64]:
    """Return dQ/dT from corners ordered (low, low), (low, high), (high, low), (high, high)."""
    low_low, low_high, high_low, high_high = temperatures
    response = np.column_stack(
        ((high_low - low_low + high_high - low_high) / 2, (low_high - low_low + high_high - high_low) / 2)
    )
    return -np.diag(np.ptp(powers, axis=0)) @ np.linalg.inv(response)


def _summarize(
    samples: Sequence[NDArray[np.float64]],
    *,
    deviations: float,
    relative_tolerance: float,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    values = np.array(samples)
    target = values.mean(axis=0)
    tolerance = deviations * values.std(axis=0, ddof=1) if len(values) > 1 else relative_tolerance * np.abs(target)
    return target, tolerance
