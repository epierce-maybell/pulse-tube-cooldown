"""Parameters, coefficients, and calibrations of the pulse-tube head model, and the files that hold them.

The package's `data/` directory holds the head designs, the parameters shared by every head, and the fit parameters.
A calibration file holds what belongs to one head and drive condition: its fitted coefficients, measured load map,
test stand, and the name of its packaged design.
"""

import os
import tomllib
from collections.abc import Mapping
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Self

import numpy as np
import tomli_w
import xarray as xr
from numpy.typing import NDArray

from pulse_tube_cooldown.core.types import Material
from pulse_tube_cooldown.library.materials.load import load_material
from pulse_tube_cooldown.library.toml_files import read_named_toml

__all__ = [
    "ClampResponseParameters",
    "FitParameters",
    "HeadCalibration",
    "HeadCoefficients",
    "HeadDesignParameters",
    "LoadMap",
    "LoadedProfileParameters",
    "RegeneratorLayer",
    "SharedParameters",
    "UnloadedProfileParameters",
    "load_fit_parameters",
    "load_head_design",
    "load_shared_parameters",
    "read_head_calibration",
    "save_head_calibration",
]

DATA_DIRECTORY = Path(__file__).parent / "data"


@dataclass(frozen=True, slots=True, kw_only=True)
class HeadCoefficients:
    """The seven coefficients fitted for one head and drive condition."""

    stage_1_warm_amplitude: float  # [W]
    stage_1_warm_exponent: float  # unitless
    stage_2_warm_amplitude: float  # [W]
    stage_2_warm_exponent: float  # unitless
    regenerator_1_exchange_coefficient: float  # [W/(kg K)]
    regenerator_2_exchange_coefficient: float  # [W/(kg K)]
    energy_coupling_coefficient: float  # unitless


@dataclass(frozen=True, slots=True, kw_only=True)
class RegeneratorLayer:
    """One material fill of a regenerator.

    Attributes:
        material: The fill material.
        mass: Mass in kg.
        length_fraction: Fraction of the regenerator's length.
    """

    material: Material
    mass: float
    length_fraction: float


@dataclass(frozen=True, slots=True, kw_only=True)
class HeadDesignParameters:
    """A head design's regenerator inventories and clamp position.

    Attributes:
        name: Design name, such as "PT1.5".
        clamp_fraction: Clamp position as a fraction of regenerator 2's length from its warm end.
        regenerator_1: Layers from warm to cold.
        regenerator_2: Layers from warm to cold.
    """

    name: str
    clamp_fraction: float
    regenerator_1: tuple[RegeneratorLayer, ...]
    regenerator_2: tuple[RegeneratorLayer, ...]


@dataclass(frozen=True, slots=True, kw_only=True)
class UnloadedProfileParameters:
    """The unloaded regenerator temperature profile.

    Attributes:
        position: Positions as fractions of the regenerator length from its warm end.
        temperature: Temperatures at `position`, in K.
        reference_stage_temperature: The profile's warm-end and cold-end temperatures, in K.
    """

    position: NDArray[np.float64]
    temperature: NDArray[np.float64]
    reference_stage_temperature: tuple[float, float]


@dataclass(frozen=True, slots=True, kw_only=True)
class ClampResponseParameters:
    """Stage and clamp temperatures against clamp load, at fixed stage heater powers.

    Attributes:
        load: Clamp loads in W.
        stage_1_temperature: Stage-1 temperatures at `load`, in K.
        stage_2_temperature: Stage-2 temperatures at `load`, in K.
        clamp_temperature: Clamp temperatures at `load`, in K.
    """

    load: NDArray[np.float64]
    stage_1_temperature: NDArray[np.float64]
    stage_2_temperature: NDArray[np.float64]
    clamp_temperature: NDArray[np.float64]


@dataclass(frozen=True, slots=True, kw_only=True)
class LoadedProfileParameters:
    """Regenerator-2 temperature profiles at several clamp loads.

    Attributes:
        load: Clamp loads in W.
        position: Positions as fractions of the regenerator length from its warm end.
        temperature: Temperatures in K, shape (loads, positions).
        source_clamp_fraction: Clamp position of the source profiles, as a fraction of the regenerator length.
        reference_stage_temperature: The profiles' warm-end and cold-end temperatures, in K.
    """

    load: NDArray[np.float64]
    position: NDArray[np.float64]
    temperature: NDArray[np.float64]
    source_clamp_fraction: float
    reference_stage_temperature: tuple[float, float]


@dataclass(frozen=True, slots=True, kw_only=True)
class SharedParameters:
    """Parameters every head shares.

    Attributes:
        energy_scale: Normalization of excess regenerator-1 energy, in J.
        cold_blend_width: Width of the warm blend, in temperature coordinates normalized by the measured spans.
        cold_tail_span_multiplier: Continuation tail lengths as multiples of the measured spans, in stage order.
        energy_temperature_range: Lowest and highest temperatures of the tabulated specific energies, in K.
        energy_knots: Number of log-spaced temperatures in the tabulated specific energies.
        quadrature_order: Gauss-Legendre order for integrating target energies over each regenerator node.
        unloaded_profile: The unloaded regenerator profile.
        clamp_response: The clamp-load temperature response.
        loaded_profiles: The loaded regenerator-2 profiles.
    """

    energy_scale: float
    cold_blend_width: float
    cold_tail_span_multiplier: tuple[float, float]
    energy_temperature_range: tuple[float, float]
    energy_knots: int
    quadrature_order: int
    unloaded_profile: UnloadedProfileParameters
    clamp_response: ClampResponseParameters
    loaded_profiles: LoadedProfileParameters


@dataclass(frozen=True, slots=True, kw_only=True)
class FitParameters:
    """Settings for fitting head coefficients to a cooldown.

    Attributes:
        stage_1_crossings: Stage-1 temperatures whose first downward crossing times are fitted, in K.
        stage_2_crossings: Stage-2 temperatures whose first downward crossing times are fitted, in K.
        stage_1_max_temperature: Measured stage-1 readings at or above this temperature, in K, are excluded.
        sample_points: Number of evenly spaced times at which simulated crossings are located.
        record_duration_multiplier: Simulated duration as a multiple of the compared record's duration.
        crossing_duration_multiplier: Simulated duration as a multiple of the latest measured crossing time.
        solver: `scipy.integrate.solve_ivp` options. Absolute tolerances apply to logarithmic temperatures.
        final_solver: Tolerances and step limits for the final evaluation, applied where tighter than `solver`.
        optimizer: `scipy.optimize.least_squares` options.
        lower_bounds: Lowest coefficients searched.
        upper_bounds: Highest coefficients searched.
        duration: Simulated duration in s, replacing the larger of the two multiplied durations.
    """

    stage_1_crossings: tuple[float, ...]
    stage_2_crossings: tuple[float, ...]
    stage_1_max_temperature: float
    sample_points: int
    record_duration_multiplier: float
    crossing_duration_multiplier: float
    solver: Mapping[str, Any]
    final_solver: Mapping[str, Any]
    optimizer: Mapping[str, Any]
    lower_bounds: HeadCoefficients
    upper_bounds: HeadCoefficients
    duration: float | None = None


@dataclass(frozen=True, slots=True, kw_only=True)
class LoadMap:
    """Measured stage temperatures against stage heater powers.

    Attributes:
        temperature: Stage temperatures in K, shape (points, 2).
        heater_power: Stage heater powers in W, shape (points, 2).
    """

    temperature: NDArray[np.float64]
    heater_power: NDArray[np.float64]

    @classmethod
    def from_dataset(cls, dataset: xr.Dataset) -> Self:
        """Takes T1 and T2 in K and Q1 and Q2 in W along `point`, as returned by `read_load_map`."""
        return cls(
            temperature=np.column_stack((dataset.T1.values, dataset.T2.values)).astype(np.float64),
            heater_power=np.column_stack((dataset.Q1.values, dataset.Q2.values)).astype(np.float64),
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class HeadCalibration:
    """A head calibrated for one drive condition on a test stand.

    Attributes:
        coefficients: The fitted coefficients.
        design: The head design.
        load_map: The head's load map on the test stand.
        test_stand: The test stand's `fridge` table: its nodes and links.
        test_stand_stages: The test stand's stage-1 and stage-2 node names.
        warm_flange_temperature: Warm-flange temperature on the test stand, in K. The test stand's other nodes and
            the start of the calibration cooldown are at this temperature.
        provenance: Sources and calibration conditions.
    """

    coefficients: HeadCoefficients
    design: HeadDesignParameters
    load_map: LoadMap
    test_stand: Mapping[str, Any]
    test_stand_stages: tuple[str, str]
    warm_flange_temperature: float
    provenance: Mapping[str, str] = field(default_factory=dict)


def load_head_design(
    name: str, /, *, directory: str | os.PathLike[str] = DATA_DIRECTORY / "heads"
) -> HeadDesignParameters:
    """Loads a head design by name, such as "PT1.5", from a directory of design files.

    Raises:
        ConfigurationError: If no file or more than one file defines the name.
    """
    table = read_named_toml(name, directory=directory, kind="head designs")
    return HeadDesignParameters(
        name=table["name"],
        clamp_fraction=table["clamp_fraction"],
        regenerator_1=regenerator_layers(table["regenerator_1"]),
        regenerator_2=regenerator_layers(table["regenerator_2"]),
    )


def regenerator_layers(tables: list[Mapping[str, Any]]) -> tuple[RegeneratorLayer, ...]:
    """Builds regenerator layers from tables of `material` name, `mass`, and `length_fraction`."""
    return tuple(
        RegeneratorLayer(
            material=load_material(table["material"]), mass=table["mass"], length_fraction=table["length_fraction"]
        )
        for table in tables
    )


def load_shared_parameters(
    path: str | os.PathLike[str] = DATA_DIRECTORY / "shared_parameters.toml",
) -> SharedParameters:
    """Loads the shared parameters. Each profile table's `file` is a CSV path relative to `path`."""
    path = Path(path)
    table = tomllib.loads(path.read_text(encoding="utf-8"))
    unloaded, clamp, loaded = table["unloaded_profile"], table["clamp_response"], table["loaded_profiles"]
    unloaded_rows = read_csv(path.parent / unloaded["file"])
    clamp_rows = read_csv(path.parent / clamp["file"])
    loaded_rows = read_csv(path.parent / loaded["file"])
    return SharedParameters(
        energy_scale=table["energy_scale"],
        cold_blend_width=table["cold_blend_width"],
        cold_tail_span_multiplier=tuple(table["cold_tail_span_multiplier"]),
        energy_temperature_range=tuple(table["energy_temperature_range"]),
        energy_knots=table["energy_knots"],
        quadrature_order=table["quadrature_order"],
        unloaded_profile=UnloadedProfileParameters(
            position=unloaded_rows[:, 0],
            temperature=unloaded_rows[:, 1],
            reference_stage_temperature=tuple(unloaded["reference_stage_temperature"]),
        ),
        clamp_response=ClampResponseParameters(
            load=clamp_rows[:, 0],
            stage_1_temperature=clamp_rows[:, 1],
            stage_2_temperature=clamp_rows[:, 2],
            clamp_temperature=clamp_rows[:, 3],
        ),
        loaded_profiles=LoadedProfileParameters(
            load=np.asarray(loaded["load"], dtype=np.float64),
            position=loaded_rows[:, 0],
            temperature=loaded_rows[:, 1:].T,
            source_clamp_fraction=loaded["source_clamp_fraction"],
            reference_stage_temperature=tuple(loaded["reference_stage_temperature"]),
        ),
    )


def read_csv(path: Path) -> NDArray[np.float64]:
    """Reads a numeric CSV table below its one header line."""
    return np.loadtxt(path, delimiter=",", skiprows=1, ndmin=2)


def load_fit_parameters(path: str | os.PathLike[str] = DATA_DIRECTORY / "fit_parameters.toml") -> FitParameters:
    """Loads fit parameters. Each `bounds` entry is a [lower, upper] pair for one coefficient."""
    table = tomllib.loads(Path(path).read_text(encoding="utf-8"))
    bounds = table["bounds"]
    return FitParameters(
        stage_1_crossings=tuple(table["stage_1_crossings"]),
        stage_2_crossings=tuple(table["stage_2_crossings"]),
        stage_1_max_temperature=table["stage_1_max_temperature"],
        sample_points=table["sample_points"],
        record_duration_multiplier=table["record_duration_multiplier"],
        crossing_duration_multiplier=table["crossing_duration_multiplier"],
        solver=table["solver"],
        final_solver=table["final_solver"],
        optimizer=table["optimizer"],
        lower_bounds=HeadCoefficients(**{name: limits[0] for name, limits in bounds.items()}),
        upper_bounds=HeadCoefficients(**{name: limits[1] for name, limits in bounds.items()}),
        duration=table.get("duration"),
    )


def read_head_calibration(
    path: str | os.PathLike[str], /, *, designs: str | os.PathLike[str] = DATA_DIRECTORY / "heads"
) -> HeadCalibration:
    """Reads a calibration file, loading its design by name from `designs`."""
    table = tomllib.loads(Path(path).read_text(encoding="utf-8"))
    rows = np.asarray(table["load_map"]["table"], dtype=np.float64)
    first, second = table["test_stand_stages"]
    return HeadCalibration(
        coefficients=HeadCoefficients(**table["coefficients"]),
        design=load_head_design(table["design"], directory=designs),
        load_map=LoadMap(temperature=rows[:, :2], heater_power=rows[:, 2:]),
        test_stand=table["test_stand"],
        test_stand_stages=(first, second),
        warm_flange_temperature=table["warm_flange_temperature"],
        provenance=table.get("provenance", {}),
    )


def save_head_calibration(calibration: HeadCalibration, path: str | os.PathLike[str], /) -> None:
    """Saves a calibration as TOML, with one load-map point per line."""
    text = tomli_w.dumps(
        {
            "design": calibration.design.name,
            "test_stand_stages": list(calibration.test_stand_stages),
            "warm_flange_temperature": calibration.warm_flange_temperature,
            "coefficients": asdict(calibration.coefficients),
            "provenance": dict(calibration.provenance),
            "test_stand": dict(calibration.test_stand),
        }
    )
    load_map = np.column_stack((calibration.load_map.temperature, calibration.load_map.heater_power))
    rows = "".join(f"    {row.tolist()},\n" for row in load_map)
    columns = ["stage_1_temperature_K", "stage_2_temperature_K", "stage_1_heater_power_W", "stage_2_heater_power_W"]
    text += f"\n[load_map]\n{tomli_w.dumps({'columns': columns})}table = [\n{rows}]\n"
    Path(path).write_text(text, encoding="utf-8")
