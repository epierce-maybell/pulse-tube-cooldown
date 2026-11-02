"""Builds a pulse-tube unit from its model and parameter files."""

import os
import tomllib
from pathlib import Path

from pulse_tube_cooldown.core.errors import ConfigurationError
from pulse_tube_cooldown.library.links.cooler.columns import load_column
from pulse_tube_cooldown.library.links.cooler.head import PulseTubeHead
from pulse_tube_cooldown.library.toml_files import read_named_toml

__all__ = ["load_pulse_tube_head"]


def load_pulse_tube_head(
    *,
    name: str,
    stages: tuple[str, str],
    model: str,
    parameters: str | os.PathLike[str],
    directory: str | os.PathLike[str] = Path(__file__).parents[2] / "heads/data",
) -> PulseTubeHead:
    """Loads one unit, expanding each physical layer into its configured matrix nodes.

    Args:
        name: Unit's name, used in generated matrix node names.
        stages: Names of the stage-1 and stage-2 nodes.
        model: Model's name in the library, such as "PT1.1".
        parameters: Per-unit TOML file with packing_constants, conductances, exchanger_effectiveness, and frequency.
        directory: Directory containing the head model TOML files.

    Returns:
        The unit, including its generated matrix nodes.

    Raises:
        ConfigurationError: A model is unknown, its column indices, fractions, or node counts are invalid,
            or its clamp position is outside the stage-2 column.
    """
    data = read_named_toml(model, directory=directory, kind="heads")
    unit = tomllib.loads(Path(parameters).read_text(encoding="utf-8"))
    columns = tuple(
        load_column(data[column], head_name=name, stage=stage)
        for stage, column in enumerate(("stage_1_column", "stage_2_column"), start=1)
    )
    clamp_fraction = float(data["clamp_fraction"])
    if not 0 <= clamp_fraction <= 1:
        raise ConfigurationError(f"{model}: clamp_fraction must lie between 0 and 1.")
    first_volume, second_volume = data["pulse_tube_volumes"]
    first_packing, second_packing = unit["packing_constants"]
    first_conductance, second_conductance = unit["conductances"]
    return PulseTubeHead(
        name=name,
        stages=stages,
        columns=(columns[0], columns[1]),
        pulse_tube_volumes=(float(first_volume), float(second_volume)),
        clamp_fraction=clamp_fraction,
        packing_constants=(float(first_packing), float(second_packing)),
        conductances=(float(first_conductance), float(second_conductance)),
        exchanger_effectiveness=float(unit["exchanger_effectiveness"]),
        frequency=float(unit["frequency"]),
    )
