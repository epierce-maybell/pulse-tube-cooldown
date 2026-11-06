"""Requirements for a manifold's drive: the charge, the mean pressure that holds it, and the supply-equation swing."""

import math
from itertools import accumulate

import numpy as np
import pytest

from pulse_tube_cooldown.library.gases import helium
from pulse_tube_cooldown.library.heads import load_head_model
from pulse_tube_cooldown.library.links.cooler.drive import (
    Compressor,
    calc_charge,
    calc_mean_pressure,
    calc_pressure_swing,
)
from pulse_tube_cooldown.library.links.cooler.head import Head, calc_gas_mass, calc_mass_swing

# Inputs from data/compressors.md: one F-70H's swept displacement rate, the middle of FP1's 18.6 to 20.0 bar relief
# swing, the test stand's 19.6 L of flexline, and a fill of 17.0 bar absolute at rest.
_DISPLACEMENT_RATE = 5.72e-3
_RELIEF_SWING = 1.93e6
_WARM_VOLUME = 19.6e-3
_FILL_PRESSURE = 1.70e6
_ROOM_TEMPERATURE = 295.0
_RELIEF_WIDTH = 5e4

# Placeholder valve resistances in Pa s/m^3 and a 1.0 Hz valve frequency, pending the fit of task 2.3.
_RESISTANCE = (1e10, 4e10)
_FREQUENCY = 1.0

_MODEL = load_head_model("PT1.5")


def _chain_temperature(model, stage_temperature, warm_temperature=_ROOM_TEMPERATURE):
    """Chain temperatures with each section's nodes linear in position between its ends' temperatures."""
    ends = (warm_temperature, *stage_temperature)
    sections = []
    for section, warm, cold in zip(model.regenerators, ends[:-1], ends[1:], strict=True):
        fractions = [
            layer.length_fraction / layer.node_count for layer in section.layers for _ in range(layer.node_count)
        ]
        midpoints = np.array(
            [end - fraction / 2 for end, fraction in zip(accumulate(fractions), fractions, strict=True)]
        )
        sections.append(warm + (cold - warm) * midpoints)
    return np.concatenate(([warm_temperature], sections[0], [ends[1]], sections[1], [ends[2]]))


def _pressures(models, stage_temperature, compressors):
    """Mean pressure and swing of a manifold whose heads share stage temperatures."""
    temperature = _chain_temperature(_MODEL, stage_temperature)
    charge = calc_charge(
        models, warm_volume=_WARM_VOLUME, fill_pressure=_FILL_PRESSURE, fill_temperature=_ROOM_TEMPERATURE
    )
    mean_pressure = calc_mean_pressure(
        lambda pressure: sum(calc_gas_mass(model, temperature, pressure) for model in models),
        charge,
        warm_volume=_WARM_VOLUME,
        warm_temperature=_ROOM_TEMPERATURE,
    )

    heads = [Head(model=model, packing=(0.0, 0.0), resistance=_RESISTANCE, frequency=_FREQUENCY) for model in models]

    def head_mass_flow(swing):
        return sum(head.frequency * calc_mass_swing(head, temperature, mean_pressure, swing).sum() for head in heads)

    swing = calc_pressure_swing(
        head_mass_flow, mean_pressure, _ROOM_TEMPERATURE, compressors=compressors, relief_width=_RELIEF_WIDTH
    )
    return mean_pressure, swing


def test_charge_fills_every_volume_at_rest():
    # HEPAK's density at 295 K and 1.70 MPa times the warm volume and both heads' pores, tubes, reservoirs, dead
    # volume, and heat-exchanger voids, from pt1_5.toml.
    pores = 2.8459e-4 + 6.0081e-5 + 1.1956e-5 + 2 * 2.3911e-5 + 1.4945e-5
    head = pores + 2.5077e-4 + 2.0086e-4 + 2 * 5.2800e-3 + 4.5e-5 + 39.3e-6 + 10.8e-6

    charge = calc_charge(
        [_MODEL, _MODEL], warm_volume=_WARM_VOLUME, fill_pressure=_FILL_PRESSURE, fill_temperature=_ROOM_TEMPERATURE
    )

    assert charge == pytest.approx(2.7524272578 * (_WARM_VOLUME + 2 * head), rel=1e-4)


def test_mean_pressure_at_rest_is_fill_pressure():
    mean_pressure, _ = _pressures(
        [_MODEL],
        (_ROOM_TEMPERATURE, _ROOM_TEMPERATURE),
        [Compressor(displacement_rate=_DISPLACEMENT_RATE, relief_swing=_RELIEF_SWING)],
    )

    assert mean_pressure == pytest.approx(_FILL_PRESSURE, rel=1e-8)


@pytest.mark.parametrize("stage_2_temperature", [4.0, 2.0])
def test_mean_pressure_holds_the_charge_when_cold(stage_2_temperature):
    # Helium at 2 K melts above about 3 MPa, so the search stays within fluid pressures.
    temperature = _chain_temperature(_MODEL, (40.0, stage_2_temperature))
    charge = calc_charge(
        [_MODEL], warm_volume=_WARM_VOLUME, fill_pressure=_FILL_PRESSURE, fill_temperature=_ROOM_TEMPERATURE
    )

    mean_pressure = calc_mean_pressure(
        lambda pressure: calc_gas_mass(_MODEL, temperature, pressure),
        charge,
        warm_volume=_WARM_VOLUME,
        warm_temperature=_ROOM_TEMPERATURE,
    )

    held = calc_gas_mass(_MODEL, temperature, mean_pressure) + _WARM_VOLUME * float(
        helium.density(_ROOM_TEMPERATURE, mean_pressure)
    )
    assert held == pytest.approx(charge, rel=1e-9)


def _linear_flow(swing_at_balance, mean_pressure, displacement_rate):
    """Head flow proportional to the swing, meeting the compressors' flow at `swing_at_balance`."""
    supply = float(helium.density(_ROOM_TEMPERATURE, mean_pressure - swing_at_balance / 2)) * displacement_rate
    return lambda swing: supply * swing / swing_at_balance


def test_swing_below_relief_meets_compressor_flow():
    mean_pressure = 1.5e6
    compressors = [Compressor(displacement_rate=_DISPLACEMENT_RATE, relief_swing=_RELIEF_SWING)]

    swing = calc_pressure_swing(
        _linear_flow(0.8e6, mean_pressure, _DISPLACEMENT_RATE),
        mean_pressure,
        _ROOM_TEMPERATURE,
        compressors=compressors,
        relief_width=_RELIEF_WIDTH,
    )

    assert swing == pytest.approx(0.8e6, rel=1e-8)


def test_swing_far_above_relief_holds_at_relief():
    compressors = [Compressor(displacement_rate=_DISPLACEMENT_RATE, relief_swing=_RELIEF_SWING)]

    swing = calc_pressure_swing(
        _linear_flow(3.0e6, 1.7e6, _DISPLACEMENT_RATE),
        1.7e6,
        _ROOM_TEMPERATURE,
        compressors=compressors,
        relief_width=_RELIEF_WIDTH,
    )

    assert swing == pytest.approx(_RELIEF_SWING, abs=1.0)


def test_compressors_add_flow_under_the_lowest_relief():
    # Unrelieved, the heads would balance both compressors' flow exactly at the lower relief swing, where the smooth
    # relief gives the relief swing less w ln 2.
    compressors = [
        Compressor(displacement_rate=_DISPLACEMENT_RATE, relief_swing=2.0e6),
        Compressor(displacement_rate=_DISPLACEMENT_RATE, relief_swing=1.9e6),
    ]

    swing = calc_pressure_swing(
        _linear_flow(1.9e6, 1.7e6, 2 * _DISPLACEMENT_RATE),
        1.7e6,
        _ROOM_TEMPERATURE,
        compressors=compressors,
        relief_width=_RELIEF_WIDTH,
    )

    assert swing == pytest.approx(1.9e6 - _RELIEF_WIDTH * math.log(2), rel=1e-8)


def _cooldown(points):
    """Mean pressures and swings of one head on one compressor, its stages cooling geometrically to 40 K and 4 K."""
    compressor = Compressor(displacement_rate=_DISPLACEMENT_RATE, relief_swing=_RELIEF_SWING)
    progress = np.linspace(0.0, 1.0, points)
    stages = np.column_stack(
        (
            _ROOM_TEMPERATURE * (40.0 / _ROOM_TEMPERATURE) ** progress,
            _ROOM_TEMPERATURE * (4.0 / _ROOM_TEMPERATURE) ** progress,
        )
    )
    mean_pressure, swing = np.array([_pressures([_MODEL], stage, [compressor]) for stage in stages]).T
    return mean_pressure, swing


@pytest.fixture(scope="module")
def cooldown():
    return _cooldown(41)


def test_warm_head_runs_at_relief(cooldown):
    _, swing = cooldown

    assert swing[0] == pytest.approx(_RELIEF_SWING, abs=_RELIEF_WIDTH)


def test_cold_head_lowers_mean_pressure_and_swing(cooldown):
    mean_pressure, swing = cooldown

    assert mean_pressure[-1] < _FILL_PRESSURE
    assert swing[-1] < _RELIEF_SWING - _RELIEF_WIDTH


def test_cooldown_pressures_are_smooth(cooldown):
    # Along a path with continuous derivatives, second differences shrink against first differences in proportion to
    # the step; at a kink they do not.
    def roughness(values):
        return np.max(np.abs(np.diff(values, 2))) / np.max(np.abs(np.diff(values)))

    fine = _cooldown(161)

    assert np.all(np.diff(fine[0]) < 0)
    for coarse_values, fine_values in zip(cooldown, fine, strict=True):
        assert roughness(fine_values) <= 0.5 * roughness(coarse_values)
