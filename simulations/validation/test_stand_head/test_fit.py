"""Requirements for the canonical PT1.5 fit to the test-stand load maps."""

import tomllib

import numpy as np
import pytest

from pulse_tube_cooldown.library.heads import load_head_model
from pulse_tube_cooldown.measurements import read_load_map
from simulations.validation.test_stand_head.fit import calc_stage_residual, fit_head


def _s0009_points():
    """S0009's 2026-05-14 load map with its run's coordinates: one compressor at 1.0 Hz."""
    points = read_load_map("data/test_stand/s0009/load_map_2026-05-14.dat")
    count = points.sizes["point"]
    return points.assign_coords(
        serial=("point", ["S0009"] * count),
        compressors=("point", [1] * count),
        valve_frequency=("point", [1.0] * count),
    )


def test_fit_recovers_the_parameters_behind_synthetic_heater_powers():
    """The fit returns the a, R, and c that generated heater powers at S0009's stage temperatures."""
    truth = {"packing": (0.05, 0.02), "resistance": (2e10, 4.5e10), "stage_coupling": 2.0}
    points = _s0009_points()
    unheated = points.assign(Q1=points.Q1 * 0, Q2=points.Q2 * 0)
    heat_into_stages = calc_stage_residual(unheated, **truth, workers=6)
    synthetic = points.assign(Q1=-heat_into_stages.Q1, Q2=-heat_into_stages.Q2)

    packing, resistance, stage_coupling = fit_head(
        synthetic, packing=(0.06, 0.022), resistance=(1.7e10, 5e10), workers=6
    )

    assert packing == pytest.approx(truth["packing"], rel=0.02)
    assert resistance == pytest.approx(truth["resistance"], rel=0.02)
    assert stage_coupling == pytest.approx(truth["stage_coupling"], rel=0.02)
    residual = calc_stage_residual(
        synthetic, packing=packing, resistance=resistance, stage_coupling=stage_coupling, workers=6
    )
    np.testing.assert_allclose(residual.Q1, 0.0, atol=0.05)
    np.testing.assert_allclose(residual.Q2, 0.0, atol=0.005)


def test_template_runs_the_canonical_pt1_5_head():
    """The test-stand template runs PT1.5 with the fitted parameter file and a nonzero stage coupling."""
    with open("data/test_stand/test_stand.toml", "rb") as file:
        (head,) = tomllib.load(file)["fridge"]["links"][0]["heads"]

    assert head["model"] == "PT1.5"
    assert head["parameters"] == "data/heads/pt1_5/parameters.toml"
    assert load_head_model("PT1.5").stage_coupling > 0
