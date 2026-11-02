"""Numerical requirements for the self-contained SN09 experiment."""

import importlib.util
import json
import tomllib
from pathlib import Path

import numpy as np
import pytest

_spec = importlib.util.spec_from_file_location("sn09_experiment", Path(__file__).with_name("fit.py"))
fit = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(fit)


def test_equilibrium_can_have_a_matrix_node_below_stage_2():
    """Preserves a real root excluded by the old plate-temperature bounds."""
    head = fit.read_head("PT1.1")
    inputs, load_map, warm, drive = fit.read_experiment()
    stages = (float(load_map.T1[0]), float(load_map.T2[0]))
    parameters = fit.Parameters(
        packing_constant=(0.1555441973538624, 0.01),
        conductance=(2.5092762169316817e-10, 1e-11),
        frequency=inputs["frequency"],
    )
    guess = np.array(
        [warm, 223.97518913, 92.18649583, stages[0], 50.2584913, 34.57765293, 6.85818033, 2.88528551, stages[1]]
    )
    column = fit.solve_columns(head, parameters, drive, warm_temperature=warm, stage_temperature=stages, guess=guess)
    np.testing.assert_allclose(column.heat_flow[[1, 2, 4, 5, 6, 7]], 0, atol=1e-7)
    assert column.temperature[7] < stages[1]
    assert fit.calc_growth_rate(head, parameters, drive, column) > 0.1


@pytest.mark.parametrize("cooling", [False, True])
def test_synthetic_load_map_recovers_its_parameters(cooling):
    """Recovers four known coefficients from a model-generated load map."""
    import xarray as xr

    head = fit.read_head("PT1.1")
    _, _, warm, drive = fit.read_experiment()
    known = fit.Parameters(packing_constant=(1.0, 1.2), conductance=(2e-10, 1e-11), frequency=1.0)
    stages = [(55.0, 4.0), (53.0, 2.4), (40.0, 4.0), (38.0, 2.3)]
    walls = fit.build_walls()
    guesses = [None] * len(stages)
    bounds = None
    if cooling:
        known = fit.Parameters(packing_constant=(0.13, 0.1), conductance=(1e-10, 4.4e-11), frequency=1.0)
        seeds = json.loads(Path(__file__).with_name("cooling_seed.json").read_text())
        guesses = [np.array(seeds["temperatures"][i]) for i in [0, 2, 3, 5]]
        bounds = (4.2e-10, 4.6e-10)
    flows = np.array(
        [
            fit.calc_stage_heat_flow(
                head, known, drive, walls, warm_temperature=warm, stage_temperature=stage, guess=guess
            )[0]
            for stage, guess in zip(stages, guesses, strict=True)
        ]
    )
    data = xr.Dataset(
        {
            "T1": ("point", np.array(stages)[:, 0]),
            "T2": ("point", np.array(stages)[:, 1]),
            "Q1": ("point", -flows[:, 0]),
            "Q2": ("point", -flows[:, 1]),
        }
    )
    initial = fit.Parameters(packing_constant=(0.7, 1.6), conductance=(1.5e-10, 1.4e-11), frequency=1.0)
    if cooling:
        initial = fit.Parameters(packing_constant=(0.1, 0.12), conductance=(0.8e-10, 5.4e-11), frequency=1.0)
    result = fit.fit_load_map(
        head, drive, walls, data, warm_temperature=warm, initial=initial, ratio_bounds=bounds, guesses=guesses
    )
    np.testing.assert_allclose(result.parameters.packing_constant, known.packing_constant, rtol=1e-5)
    np.testing.assert_allclose(result.parameters.conductance, known.conductance, rtol=1e-5)
    np.testing.assert_allclose(result.residual, 0, atol=1e-6)


def test_candidate_balances_with_an_independently_solved_drive():
    """Checks the archived candidate against the library equations and its reported residuals."""
    parameters = tomllib.loads(Path(__file__).with_name("candidate_parameters.toml").read_text())
    report = json.loads(Path(__file__).with_name("results.json").read_text())
    head = fit.read_head("PT1.1")
    _, load_map, warm, drive = fit.read_experiment()
    candidate = fit.Parameters(
        packing_constant=tuple(parameters["packing_constants"]),
        conductance=tuple(parameters["conductances"]),
        frequency=parameters["frequency"],
    )
    for point, row in enumerate(report["points"]):
        t = np.array(row["temperature_K"])
        mean, swing = fit.calc_pressures(head, candidate, drive, t)
        column = fit.evaluate_head(head, candidate, drive, t, mean, swing)
        np.testing.assert_allclose(column.heat_flow[fit.column_indices(head)], 0, atol=1e-7)
        wall = fit.calc_link_heat_flow(fit.build_walls(), {fit.WARM_FLANGE: warm, fit.STAGE_1: t[3], fit.STAGE_2: t[8]})
        residual = column.heat_flow[[3, 8]] + [
            wall[fit.STAGE_1] + float(load_map.Q1[point]),
            wall[fit.STAGE_2] + float(load_map.Q2[point]),
        ]
        np.testing.assert_allclose(residual, row["residual_W"], atol=1e-7)
        np.testing.assert_allclose(column.utilization, row["utilization"], rtol=1e-9)
        for recomputed, archived in zip(column.capacity_ratio, row["capacity_ratio"], strict=True):
            np.testing.assert_allclose(recomputed, archived, rtol=1e-9)
        assert fit.calc_growth_rate(head, candidate, drive, column) < 0
        changed = fit.Parameters(
            packing_constant=(candidate.packing_constant[0] * 2, candidate.packing_constant[1]),
            conductance=(candidate.conductance[0] * 2, candidate.conductance[1]),
            frequency=1.0,
        )
        repeat = fit.solve_columns(head, changed, drive, warm_temperature=warm, stage_temperature=(t[3], t[8]), guess=t)
        np.testing.assert_allclose(repeat.temperature, t, atol=1e-7)
        assert repeat.heat_flow[-1] == pytest.approx(column.heat_flow[-1], abs=1e-7)


def test_archived_capacity_ratio_matches_coupling_findings():
    """The May-fit profiles' pore-gas to matrix heat-capacity ratio, estimated in coupling_findings.txt."""
    report = json.loads(Path(__file__).with_name("results.json").read_text())
    middle_load = [row for row in report["points"] if row["heater_power_W"][1] == 0.5]
    assert len(middle_load) == 2
    for row in middle_load:
        stage_2_column = np.array(row["capacity_ratio"][1])
        assert 0.04 <= stage_2_column.min() <= 0.08
        assert 0.38 <= stage_2_column.max() <= 0.59
