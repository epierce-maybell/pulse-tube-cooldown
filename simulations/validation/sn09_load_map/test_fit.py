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
    """A stable fitted root has a matrix minimum below the measured plate."""
    head = fit.read_head("PT1.1")
    inputs, load_map, warm, drive = fit.read_experiment()
    stages = (float(load_map.T1[0]), float(load_map.T2[0]))
    report = json.loads(Path(__file__).with_name("results.json").read_text())
    parameters = fit.Parameters(
        packing_constant=tuple(report["packing_constants"]),
        conductance=tuple(report["conductances_W_per_Pa2"]),
        exchanger_effectiveness=report["exchanger_effectiveness"],
        frequency=inputs["frequency"],
    )
    guess = np.array(report["points"][0]["temperature_K"])
    column = fit.solve_columns(head, parameters, drive, warm_temperature=warm, stage_temperature=stages, guess=guess)
    np.testing.assert_allclose(column.heat_flow[fit.column_indices(head)], 0, atol=1e-7)
    assert column.temperature[-2] < stages[1]
    assert fit.calc_growth_rate(head, parameters, drive, column) < 0


def test_candidate_balances_with_an_independently_solved_drive():
    """Checks the archived candidate against the library equations and its reported residuals."""
    parameters = tomllib.loads(Path(__file__).with_name("candidate_parameters.toml").read_text())
    report = json.loads(Path(__file__).with_name("results.json").read_text())
    head = fit.read_head("PT1.1")
    _, load_map, warm, drive = fit.read_experiment()
    candidate = fit.Parameters(
        packing_constant=tuple(parameters["packing_constants"]),
        conductance=tuple(parameters["conductances"]),
        exchanger_effectiveness=parameters["exchanger_effectiveness"],
        frequency=parameters["frequency"],
    )
    for point, row in enumerate(report["points"]):
        t = np.array(row["temperature_K"])
        mean, swing = fit.calc_pressures(head, candidate, drive, t)
        column = fit.evaluate_head(head, candidate, drive, t, mean, swing)
        np.testing.assert_allclose(column.heat_flow[fit.column_indices(head)], 0, atol=1e-7)
        wall = fit.calc_link_heat_flow(
            fit.build_walls(), {fit.WARM_FLANGE: warm, fit.STAGE_1: t[head.stage_1], fit.STAGE_2: t[head.stage_2]}
        )
        residual = column.heat_flow[[head.stage_1, head.stage_2]] + [
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
            exchanger_effectiveness=candidate.exchanger_effectiveness,
            frequency=1.0,
        )
        repeat = fit.solve_columns(
            head, changed, drive, warm_temperature=warm, stage_temperature=(t[head.stage_1], t[head.stage_2]), guess=t
        )
        np.testing.assert_allclose(repeat.temperature, t, atol=1e-7)
        assert repeat.heat_flow[-1] == pytest.approx(column.heat_flow[-1], abs=1e-7)


def test_archived_capacity_ratio_matches_coupling_findings():
    """The May-fit profiles' pore-gas to matrix heat-capacity ratio, estimated in coupling_findings.txt."""
    report = json.loads((Path(__file__).parent / "baseline/results.json").read_text())
    middle_load = [row for row in report["points"] if row["heater_power_W"][1] == 0.5]
    assert len(middle_load) == 2
    for row in middle_load:
        stage_2_column = np.array(row["capacity_ratio"][1])
        assert 0.04 <= stage_2_column.min() <= 0.08
        assert 0.38 <= stage_2_column.max() <= 0.59


def test_refit_improves_stage_2_and_retains_comparison():
    """The accepted May refit beats the preserved stage-2 RMS and includes the uncoupled control."""
    folder = Path(__file__).parent
    baseline = json.loads((folder / "baseline/results.json").read_text())
    result = json.loads((folder / "results.json").read_text())
    control = json.loads((folder / "uncoupled/results.json").read_text())
    assert result["rms_residual_W"][1] < baseline["rms_residual_W"][1]
    assert control["cold_exchanger_coupling"] is False
    assert result["cold_exchanger_coupling"] is True
    assert result["exchanger_effectiveness"] > 0


def test_decomposition_reconstructs_stage_power():
    """Each archived component balance reconstructs the available power independently of the fitter."""
    folder = Path(__file__).parent
    for variant in (folder, folder / "uncoupled"):
        result = json.loads((variant / "results.json").read_text())
        decomposition = json.loads((variant / "decomposition.json").read_text())
        for point, terms in zip(result["points"], decomposition, strict=True):
            np.testing.assert_allclose(
                [terms["Q1_W"], terms["Q2_W"]],
                np.array(point["heater_power_W"]) - point["residual_W"],
                atol=1e-7,
                rtol=0,
            )


def test_uncoupled_control_matches_hand_computed_finite_capacity_head():
    """The task-6.2 hand calculation fixes the uncoupled control's plate boundary."""
    for effectiveness in (0.0, 0.8):
        heat_flow = fit.calc_uncoupled_heat_flow(
            np.array([300.0, 60.0, 40.0, 20.0, 8.0, 4.0]),
            (np.array([800.0]), np.array([40.0, 15.0])),
            (np.array([7.163822e-3]), np.array([0.1480389, 0.9654703])),
            np.array([5.503374e-3, 4.416387e-3, 4.106599e-3, 3.106279e-3, 1.748532e-3]),
            (30.0, 1.5),
            (0.8, 1.3),
            effectiveness,
            1.4e6,
            1.2,
        )
        np.testing.assert_allclose(heat_flow, [-295.4320, 308.6979, -565.9027, 350.9011, 171.4323, 30.30343], rtol=1e-6)
        assert heat_flow.sum() == pytest.approx(0.0, abs=1e-9)
