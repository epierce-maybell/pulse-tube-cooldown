"""Behavioral contract for the provisional SN09 bare-head experiment."""

import importlib.util
from pathlib import Path

import numpy as np
import pytest
import xarray as xr

from pulse_tube_cooldown.engine.solver import simulate_cooldown

_PATH = Path(__file__).parents[2] / "simulations/validation/sn09_bare_head/run.py"
_spec = importlib.util.spec_from_file_location("sn09_bare_head", _PATH)
stand = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(stand)


def test_network_has_inventory_matrix_only_columns_and_four_walls():
    model = stand.build_model("may", max_time=1.0)
    graph = model.graphs[0]
    cooler, *walls = graph.links
    unit = cooler.heads[0]
    composition = {node.name: {material.name: mass for material, mass in node.composition} for node in model.nodes}
    assert composition["Stage 1"] == {
        "OFHC Cu RRR100": pytest.approx(0.3928 + 1.191),
        "Al 6061-T6": 0.4536,
        "304 SS": 0.0709,
    }
    assert composition["Stage 2"] == {"304 SS": 0.4218, "OFHC Cu RRR100": 0.0397}
    assert len(model.nodes) == 3 + len(unit.matrix_nodes)
    for node in unit.matrix_nodes:
        assert node in model.nodes
        assert len(node.composition) == 1
    assert unit.clamp_node in composition
    assert model.fixed_boundary_nodes == (graph.get_node_index("Warm flange"),)
    assert len(walls) == 4
    # CAD: outer radius 32.11 mm, length 180 mm; estimated wall 0.5 mm.
    wall = next(wall for wall in walls if wall.name == "stage1_regenerator wall")
    assert wall.area == pytest.approx(np.pi * (0.03211**2 - 0.03161**2))
    assert wall.length == pytest.approx(0.180)


def test_august_changes_compressor_count_with_same_head_parameters():
    may = stand.build_model("may", max_time=1.0).graphs[0].links[0]
    august = stand.build_model("august", max_time=1.0).graphs[0].links[0]
    assert len(may.compressors) == 1
    assert len(august.compressors) == 2
    assert august.compressors == may.compressors * 2
    assert august.heads[0].packing_constants == may.heads[0].packing_constants
    assert august.heads[0].conductances == may.heads[0].conductances
    assert august.heads[0].frequency == may.heads[0].frequency == 1.0


def test_measurements_exclude_appended_may_run_and_use_timestamp_clock():
    data = stand.read_measurements("may")
    assert data.time[0] == 0.0
    assert str(data.attrs["origin"]) == "2026-05-14T10:16:45"
    assert data.time[-1] < (4 * 3600 + 2 * 60 + 32)
    assert data.attrs["end_exclusive"] == "2026-05-14T14:19:17"
    np.testing.assert_array_equal(np.diff(data.time), 2.0)


def test_comparison_uses_valid_t1_below_219_and_all_finite_t2_without_tpr():
    measured = xr.Dataset(
        {
            "T1": ("time", [220.0, 218.0, np.nan, 100.0]),
            "T2": ("time", [300.0, np.nan, 20.0, 4.0]),
            "TPR": ("time", [999.0] * 4),
        },
        coords={"time": [0.0, 1.0, 2.0, 3.0]},
    )
    temperature = xr.DataArray(
        [[110.0, 150.0], [50.0, 2.0]],
        dims=("time", "node"),
        coords={"time": [0.0, 3.0], "node": ["Stage 1", "Stage 2"]},
    )
    comparison = stand.compare_cooldown(measured, temperature)
    assert set(comparison) == {"T1", "T2"}
    assert comparison["T1"]["samples"] == 2
    assert comparison["T2"]["samples"] == 3
    assert comparison["T1"]["sum_squared_log_residual"] == pytest.approx(
        np.log(90.0 / 218.0) ** 2 + np.log(50.0 / 100.0) ** 2
    )
    # Predictions are interpolated to the measured timestamps.
    assert comparison["T2"]["sum_squared_log_residual"] == pytest.approx(
        np.log(150.0 / 300.0) ** 2 + np.log((150.0 - 148.0 * 2 / 3) / 20.0) ** 2 + np.log(2.0 / 4.0) ** 2
    )


def test_network_runs_with_fixed_flange_and_conservative_links():
    model = stand.build_model("may", max_time=10.0)
    result = simulate_cooldown(model)
    assert result.status == "max_time"
    temperature = result.data.temperature
    np.testing.assert_allclose(temperature.sel(node="Warm flange"), model.initial_temperature[0])
    assert float(temperature.sel(node="Stage 1")[-1]) < model.initial_temperature[0]
    assert float(temperature.sel(node="Stage 2")[-1]) < model.initial_temperature[0]
    np.testing.assert_allclose(result.data.heat_flow.sum("node"), 0.0, atol=1e-9)
