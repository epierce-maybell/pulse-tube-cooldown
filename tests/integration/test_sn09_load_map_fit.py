"""The SN09 load-map fit: steady-state columns and recovery of a head's parameters from a load map it generated."""

import importlib.util
from pathlib import Path

import numpy as np
import pytest
import xarray as xr

_PATH = Path(__file__).parents[2] / "simulations" / "validation" / "sn09_load_map" / "fit.py"
_spec = importlib.util.spec_from_file_location("sn09_load_map_fit", _PATH)
fit = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(fit)

_WARM_TEMPERATURE = 294.2
# One F-70H, from data/compressors.md, with 80 g of helium over the head and the stand's 29.6 L of lines and buffers.
_DRIVE = fit.Drive(charge=0.08, warm_volume=29.6e-3, displacement_rate=(5.72e-3,), relief_swing=(1.93e6,))
_TRUE = fit.Parameters(
    packing_constant=(1.0, 1.2), conductance=(2.0e-10, 1.0e-11), exchanger_effectiveness=0.3, frequency=1.0
)
# Stage temperatures spanning the May load map.
_STAGE_TEMPERATURE = [(55.0, 4.0), (40.0, 4.0), (38.0, 2.3)]


@pytest.fixture(scope="module")
def head():
    return fit.read_head("PT1.1", directory=_PATH.parent / "baseline/head_model")


@pytest.fixture(scope="module")
def walls():
    return fit.build_walls()


def test_column_nodes_are_at_steady_state(head):
    solution = fit.solve_columns(head, _TRUE, _DRIVE, warm_temperature=_WARM_TEMPERATURE, stage_temperature=(40.0, 3.0))

    column_nodes = np.setdiff1d(np.arange(len(solution.temperature)), [0, head.stage_1, head.stage_2])
    np.testing.assert_allclose(solution.heat_flow[column_nodes], 0.0, atol=1e-6)
    assert solution.temperature[0] == _WARM_TEMPERATURE
    assert solution.temperature[[head.stage_1, head.stage_2]].tolist() == [40.0, 3.0]
    assert np.all(np.diff(solution.temperature) < 0)


def test_fit_recovers_parameters_of_its_own_load_map(head, walls):
    heat_flow = np.array(
        [
            fit.calc_stage_heat_flow(
                head, _TRUE, _DRIVE, walls, warm_temperature=_WARM_TEMPERATURE, stage_temperature=stages
            )[0]
            for stages in _STAGE_TEMPERATURE
        ]
    )
    stage_1, stage_2 = np.array(_STAGE_TEMPERATURE).T
    load_map = xr.Dataset(
        {
            "Q1": ("point", -heat_flow[:, 0]),
            "Q2": ("point", -heat_flow[:, 1]),
            "T1": ("point", stage_1),
            "TPR": ("point", np.full(len(stage_1), np.nan)),
            "T2": ("point", stage_2),
        }
    )
    initial = fit.Parameters(
        packing_constant=(0.98, 1.22), conductance=(2.01e-10, 1.01e-11), exchanger_effectiveness=0.33, frequency=1.0
    )

    result = fit.fit_load_map(head, _DRIVE, walls, load_map, warm_temperature=_WARM_TEMPERATURE, initial=initial)

    np.testing.assert_allclose(result.parameters.packing_constant, _TRUE.packing_constant, rtol=1e-3)
    np.testing.assert_allclose(result.parameters.conductance, _TRUE.conductance, rtol=1e-3)
    assert result.parameters.exchanger_effectiveness == pytest.approx(_TRUE.exchanger_effectiveness, rel=1e-3)
    assert result.parameters.frequency == _TRUE.frequency
    np.testing.assert_allclose(result.residual, 0.0, atol=1e-6)
