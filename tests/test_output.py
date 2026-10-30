"""Tests for saving a Result."""

import dataclasses

import numpy as np
import xarray as xr
from numpy.typing import NDArray

from pulse_tube_cooldown.core.types import Material, Model, Node, Phase, PhaseSolution, Procedure
from pulse_tube_cooldown.engine.solver import build_result
from pulse_tube_cooldown.output import save_output

_CONSTANT = Material(
    name="Constant",
    density=1.0,
    specific_heat=lambda t: 1.0 + 0.0 * t,
    thermal_conductivity=lambda t: 0.0 * t,
    conductivity_integral=lambda t: 0.0 * t,
)


@dataclasses.dataclass(frozen=True, slots=True)
class _Conductance:
    name: str
    nodes: tuple[str, ...]
    conductance: float

    def calc_heat_flow(self, temperature: NDArray[np.float64]) -> NDArray[np.float64]:
        flow = self.conductance * (temperature[1] - temperature[0])
        return np.array([flow, -flow])


_META = {"experiment": "Toy switch", "version": "abc1234"}


def _switching_result():
    # Precool runs a 0.5 W/K strap and a 0.1 W/K pump; cooldown drops the pump and swaps in a 1 W/K strap.
    strap = _Conductance("strap", ("stage", "bath"), 0.5)
    pump = _Conductance("pump", ("stage", "bath"), 0.1)
    precool = Phase("precool", (strap, pump), lambda time, temperature: temperature[0] - 50.0)
    cooldown = Phase("cooldown", (_Conductance("strap", ("stage", "bath"), 1.0),), None)
    nodes = (Node("stage", ((_CONSTANT, 2.0),)), Node("bath", ((_CONSTANT, 1.0),)))
    model = Model(
        _META,
        nodes,
        Procedure(phases=(precool, cooldown), stop_condition=None),
        np.array([300.0, 4.0]),
        (1,),
        {"max_time": 5.0, "save_every": 2.0},
    )
    solutions = [
        PhaseSolution(
            phase=precool,
            end="until",
            time=np.array([0.0, 2.0, 3.0]),
            temperature=np.array([[300.0, 4.0], [100.0, 4.0], [50.0, 4.0]]),
            message="precool ended",
        ),
        PhaseSolution(
            phase=cooldown,
            end="max_time",
            time=np.array([3.0, 5.0]),
            temperature=np.array([[50.0, 4.0], [10.0, 4.0]]),
            message="cooldown ended",
        ),
    ]
    return build_result(model, solutions)


def test_saved_output_reopens_identical_to_the_result(tmp_path):
    result = _switching_result()
    path = tmp_path / "run.nc"

    save_output(result, path)

    with xr.open_dataset(path) as saved:
        xr.testing.assert_identical(saved.load(), result.data)


def test_saved_output_keeps_both_sides_of_a_phase_switch(tmp_path):
    path = tmp_path / "run.nc"

    save_output(_switching_result(), path)

    with xr.open_dataset(path) as saved:
        np.testing.assert_array_equal(saved.time, [0.0, 2.0, 3.0, 3.0, 5.0])
        assert list(saved.phase.values) == ["precool", "precool", "precool", "cooldown", "cooldown"]
        # At 50 K: 0.5 W/K * (4 K - 50 K) before the switch and 1 W/K * (4 K - 50 K) after.
        np.testing.assert_allclose(saved.heat_flow.isel(time=[2, 3]).sel(link="strap", node="stage"), [-23.0, -46.0])
        assert np.isnan(saved.heat_flow.isel(time=3).sel(link="pump")).all()


def test_saved_output_records_status_and_meta(tmp_path):
    path = tmp_path / "run.nc"

    save_output(_switching_result(), path)

    with xr.open_dataset(path) as saved:
        assert saved.attrs["status"] == "max_time"
        assert saved.attrs["message"] == "cooldown ended"
        assert saved.attrs["experiment"] == "Toy switch"
        assert saved.attrs["version"] == "abc1234"
        assert saved.connected.dtype == bool
