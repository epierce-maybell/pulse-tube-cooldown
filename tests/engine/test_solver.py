"""Tests for integrating a model's phases."""

import dataclasses

import numpy as np
import pytest
from numpy.typing import NDArray

from pulse_tube_cooldown.core.errors import ConfigurationError
from pulse_tube_cooldown.core.types import Material, Model, Node, Phase, PhaseSolution, Procedure
from pulse_tube_cooldown.engine.solver import integrate_phase

_CONSTANT = Material(
    name="Constant",
    specific_heat=lambda t: 1.0 + 0.0 * t,
    thermal_conductivity=lambda t: 0.0 * t,
    conductivity_integral=lambda t: 0.0 * t,
)


@dataclasses.dataclass(frozen=True, slots=True)
class _Conductance:
    """Constant conductance G from the first node to the second."""

    name: str
    nodes: tuple[str, str]
    conductance: float

    def calc_heat_flow(self, temperature: NDArray[np.float64]) -> NDArray[np.float64]:
        flow = self.conductance * (temperature[1] - temperature[0])
        return np.array([flow, -flow])


@dataclasses.dataclass(frozen=True, slots=True)
class _Runaway:
    """Heat input k * T^2 into a single node, which diverges in finite time."""

    name: str
    nodes: tuple[str]
    coefficient: float

    def calc_heat_flow(self, temperature: NDArray[np.float64]) -> NDArray[np.float64]:
        return self.coefficient * temperature**2


# A 2 J/K stage relaxes toward a 4 K bath through 0.5 W/K, so tau = C / G = 4 s and
# T(t) = 4 + (T0 - 4) * exp(-(t - t0) / 4).
_STAGE = Node("stage", ((_CONSTANT, 2.0),))
_BATH = Node("bath", ((_CONSTANT, 1.0),))
_LINK = _Conductance("strap", ("stage", "bath"), 0.5)
_TAU = 4.0
_OPTIONS = {"solver": "Radau", "max_time": 20.0, "save_every": 2.0, "rtol": 1e-10, "atol": 1e-12}


def _relaxation(time: NDArray[np.float64], start_time: float, start_temperature: float) -> NDArray[np.float64]:
    return 4.0 + (start_temperature - 4.0) * np.exp(-(time - start_time) / _TAU)


def _stage_below(threshold: float):
    def condition(time: float, temperature: NDArray[np.float64]) -> float:
        return float(temperature[0] - threshold)

    return condition


def _after(cutoff: float):
    def condition(time: float, temperature: NDArray[np.float64]) -> float:
        return cutoff - time

    return condition


def _model(phases: tuple[Phase, ...], stop_condition=None, options=_OPTIONS) -> Model:
    return Model(
        nodes=(_STAGE, _BATH),
        procedure=Procedure(phases=phases, stop_condition=stop_condition),
        initial_temperature=np.array([300.0, 4.0]),
        fixed_boundary_nodes=(1,),
        solver_options=options,
    )


def _final_phase() -> Phase:
    return Phase("cooldown", (_LINK,), None)


def test_relaxation_matches_analytic_solution_until_max_time():
    phase = _final_phase()
    model = _model((phase,))

    solution = integrate_phase(model, phase, 0.0, model.initial_temperature)

    assert isinstance(solution, PhaseSolution)
    assert solution.end == "max_time"
    assert solution.phase is phase
    np.testing.assert_allclose(solution.time, np.arange(0.0, 20.1, 2.0))
    assert solution.temperature.shape == (11, 2)
    np.testing.assert_allclose(solution.temperature[:, 0], _relaxation(solution.time, 0.0, 300.0), rtol=1e-6)
    np.testing.assert_array_equal(solution.temperature[:, 1], 4.0)


def test_phase_starting_later_samples_from_its_start_and_ends_at_run_max_time():
    phase = _final_phase()
    model = _model((phase,))

    solution = integrate_phase(model, phase, 5.0, np.array([100.0, 4.0]))

    assert solution.end == "max_time"
    np.testing.assert_allclose(solution.time, [5.0, 7.0, 9.0, 11.0, 13.0, 15.0, 17.0, 19.0, 20.0])
    np.testing.assert_allclose(solution.temperature[:, 0], _relaxation(solution.time, 5.0, 100.0), rtol=1e-6)


def test_until_ends_the_phase_at_the_crossing_with_its_state_as_last_sample():
    # T = 100 K when exp(-t / 4) = 96 / 296, at t = 4 ln(296 / 96).
    crossing = _TAU * np.log(296.0 / 96.0)
    phase = Phase("precool", (_LINK,), _stage_below(100.0))
    model = _model((phase, _final_phase()))

    solution = integrate_phase(model, phase, 0.0, model.initial_temperature)

    assert solution.end == "until"
    np.testing.assert_allclose(solution.time, [0.0, 2.0, 4.0, crossing], rtol=1e-6)
    np.testing.assert_allclose(solution.temperature[-1], [100.0, 4.0], rtol=1e-6)


def test_stop_condition_ends_the_phase():
    phase = _final_phase()
    model = _model((phase,), stop_condition=_after(7.0))

    solution = integrate_phase(model, phase, 0.0, model.initial_temperature)

    assert solution.end == "stop_condition"
    np.testing.assert_allclose(solution.time[-1], 7.0)
    np.testing.assert_allclose(solution.temperature[-1, 0], _relaxation(7.0, 0.0, 300.0), rtol=1e-6)


@pytest.mark.parametrize(("until_at", "stop_at", "expected"), [(3.0, 9.0, "until"), (9.0, 3.0, "stop_condition")])
def test_earlier_of_until_and_stop_condition_ends_the_phase(until_at, stop_at, expected):
    phase = Phase("precool", (_LINK,), _after(until_at))
    model = _model((phase, _final_phase()), stop_condition=_after(stop_at))

    solution = integrate_phase(model, phase, 0.0, model.initial_temperature)

    assert solution.end == expected
    np.testing.assert_allclose(solution.time[-1], 3.0)


def test_integration_uses_only_the_given_phases_links():
    isolated = Phase("isolated", (), _after(10.0))
    model = _model((isolated, _final_phase()))

    solution = integrate_phase(model, isolated, 0.0, model.initial_temperature)

    np.testing.assert_array_equal(solution.temperature[:, 0], 300.0)


def test_solver_failure_is_reported_with_its_message():
    # dT/dt = T^2 / 2 from T = 300 K diverges at t = 2 / 300 s.
    runaway = Phase("runaway", (_Runaway("heater", ("stage",), 1.0),), None)
    model = _model((runaway,))

    solution = integrate_phase(model, runaway, 0.0, model.initial_temperature)

    assert solution.end == "solver_failure"
    assert solution.message
    assert solution.time[-1] < 2.0 / 300.0


@pytest.mark.parametrize("missing", ["max_time", "save_every"])
def test_missing_required_solver_option_is_a_configuration_error(missing):
    phase = _final_phase()
    options = {key: value for key, value in _OPTIONS.items() if key != missing}
    model = _model((phase,), options=options)

    with pytest.raises(ConfigurationError, match=missing):
        integrate_phase(model, phase, 0.0, model.initial_temperature)
