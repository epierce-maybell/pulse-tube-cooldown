"""Tests for integrating a model's phases."""

import dataclasses

import numpy as np
import pytest
from numpy.typing import NDArray

from pulse_tube_cooldown.core.errors import ConfigurationError
from pulse_tube_cooldown.core.types import Graph, Material, Model, Node, Phase, PhaseSolution, Procedure
from pulse_tube_cooldown.engine.solver import (
    build_result,
    integrate_phase,
    reconstruct_heat_flow,
    run_procedure,
    simulate_cooldown,
)

_CONSTANT = Material(
    name="Constant",
    density=1.0,
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
        meta={},
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


def test_two_phase_run_switches_at_until_and_continues_from_the_switch_state():
    # Precool relaxes with tau = 4 s until the stage reaches 100 K at t = 4 ln(296 / 96). The cooldown phase's
    # 1 W/K strap then relaxes it with tau = 2 s from 100 K at that time.
    crossing = _TAU * np.log(296.0 / 96.0)
    precool = Phase("precool", (_LINK,), _stage_below(100.0))
    cooldown = Phase("cooldown", (_Conductance("strap", ("stage", "bath"), 1.0),), None)
    model = _model((precool, cooldown))

    first, second = run_procedure(model)

    assert (first.phase, first.end) == (precool, "until")
    assert (second.phase, second.end) == (cooldown, "max_time")
    np.testing.assert_allclose(first.time[-1], crossing, rtol=1e-6)
    np.testing.assert_allclose(second.time[0], first.time[-1])
    np.testing.assert_allclose(second.temperature[0], first.temperature[-1])
    np.testing.assert_allclose(second.time[-1], 20.0)
    expected = 4.0 + 96.0 * np.exp(-(second.time - crossing) / 2.0)
    np.testing.assert_allclose(second.temperature[:, 0], expected, rtol=1e-6)


@pytest.mark.parametrize("threshold", [300.0, 400.0])
def test_phase_whose_until_is_met_at_its_start_ends_with_zero_length(threshold):
    precool = Phase("precool", (_LINK,), _stage_below(threshold))
    cooldown = _final_phase()
    model = _model((precool, cooldown))

    first, second = run_procedure(model)

    assert (first.phase, first.end) == (precool, "until")
    np.testing.assert_array_equal(first.time, [0.0])
    np.testing.assert_array_equal(first.temperature, [model.initial_temperature])
    assert (second.phase, second.end) == (cooldown, "max_time")
    assert second.time[0] == 0.0
    np.testing.assert_array_equal(second.temperature[0], model.initial_temperature)


def test_run_starts_at_zero_from_the_initial_temperature():
    model = _model((_final_phase(),))

    (solution,) = run_procedure(model)

    assert solution.time[0] == 0.0
    np.testing.assert_array_equal(solution.temperature[0], model.initial_temperature)


@pytest.mark.parametrize(
    ("first_phase", "stop_condition", "expected"),
    [
        (Phase("precool", (_LINK,), _after(10.0)), _after(3.0), "stop_condition"),
        (Phase("precool", (_LINK,), _after(50.0)), None, "max_time"),
        (Phase("runaway", (_Runaway("heater", ("stage",), 1.0),), _after(10.0)), None, "solver_failure"),
    ],
)
def test_run_ends_without_further_phases_unless_until_ends_the_phase(first_phase, stop_condition, expected):
    model = _model((first_phase, _final_phase()), stop_condition=stop_condition)

    solutions = run_procedure(model)

    assert [(solution.phase, solution.end) for solution in solutions] == [(first_phase, expected)]


_SAMPLES = np.array([[300.0, 4.0], [100.0, 4.0], [10.0, 4.0]])
_HEATER = _Runaway("heater", ("stage",), 1e-3)


def test_reconstructed_flows_match_each_link_at_each_sample():
    graph = Graph(nodes=(_STAGE, _BATH), links=(_LINK, _HEATER))

    flows = reconstruct_heat_flow(graph, _SAMPLES)

    # Strap: 0.5 W/K * (4 K - T) into the stage and the opposite into the bath. Heater: 1e-3 * T^2 into the stage.
    expected = [
        [[-148.0, 148.0], [90.0, np.nan]],
        [[-48.0, 48.0], [10.0, np.nan]],
        [[-3.0, 3.0], [0.1, np.nan]],
    ]
    np.testing.assert_allclose(flows, expected)


def test_phase_without_links_reconstructs_an_empty_link_axis():
    graph = Graph(nodes=(_STAGE, _BATH), links=())

    assert reconstruct_heat_flow(graph, _SAMPLES).shape == (3, 0, 2)


def _switching_run(end="max_time"):
    # Precool runs the 0.5 W/K strap and the heater; cooldown swaps in a 1 W/K strap and drops the heater.
    precool = Phase("precool", (_LINK, _HEATER), _stage_below(50.0))
    cooldown = Phase("cooldown", (_Conductance("strap", ("stage", "bath"), 1.0),), None)
    model = _model((precool, cooldown))
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
            end=end,
            time=np.array([3.0, 5.0]),
            temperature=np.array([[50.0, 4.0], [10.0, 4.0]]),
            message="cooldown ended",
        ),
    ]
    return model, solutions


def test_result_keeps_both_sides_of_a_phase_switch():
    model, solutions = _switching_run()

    data = build_result(model, solutions).data

    np.testing.assert_array_equal(data.time, [0.0, 2.0, 3.0, 3.0, 5.0])
    assert list(data.phase.values) == ["precool", "precool", "precool", "cooldown", "cooldown"]
    np.testing.assert_array_equal(data.temperature, [[300.0, 4.0], [100.0, 4.0], [50.0, 4.0], [50.0, 4.0], [10.0, 4.0]])
    # At 50 K: 0.5 W/K * (4 K - 50 K) before the switch and 1 W/K * (4 K - 50 K) after.
    np.testing.assert_allclose(data.heat_flow.isel(time=2).sel(link="strap"), [-23.0, 23.0])
    np.testing.assert_allclose(data.heat_flow.isel(time=3).sel(link="strap"), [-46.0, 46.0])


def test_result_flows_are_nan_where_a_link_is_inactive_or_untouched():
    model, solutions = _switching_run()

    heater = build_result(model, solutions).data.heat_flow.sel(link="heater")

    # Heater: 1e-3 * T^2 into the stage during precool; it never touches the bath.
    np.testing.assert_allclose(heater.sel(node="stage")[:3], [90.0, 10.0, 2.5])
    assert np.isnan(heater.sel(node="stage")[3:]).all()
    assert np.isnan(heater.sel(node="bath")).all()


def test_result_coordinates_name_nodes_and_every_procedure_link():
    model, solutions = _switching_run()

    data = build_result(model, solutions).data

    assert list(data.node.values) == ["stage", "bath"]
    assert list(data.link.values) == ["strap", "heater"]
    assert data.heat_flow.dims == ("time", "link", "node")
    np.testing.assert_array_equal(data.connected.sel(link="strap"), [True, True])
    np.testing.assert_array_equal(data.connected.sel(link="heater"), [True, False])


def test_result_includes_links_of_phases_that_never_ran():
    model, solutions = _switching_run()
    precool = dataclasses.replace(solutions[0], end="stop_condition")

    data = build_result(model, [precool]).data

    assert list(data.link.values) == ["strap", "heater"]
    assert list(data.phase.values) == ["precool"] * 3


@pytest.mark.parametrize("end", ["stop_condition", "max_time", "solver_failure"])
def test_result_status_is_how_the_final_phase_ended(end):
    model, solutions = _switching_run(end)

    result = build_result(model, solutions)

    assert result.status == end
    assert result.model is model
    assert result.data.attrs["status"] == end
    assert result.data.attrs["message"] == "cooldown ended"


def test_result_attributes_carry_model_meta():
    model, solutions = _switching_run()
    model = dataclasses.replace(model, meta={"experiment": "Toy switch", "version": "abc1234"})

    attrs = build_result(model, solutions).data.attrs

    assert attrs == {
        "experiment": "Toy switch",
        "version": "abc1234",
        "status": "max_time",
        "message": "cooldown ended",
    }


def test_result_records_units():
    model, solutions = _switching_run()

    data = build_result(model, solutions).data

    assert data.time.attrs["units"] == "s"
    assert data.temperature.attrs["units"] == "K"
    assert data.heat_flow.attrs["units"] == "W"


def test_result_flows_are_nan_throughout_a_phase_without_links():
    isolated = Phase("isolated", (), _after(1.0))
    model = _model((isolated, _final_phase()))
    solution = PhaseSolution(
        phase=isolated,
        end="stop_condition",
        time=np.array([0.0, 1.0]),
        temperature=np.array([[300.0, 4.0]] * 2),
        message="",
    )

    data = build_result(model, [solution]).data

    assert list(data.link.values) == ["strap"]
    assert np.isnan(data.heat_flow).all()


def _simulated_switch(stop_condition=None):
    # Precool's 0.5 W/K strap cools the stage with tau = 4 s until 100 K at t = 4 ln(296 / 96); cooldown's 1 W/K strap
    # then cools it with tau = 2 s from that state until 20 s.
    precool = Phase("precool", (_LINK,), _stage_below(100.0))
    cooldown = Phase("cooldown", (_Conductance("strap", ("stage", "bath"), 1.0),), None)
    return _model((precool, cooldown), stop_condition=stop_condition)


def test_simulation_follows_both_phases_and_keeps_the_switch():
    crossing = _TAU * np.log(296.0 / 96.0)
    model = _simulated_switch()

    result = simulate_cooldown(model)

    assert result.model is model
    assert result.status == "max_time"
    data = result.data
    precool_time = np.append(np.arange(0.0, crossing, 2.0), crossing)
    cooldown_time = np.append(np.arange(crossing, 20.0, 2.0), 20.0)
    np.testing.assert_allclose(data.time, np.concatenate([precool_time, cooldown_time]), rtol=1e-6)
    assert list(data.phase.values) == ["precool"] * precool_time.size + ["cooldown"] * cooldown_time.size
    expected = np.concatenate(
        [_relaxation(precool_time, 0.0, 300.0), 4.0 + 96.0 * np.exp(-(cooldown_time - crossing) / 2.0)]
    )
    np.testing.assert_allclose(data.temperature.sel(node="stage"), expected, rtol=1e-6)
    np.testing.assert_array_equal(data.temperature.sel(node="bath"), 4.0)


def test_simulation_heat_flows_follow_each_phases_strap():
    data = simulate_cooldown(_simulated_switch()).data

    conductance = np.where(data.phase.values == "precool", 0.5, 1.0)
    into_stage = conductance * (4.0 - data.temperature.sel(node="stage").values)
    np.testing.assert_allclose(data.heat_flow.sel(link="strap", node="stage"), into_stage, rtol=1e-6)
    np.testing.assert_allclose(data.heat_flow.sel(link="strap", node="bath"), -into_stage, rtol=1e-6)
    # At the 100 K switch, 0.5 W/K and then 1 W/K carry heat from the stage to the 4 K bath.
    switch = np.flatnonzero(np.isclose(data.time.values, _TAU * np.log(296.0 / 96.0)))
    np.testing.assert_allclose(data.heat_flow.isel(time=switch).sel(link="strap", node="stage"), [-48.0, -96.0])


def test_simulation_status_reports_a_procedure_stop():
    result = simulate_cooldown(_simulated_switch(stop_condition=_after(1.0)))

    assert result.status == "stop_condition"
    assert set(result.data.phase.values) == {"precool"}
    np.testing.assert_allclose(result.data.time[-1], 1.0)
