"""Test an archetypal head fitted to one-compressor test-stand records."""

import tomllib
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest
import xarray as xr

from pulse_tube_cooldown import configure_model, simulate_cooldown
from pulse_tube_cooldown.engine.steady_state import solve_steady_state
from pulse_tube_cooldown.library.links.pulse_tube import (
    HeadCalibration,
    HeadCoefficients,
    LoadMap,
    build_pulse_tube_head,
    fit_pulse_tube_head,
    load_fit_parameters,
    load_head_design,
)
from tests.acceptance.fixtures.pulse_tube_cooler.compile_records import compile_records


@pytest.fixture(scope="session")
def objectives():
    """Return one-compressor measurement objectives and absolute tolerances."""
    return compile_records()[1]


@pytest.fixture(scope="session")
def model(objectives):
    """Return the unloaded test stand with a fitted head, initially at objectives.room_temperature."""
    config_path = Path(__file__).parent / "fixtures" / "pulse_tube_cooler" / "test_stand.toml"
    config = tomllib.loads(config_path.read_text(encoding="utf-8"))
    stand = configure_model(config_path)
    graph = stand.graphs[0]
    settings = config["calibration"]
    initial = HeadCalibration(
        coefficients=HeadCoefficients(**settings["coefficients"]),
        design=load_head_design(settings["design"]),
        load_map=LoadMap(
            temperature=np.array([objective.target for objective in objectives.load_map]).reshape(-1, 2),
            heater_power=np.array([objective.heater_power for objective in objectives.load_map[::2]]),
        ),
        test_stand=config["fridge"],
        test_stand_stages=tuple(settings["stages"]),
        warm_flange_temperature=objectives.room_temperature,
    )
    time = np.union1d(
        np.arange(0.0, max(objective.target for objective in objectives.cooldown), stand.solver_options["save_every"]),
        [objective.target for objective in objectives.cooldown],
    )
    traces = {}
    for stage in (1, 2):
        crossings = [objective for objective in objectives.cooldown if objective.stage == stage]
        traces[f"T{stage}"] = (
            "time",
            np.interp(
                time,
                [0.0, *(objective.target for objective in crossings)],
                [objectives.room_temperature, *(objective.temperature for objective in crossings)],
            ),
        )
    cooldown = xr.Dataset(traces, coords={"time": (time * 1e9).astype("datetime64[ns]")})
    defaults = load_fit_parameters()
    fit = settings["fit"]
    parameters = replace(
        defaults,
        sample_points=fit["sample_points"],
        solver={**defaults.solver, **fit["solver"]},
        optimizer={**defaults.optimizer, **fit["optimizer"]},
    )
    fitted = fit_pulse_tube_head(
        initial,
        cooldown,
        test_stand=graph,
        start=cooldown.time.values[0],
        end=cooldown.time.values[-1],
        parameters=parameters,
    )
    head = build_pulse_tube_head(
        fitted, graph, name="Head", stages=fitted.test_stand_stages, warm_flange_temperature=objectives.room_temperature
    )
    nodes = (*stand.nodes, *head.internal_nodes)
    phase = replace(stand.procedure.phases[0], links=(*graph.links, head))
    return replace(
        stand,
        nodes=nodes,
        procedure=replace(stand.procedure, phases=(phase,)),
        initial_temperature=np.full(len(nodes), objectives.room_temperature),
    )


@pytest.fixture(scope="module")
def load_map(model, result, objectives):
    """Return settled stage temperatures in K, keyed by the pair of stage heater powers in W."""
    graph = model.graphs[0]
    stages = [graph.get_node_index(name) for name in ("Stage 1", "Stage 2")]
    temperatures = {}
    temperature = result.data.temperature.values[-1]
    for powers in sorted({objective.heater_power for objective in objectives.load_map}):
        temperature = solve_steady_state(
            graph,
            temperature,
            fixed_boundary_nodes=model.fixed_boundary_nodes,
            heat_input={"Stage 1": powers[0], "Stage 2": powers[1]},
        )
        temperatures[powers] = temperature[stages]
    return temperatures


@pytest.fixture(scope="module")
def result(model):
    """Return the fitted head's unloaded test-stand cooldown."""
    return simulate_cooldown(model)


def test_load_map_matches_measurements(load_map, objectives):
    """Verify both settled stage temperatures at every heater-power pair against the measured ranges."""
    for objective in objectives.load_map:
        temperature = load_map[objective.heater_power][objective.stage - 1]

        assert temperature == pytest.approx(objective.target, rel=0.0, abs=objective.tolerance), objective


def test_stage_coupling_matches_measurements(load_map, objectives):
    """Verify all four cell-derived heat-flow responses in W/K against the measured ranges."""
    powers = np.array(sorted(load_map))
    low_low, low_high, high_low, high_high = [load_map[tuple(power)] for power in powers]
    temperature_response = np.column_stack(
        ((high_low - low_low + high_high - low_high) / 2, (low_high - low_low + high_high - high_low) / 2)
    )
    heat_flow_response = np.linalg.solve(temperature_response.T, -np.diag(np.ptp(powers, axis=0))).T

    for objective in objectives.coupling:
        response = heat_flow_response[objective.stage - 1, objective.temperature_stage - 1]

        assert response == pytest.approx(objective.target, rel=0.0, abs=objective.tolerance), objective


def test_cooldown_crossings_match_measurements(result, objectives):
    """Verify each stage's first downward temperature-crossing time in s against its measured range."""
    time = result.data.time.values
    for objective in objectives.cooldown:
        temperature = result.data.temperature.sel(node=f"Stage {objective.stage}").values
        crossings = np.flatnonzero(
            (temperature[:-1] > objective.temperature) & (temperature[1:] <= objective.temperature)
        )
        assert crossings.size > 0, objective
        index = crossings[0]
        crossing_time = np.interp(
            objective.temperature, temperature[index : index + 2][::-1], time[index : index + 2][::-1]
        )

        assert crossing_time == pytest.approx(objective.target, rel=0.0, abs=objective.tolerance), objective
