"""Test physical invariants of a full cooldown simulation, for any fridge and configuration."""

import tomllib
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from pulse_tube_cooldown import configure_model, simulate_cooldown
from pulse_tube_cooldown.core.types import Model, Node, Phase, Result
from pulse_tube_cooldown.engine.steady_state import solve_steady_state

pytestmark = [pytest.mark.slow, pytest.mark.filterwarnings("error")]


@pytest.fixture(scope="module")
def model():
    """Return the synthetic fridge, which uses every link type."""
    return configure_model(Path(__file__).parent / "fixtures" / "synthetic_fridge.toml")


@pytest.fixture(scope="module")
def result(model):
    """Return the synthetic fridge's cooldown."""
    return simulate_cooldown(model)


@pytest.fixture(scope="module")
def passive_model(model):
    """Return the synthetic fridge with only links that carry heat from hot to cold, starting between 30 and 295 K."""
    with (Path(__file__).parent / "fixtures" / "synthetic_fridge.toml").open("rb") as file:
        links = tomllib.load(file)["fridge"]["links"]
    passive_types = {
        "conduction",
        "exchange_gas",
        "gas_gap_switch",
        "joint_contact",
        "multi_layer_insulation",
        "radiation",
    }
    passive = {link["name"] for link in links if link["type"] in passive_types}
    phases = tuple(
        Phase(phase.name, tuple(link for link in phase.links if link.name in passive), phase.until)
        for phase in model.procedure.phases
    )
    return replace(
        model,
        procedure=replace(model.procedure, phases=phases),
        initial_temperature=np.linspace(30.0, 295.0, len(model.nodes)),
    )


@pytest.fixture(scope="module")
def passive_result(passive_model):
    """Return the passive fridge's cooldown."""
    return simulate_cooldown(passive_model)


def net_heat_flow(result: Result) -> np.ndarray:
    """Return the net heat flow in W into each node at each sample, indexed by time and node."""
    return result.data.heat_flow.sum("link", skipna=True).values


def free_nodes(model: Model) -> np.ndarray:
    """Return the indices of nodes not held at their initial temperatures."""
    return np.setdiff1d(np.arange(len(model.nodes)), model.fixed_boundary_nodes)


def stored_energy(model: Model, temperature: np.ndarray) -> float:
    """Return the energy in J stored in the free nodes' materials, from a fixed reference temperature."""
    return sum(
        mass * float(material.specific_heat_integral(np.array([temperature[index]]))[0])
        for index in free_nodes(model)
        for material, mass in model.nodes[index].composition
    )


def scale_masses(model: Model, factor: float) -> Model:
    """Return the model with every node's material masses multiplied by `factor`."""
    nodes = tuple(
        Node(node.name, tuple((material, factor * mass) for material, mass in node.composition)) for node in model.nodes
    )
    return replace(model, nodes=nodes)


def test_uniform_temperature_stays_put(passive_model):
    """Verify that a passive fridge with every node at 77 K stays at 77 K."""
    uniform = replace(passive_model, initial_temperature=np.full(len(passive_model.nodes), 77.0))

    result = simulate_cooldown(uniform)

    np.testing.assert_allclose(result.data.temperature.values, 77.0, rtol=0.0, atol=1e-9)


def test_total_entropy_never_decreases(passive_result):
    """Verify that the entropy production rate, the sum of net heat flow over temperature, is never negative."""
    entropy_rate = net_heat_flow(passive_result) / passive_result.data.temperature.values
    scale = np.abs(entropy_rate).sum(axis=1)

    assert np.all(entropy_rate.sum(axis=1) >= -1e-9 * scale)


def test_temperatures_stay_within_initial_and_boundary_extremes(passive_result):
    """Verify that no node leaves the range of the initial and boundary temperatures."""
    initial = passive_result.model.initial_temperature
    temperature = passive_result.data.temperature.values

    assert temperature.min() >= initial.min() * (1 - 1e-9)
    assert temperature.max() <= initial.max() * (1 + 1e-9)


def test_cooldown_conserves_energy(model):
    """Verify that the change in stored energy equals the time-integrated net heat flow into the free nodes."""
    sampled = replace(model, solver_options={**model.solver_options, "save_every": 10.0})
    result = simulate_cooldown(sampled)
    temperature = result.data.temperature.values
    heat_flow = net_heat_flow(result)[:, free_nodes(model)].sum(axis=1)

    energy_change = stored_energy(model, temperature[-1]) - stored_energy(model, temperature[0])
    heat_in = np.trapezoid(heat_flow, result.data.time.values)

    assert heat_in == pytest.approx(energy_change, rel=1e-3)


def test_long_cooldown_reaches_the_steady_state(model):
    """Verify that a cooldown run without a stop condition ends at the steady state of its final phase."""
    unstopped = replace(
        model,
        procedure=replace(model.procedure, stop_condition=None),
        solver_options={**model.solver_options, "max_time": 1.0e6},
    )
    result = simulate_cooldown(unstopped)

    steady = solve_steady_state(
        unstopped.graphs[-1], model.initial_temperature, fixed_boundary_nodes=model.fixed_boundary_nodes
    )

    np.testing.assert_allclose(result.data.temperature.values[-1], steady, rtol=1e-3)


def test_more_thermal_mass_takes_longer_to_cool(model):
    """Verify that cooldown time increases as every node's mass is scaled by 0.5, 1, and 2."""
    times = [float(simulate_cooldown(scale_masses(model, factor)).data.time[-1]) for factor in (0.5, 1.0, 2.0)]

    assert times[0] < times[1] < times[2]


def test_tighter_solver_tolerance_preserves_the_cooldown(model, result):
    """Verify that tightening the solver tolerances by 1000 changes cooldown time and final temperatures by < 0.1%."""
    tight = replace(model, solver_options={**model.solver_options, "rtol": 1e-9, "atol": 1e-11})

    tight_result = simulate_cooldown(tight)

    assert float(tight_result.data.time[-1]) == pytest.approx(float(result.data.time[-1]), rel=1e-3)
    np.testing.assert_allclose(tight_result.data.temperature.values[-1], result.data.temperature.values[-1], rtol=1e-3)
