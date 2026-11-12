"""Test scalar solver traversal of each packaged material property."""

import tomllib
from pathlib import Path

import numpy as np
import pytest
from scipy.integrate import quad, solve_ivp

from pulse_tube_cooldown.library import materials
from pulse_tube_cooldown.library.materials.load import load_material


@pytest.fixture(
    params=sorted((Path(materials.__file__).parent / "data").glob("*.toml")),
    ids=lambda path: path.stem,
)
def material_entry(request):
    """Return a packaged material entry, parametrized by its data file."""
    return tomllib.loads(request.param.read_text(encoding="utf-8"))


@pytest.fixture
def material(material_entry):
    """Return the material through its normal loader."""
    return load_material(material_entry["name"])


def _cool(rate, temperatures, max_time):
    """Integrate a scalar cooling rate and record downward temperature crossings."""
    events = [lambda time, temperature, target=target: temperature[0] - target for target in temperatures[1:]]
    for event in events:
        event.direction = -1.0
    events[-1].terminal = True
    return solve_ivp(
        rate,
        (0.0, max_time),
        [temperatures[0]],
        method="Radau",
        events=events,
        rtol=1e-8,
        atol=1e-10,
    )


def test_specific_heat_traversal_matches_constant_power_cooling(material, material_entry):
    """Verify constant-power cooling across the material's specific-heat range.

    The reference time is removed energy / cooling power. The solver integrates
    temperature at a rate of -cooling power / heat capacity.
    """
    mass = 1.0  # kg
    cooling_power = 1.0  # W
    spec = material_entry["specific_heat"]
    final_temperature, initial_temperature = spec["range"]
    table_temperatures = [row[0] for row in spec.get("table", spec.get("measured", []))]
    temperatures = np.union1d(np.geomspace(final_temperature, initial_temperature, 13), table_temperatures)
    temperatures = temperatures[(final_temperature <= temperatures) & (temperatures <= initial_temperature)][::-1]

    expected_times = []
    for hot, cold in zip(temperatures[:-1], temperatures[1:], strict=True):
        specific_energy, _ = quad(
            lambda temperature: material.specific_heat(np.array(temperature)).item(),
            cold,
            hot,
            epsabs=1e-10,
            epsrel=1e-10,
        )
        expected_times.append(mass * specific_energy / cooling_power)

    def cooling_rate(time, temperature):
        heat_capacity = mass * material.specific_heat(np.maximum(temperature, final_temperature))
        return -cooling_power / heat_capacity

    solution = _cool(cooling_rate, temperatures, max_time=1.1 * sum(expected_times))

    assert solution.success, solution.message
    assert np.all(np.isfinite(solution.y))
    assert np.all(np.diff(solution.y[0]) < 0.0)
    assert all(len(times) == 1 for times in solution.t_events)
    elapsed = np.diff(np.concatenate(([0.0], *solution.t_events)))
    np.testing.assert_allclose(elapsed, expected_times, rtol=1e-4, atol=1e-8)


def test_conductivity_traversal_matches_cooling_to_a_fixed_bath(material, material_entry):
    """Verify cooling of a constant-heat-capacity body through a rod to a fixed bath.

    The solver obtains heat flow from the loaded conductivity integral. The
    reference integrates conductivity independently, then integrates heat
    capacity / heat flow to obtain the time between temperature checkpoints.
    """
    if material.thermal_conductivity is None:
        pytest.skip("Material supplies no thermal conductivity.")

    heat_capacity = 1.0  # J/K
    rod_area = 1e-4  # m^2
    rod_length = 0.1  # m
    spec = material_entry["thermal_conductivity"]
    bath_temperature, initial_temperature = spec["range"]
    final_temperature = bath_temperature + 1.0  # K; equilibrium at the bath takes infinite time.
    knots = np.array([row[0] for row in spec.get("table", spec.get("measured", []))])
    temperatures = np.union1d(np.geomspace(final_temperature, initial_temperature, 13), knots)
    temperatures = temperatures[(final_temperature <= temperatures) & (temperatures <= initial_temperature)][::-1]

    def reference_heat_flow(temperature):
        integrated_conductivity, _ = quad(
            lambda temperature: material.thermal_conductivity(np.array(temperature)).item(),
            bath_temperature,
            temperature,
            points=knots[(bath_temperature < knots) & (knots < temperature)],
            limit=max(100, len(knots) + 1),
            epsabs=1e-10,
            epsrel=1e-10,
        )
        return rod_area / rod_length * integrated_conductivity

    expected_times = []
    for hot, cold in zip(temperatures[:-1], temperatures[1:], strict=True):
        duration, _ = quad(
            lambda temperature: heat_capacity / reference_heat_flow(temperature),
            cold,
            hot,
            epsabs=1e-10,
            epsrel=1e-10,
        )
        expected_times.append(duration)

    bath_integral = material.conductivity_integral(np.array(bath_temperature))

    def cooling_rate(time, temperature):
        integrated_conductivity = material.conductivity_integral(temperature) - bath_integral
        heat_flow = rod_area / rod_length * integrated_conductivity
        return -heat_flow / heat_capacity

    solution = _cool(cooling_rate, temperatures, max_time=1.1 * sum(expected_times))

    assert solution.success, solution.message
    assert np.all(np.isfinite(solution.y))
    assert np.all(np.diff(solution.y[0]) < 0.0)
    assert all(len(times) == 1 for times in solution.t_events)
    elapsed = np.diff(np.concatenate(([0.0], *solution.t_events)))
    np.testing.assert_allclose(elapsed, expected_times, rtol=1e-4, atol=1e-8)
