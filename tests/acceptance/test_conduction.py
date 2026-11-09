"""Test physical invariants and known solutions for conduction links."""

from dataclasses import replace

import numpy as np
import pytest

from pulse_tube_cooldown.core.types import Material
from pulse_tube_cooldown.library.links.conduction import Conduction


@pytest.fixture
def constant_material():
    """Return a material with constant conductivity of 2 W/(m K)."""
    return Material(
        name="Constant conductivity",
        density=1000.0,
        specific_heat=lambda temperature: np.full_like(temperature, 100.0),
        specific_heat_integral=lambda temperature: 100.0 * temperature,
        thermal_conductivity=lambda temperature: np.full_like(temperature, 2.0),
        conductivity_integral=lambda temperature: 2.0 * temperature,
    )


@pytest.fixture
def quadratic_material(constant_material):
    """Return a material with conductivity (T / 10)^2 W/(m K), for temperature T in K."""
    return replace(
        constant_material,
        name="Quadratic conductivity",
        thermal_conductivity=lambda temperature: (temperature / 10.0) ** 2,
        conductivity_integral=lambda temperature: temperature**3 / 300.0,
    )


@pytest.fixture
def link(quadratic_material):
    """Return a rod with area 3e-4 m^2 and length 0.1 m."""
    return Conduction(
        name="Rod",
        nodes=("Plate", "Bath"),
        material=quadratic_material,
        area=3.0e-4,
        length=0.1,
    )


@pytest.mark.parametrize("temperature", [(4.0, 300.0), (300.0, 4.0), (50.0, 77.0), (77.0, 77.0)])
def test_conduction_conserves_energy(link, temperature):
    """Verify that heat entering one node equals heat leaving the other."""
    terms = link.calc_energy_balance_terms(np.array(temperature))

    assert terms.heat_flow.sum() == pytest.approx(0.0, abs=1e-12)


@pytest.mark.parametrize("temperature", [4.0, 77.0, 300.0])
def test_equal_temperatures_carry_no_heat(link, temperature):
    """Verify that equal endpoint temperatures produce zero heat flow."""
    terms = link.calc_energy_balance_terms(np.array([temperature, temperature]))

    np.testing.assert_array_equal(terms.heat_flow, [0.0, 0.0])


@pytest.mark.parametrize("temperature", [(4.0, 300.0), (300.0, 4.0), (50.0, 77.0), (77.0, 50.0)])
def test_heat_flows_from_hot_to_cold(link, temperature):
    """Verify that the colder node gains heat and the hotter node loses heat."""
    temperature = np.array(temperature)

    terms = link.calc_energy_balance_terms(temperature)

    assert terms.heat_flow[np.argmin(temperature)] > 0.0
    assert terms.heat_flow[np.argmax(temperature)] < 0.0


def test_swapping_temperatures_swaps_heat_flows(link):
    """Verify that exchanging endpoint temperatures exchanges their heat flows."""
    forward = link.calc_energy_balance_terms(np.array([40.0, 250.0]))
    reverse = link.calc_energy_balance_terms(np.array([250.0, 40.0]))

    np.testing.assert_allclose(reverse.heat_flow, forward.heat_flow[::-1], rtol=1e-12, atol=0.0)


@pytest.mark.parametrize("area, length, scale", [(6.0e-4, 0.1, 2.0), (3.0e-4, 0.2, 0.5)])
def test_heat_flow_scales_with_area_over_length(link, area, length, scale):
    """Verify that doubling area doubles flow and doubling length halves flow."""
    resized_link = replace(link, area=area, length=length)
    temperature = np.array([10.0, 20.0])

    original = link.calc_energy_balance_terms(temperature)
    resized = resized_link.calc_energy_balance_terms(temperature)

    np.testing.assert_allclose(resized.heat_flow, scale * original.heat_flow, rtol=1e-12, atol=0.0)


@pytest.mark.parametrize("temperature", [(4.0, 300.0), (300.0, 4.0), (77.0, 77.0)])
def test_conduction_adds_no_heat_capacity(link, temperature):
    """Verify that the conduction link contributes zero heat capacity."""
    terms = link.calc_energy_balance_terms(np.array(temperature))

    assert np.all(terms.heat_capacity == 0.0)


@pytest.mark.parametrize("temperature, expected", [((10.0, 300.0), [1.74, -1.74]), ((300.0, 10.0), [-1.74, 1.74])])
def test_constant_conductivity_matches_fouriers_law(link, constant_material, temperature, expected):
    """Verify that constant conductivity gives the known Fourier heat flow in W."""
    link = replace(link, material=constant_material)

    terms = link.calc_energy_balance_terms(np.array(temperature))

    np.testing.assert_allclose(terms.heat_flow, expected, rtol=1e-12, atol=0.0, strict=True)


@pytest.mark.parametrize(
    "temperature, expected",
    [((10.0, 20.0), [0.07, -0.07]), ((20.0, 40.0), [0.56, -0.56]), ((40.0, 20.0), [-0.56, 0.56])],
)
def test_quadratic_conductivity_matches_integrated_fouriers_law(link, temperature, expected):
    """Verify that quadratic conductivity gives the known integrated heat flow in W."""
    terms = link.calc_energy_balance_terms(np.array(temperature))

    np.testing.assert_allclose(terms.heat_flow, expected, rtol=1e-12, atol=0.0, strict=True)


def test_heat_flow_is_independent_of_integral_reference(link, quadratic_material):
    """Verify that shifting the conductivity integral's reference preserves heat flow."""
    shifted_material = replace(
        quadratic_material,
        conductivity_integral=lambda temperature: temperature**3 / 300.0 + 1000.0,
    )
    shifted_link = replace(link, material=shifted_material)
    temperature = np.array([10.0, 20.0])

    original = link.calc_energy_balance_terms(temperature)
    shifted = shifted_link.calc_energy_balance_terms(temperature)

    np.testing.assert_allclose(shifted.heat_flow, original.heat_flow, rtol=1e-12, atol=0.0)
