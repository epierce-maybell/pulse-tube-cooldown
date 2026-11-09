"""Test physical invariants and known solutions for radiation links."""

from dataclasses import replace

import numpy as np
import pytest

from pulse_tube_cooldown.core.types import Surface
from pulse_tube_cooldown.library.links.radiation import Radiation


def constant_surface(emissivity):
    """Return a surface with constant emissivity."""
    return Surface(
        name=f"Constant emissivity {emissivity}",
        emissivity=lambda temperature: np.full_like(temperature, emissivity),
    )


@pytest.fixture
def black_surface():
    """Return a black surface with emissivity one."""
    return constant_surface(1.0)


@pytest.fixture
def link():
    """Return enclosed surfaces with areas 0.1 and 0.2 m^2 and emissivities 0.5 and 0.25."""
    return Radiation(
        name="Shield",
        nodes=("Inner", "Outer"),
        area=(0.1, 0.2),
        surface=(constant_surface(0.5), constant_surface(0.25)),
    )


@pytest.mark.parametrize("temperature", [(4.0, 300.0), (300.0, 4.0), (50.0, 77.0), (77.0, 77.0)])
def test_radiation_conserves_energy(link, temperature):
    """Verify that heat entering one node equals heat leaving the other."""
    terms = link.calc_energy_balance_terms(np.array(temperature))

    assert terms.heat_flow.sum() == pytest.approx(0.0, abs=1e-12)


@pytest.mark.parametrize("temperature", [4.0, 77.0, 300.0])
def test_equal_temperatures_carry_no_heat(link, temperature):
    """Verify that equal temperatures produce zero flow for different areas and finishes."""
    terms = link.calc_energy_balance_terms(np.array([temperature, temperature]))

    np.testing.assert_array_equal(terms.heat_flow, [0.0, 0.0])


@pytest.mark.parametrize("temperature", [(4.0, 300.0), (300.0, 4.0), (50.0, 77.0), (77.0, 50.0)])
def test_heat_flows_from_hot_to_cold(link, temperature):
    """Verify that the colder node gains heat and the hotter node loses heat."""
    temperature = np.array(temperature)

    terms = link.calc_energy_balance_terms(temperature)

    assert terms.heat_flow[np.argmin(temperature)] > 0.0
    assert terms.heat_flow[np.argmax(temperature)] < 0.0


@pytest.mark.parametrize("temperature", [(4.0, 300.0), (300.0, 4.0), (77.0, 77.0)])
def test_radiation_adds_no_heat_capacity(link, temperature):
    """Verify that the radiation link contributes zero heat capacity."""
    terms = link.calc_energy_balance_terms(np.array(temperature))

    assert np.all(terms.heat_capacity == 0.0)


def test_swapping_temperatures_swaps_heat_flows_for_constant_emissivities(link):
    """Verify that constant emissivities preserve flow magnitude under temperature reversal."""
    forward = link.calc_energy_balance_terms(np.array([40.0, 250.0]))
    reverse = link.calc_energy_balance_terms(np.array([250.0, 40.0]))

    np.testing.assert_allclose(reverse.heat_flow, forward.heat_flow[::-1], rtol=1e-12, atol=0.0)


@pytest.mark.parametrize("area, scale", [((0.2, 0.4), 2.0), ((0.05, 0.1), 0.5)])
def test_heat_flow_scales_with_both_areas(link, area, scale):
    """Verify that scaling both areas by the same factor scales heat flow by that factor."""
    resized_link = replace(link, area=area)
    temperature = np.array([40.0, 250.0])

    original = link.calc_energy_balance_terms(temperature)
    resized = resized_link.calc_energy_balance_terms(temperature)

    np.testing.assert_allclose(resized.heat_flow, scale * original.heat_flow, rtol=1e-12, atol=0.0)


@pytest.mark.parametrize("temperature", [(40.0, 100.0), (100.0, 40.0)])
def test_doubling_temperatures_multiplies_flow_by_sixteen_for_constant_emissivities(link, temperature):
    """Verify fourth-power temperature scaling with constant emissivities."""
    temperature = np.array(temperature)

    original = link.calc_energy_balance_terms(temperature)
    doubled = link.calc_energy_balance_terms(2.0 * temperature)

    np.testing.assert_allclose(doubled.heat_flow, 16.0 * original.heat_flow, rtol=1e-12, atol=0.0)


@pytest.mark.parametrize("emissivity", [(0.75, 0.25), (0.5, 0.5)])
@pytest.mark.parametrize("temperature", [(40.0, 250.0), (250.0, 40.0)])
def test_increasing_either_emissivity_increases_heat_exchange(link, emissivity, temperature):
    """Verify that increasing either surface's emissivity increases flow magnitude."""
    higher_emissivity_link = replace(link, surface=tuple(constant_surface(value) for value in emissivity))
    temperature = np.array(temperature)

    original = link.calc_energy_balance_terms(temperature)
    increased = higher_emissivity_link.calc_energy_balance_terms(temperature)

    assert np.all(np.abs(increased.heat_flow) > np.abs(original.heat_flow))


@pytest.mark.parametrize(
    "temperature, expected",
    [
        ((100.0, 200.0), [8.505561628776644, -8.505561628776644]),
        ((200.0, 100.0), [-8.505561628776644, 8.505561628776644]),
    ],
)
def test_black_surfaces_match_stefan_boltzmann_law(link, black_surface, temperature, expected):
    """Verify the solution matches independently-derived result from black-body radiation law."""
    black_link = replace(link, surface=(black_surface, black_surface))

    terms = black_link.calc_energy_balance_terms(np.array(temperature))

    np.testing.assert_allclose(terms.heat_flow, expected, rtol=1e-12, atol=0.0, strict=True)


@pytest.mark.parametrize("outer_area", [0.1, 0.4, 1.0])
def test_black_surface_exchange_is_independent_of_outer_area(link, black_surface, outer_area):
    """Verify that changing a black enclosure's area preserves heat flow at fixed inner area."""
    black_link = replace(link, surface=(black_surface, black_surface))
    resized_link = replace(black_link, area=(0.1, outer_area))
    temperature = np.array([40.0, 250.0])

    original = black_link.calc_energy_balance_terms(temperature)
    resized = resized_link.calc_energy_balance_terms(temperature)

    np.testing.assert_allclose(resized.heat_flow, original.heat_flow, rtol=1e-12, atol=0.0)


@pytest.mark.parametrize("temperature", [(40.0, 250.0), (250.0, 40.0)])
def test_equal_area_half_emissivity_surfaces_exchange_one_third_of_black_body_flow(link, black_surface, temperature):
    """Verify that equal-area surfaces with emissivity 0.5 exchange one third of black-body flow."""
    gray_link = replace(link, area=(0.1, 0.1), surface=(link.surface[0], link.surface[0]))
    black_link = replace(gray_link, surface=(black_surface, black_surface))
    temperature = np.array(temperature)

    gray = gray_link.calc_energy_balance_terms(temperature)
    black = black_link.calc_energy_balance_terms(temperature)

    np.testing.assert_allclose(gray.heat_flow, black.heat_flow / 3.0, rtol=1e-12, atol=0.0)


@pytest.mark.parametrize("temperature", [(40.0, 250.0), (250.0, 40.0)])
def test_unequal_area_gray_surfaces_exchange_two_sevenths_of_black_body_flow(link, black_surface, temperature):
    """Verify the flow ratio for area ratio 1:2 and inner/outer emissivities 0.5/0.25."""
    black_link = replace(link, surface=(black_surface, black_surface))
    temperature = np.array(temperature)

    gray = link.calc_energy_balance_terms(temperature)
    black = black_link.calc_energy_balance_terms(temperature)

    np.testing.assert_allclose(gray.heat_flow, 2.0 * black.heat_flow / 7.0, rtol=1e-12, atol=0.0)


@pytest.mark.parametrize(
    "temperature, emissivity",
    [((50.0, 200.0), (0.05, 0.3)), ((200.0, 50.0), (0.2, 0.45))],
)
def test_temperature_dependent_finishes_match_constants_at_each_surface_temperature(link, temperature, emissivity):
    """Verify equivalence to constant finishes evaluated at each surface's own temperature."""
    varying_link = replace(
        link,
        surface=(
            Surface(name="Rising emissivity", emissivity=lambda temperature: temperature / 1000.0),
            Surface(name="Falling emissivity", emissivity=lambda temperature: 0.5 - temperature / 1000.0),
        ),
    )
    constant_link = replace(link, surface=tuple(constant_surface(value) for value in emissivity))
    temperature = np.array(temperature)

    varying = varying_link.calc_energy_balance_terms(temperature)
    constant = constant_link.calc_energy_balance_terms(temperature)

    np.testing.assert_allclose(varying.heat_flow, constant.heat_flow, rtol=1e-12, atol=0.0)
