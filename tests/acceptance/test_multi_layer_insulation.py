"""Test physical invariants and scaling relationships for multi-layer insulation links."""

from dataclasses import replace

import numpy as np
import pytest

from pulse_tube_cooldown.library.links.multi_layer_insulation import MultiLayerInsulation


@pytest.fixture
def combined_link():
    """Return a synthetic blanket with four layers at a density of 1000 layers/m."""
    return MultiLayerInsulation(
        name="Blanket",
        nodes=("Inner", "Outer"),
        area=0.1,
        layers=4,
        layer_density=1000.0,
        spacer_coefficient=1.0e-12,
        layer_density_exponent=2.56,
        radiation_coefficient=5.0e-10,
        emissivity=0.04,
        source="Synthetic blanket",
    )


@pytest.fixture
def spacer_link(combined_link):
    """Return a blanket carrying heat only through its spacers."""
    return replace(combined_link, radiation_coefficient=0.0)


@pytest.fixture
def radiation_link(combined_link):
    """Return a blanket carrying heat only through radiation."""
    return replace(combined_link, spacer_coefficient=0.0)


@pytest.fixture(params=["spacers", "radiation", "combined"])
def link(request, spacer_link, radiation_link, combined_link):
    """Return spacer-only, radiation-only, and combined blankets for shared invariants."""
    return {"spacers": spacer_link, "radiation": radiation_link, "combined": combined_link}[request.param]


@pytest.mark.parametrize("temperature", [(4.0, 300.0), (300.0, 4.0), (50.0, 77.0), (77.0, 77.0)])
def test_mli_conserves_energy(link, temperature):
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


@pytest.mark.parametrize("temperature", [(4.0, 300.0), (300.0, 4.0), (77.0, 77.0)])
def test_mli_adds_no_heat_capacity(link, temperature):
    """Verify that the MLI link contributes zero heat capacity."""
    terms = link.calc_energy_balance_terms(np.array(temperature))

    assert np.all(terms.heat_capacity == 0.0)


@pytest.mark.parametrize("area, scale", [(0.2, 2.0), (0.05, 0.5)])
def test_heat_flow_scales_with_blanket_area(link, area, scale):
    """Verify that scaling blanket area scales heat flow by the same factor."""
    resized_link = replace(link, area=area)
    temperature = np.array([50.0, 300.0])

    original = link.calc_energy_balance_terms(temperature)
    resized = resized_link.calc_energy_balance_terms(temperature)

    np.testing.assert_allclose(resized.heat_flow, scale * original.heat_flow, rtol=1e-12, atol=0.0)


@pytest.mark.parametrize("temperature", [(50.0, 300.0), (300.0, 50.0)])
def test_combined_flow_equals_sum_of_spacer_and_radiation_flows(
    combined_link, spacer_link, radiation_link, temperature
):
    """Verify that spacer conduction and radiation contribute additively."""
    temperature = np.array(temperature)

    combined = combined_link.calc_energy_balance_terms(temperature)
    spacers = spacer_link.calc_energy_balance_terms(temperature)
    radiation = radiation_link.calc_energy_balance_terms(temperature)

    np.testing.assert_allclose(combined.heat_flow, spacers.heat_flow + radiation.heat_flow, rtol=1e-12, atol=0.0)


@pytest.mark.parametrize("spacer_coefficient, scale", [(2.0e-12, 2.0), (5.0e-13, 0.5)])
def test_spacer_flow_scales_with_spacer_coefficient(spacer_link, spacer_coefficient, scale):
    """Verify linear scaling of spacer-only flow with its coefficient."""
    adjusted_link = replace(spacer_link, spacer_coefficient=spacer_coefficient)
    temperature = np.array([50.0, 300.0])

    original = spacer_link.calc_energy_balance_terms(temperature)
    adjusted = adjusted_link.calc_energy_balance_terms(temperature)

    np.testing.assert_allclose(adjusted.heat_flow, scale * original.heat_flow, rtol=1e-12, atol=0.0)


@pytest.mark.parametrize("radiation_coefficient, scale", [(1.0e-9, 2.0), (2.5e-10, 0.5)])
def test_radiation_flow_scales_with_radiation_coefficient(radiation_link, radiation_coefficient, scale):
    """Verify linear scaling of radiation-only flow with its coefficient."""
    adjusted_link = replace(radiation_link, radiation_coefficient=radiation_coefficient)
    temperature = np.array([50.0, 300.0])

    original = radiation_link.calc_energy_balance_terms(temperature)
    adjusted = adjusted_link.calc_energy_balance_terms(temperature)

    np.testing.assert_allclose(adjusted.heat_flow, scale * original.heat_flow, rtol=1e-12, atol=0.0)


@pytest.mark.parametrize("emissivity, scale", [(0.08, 2.0), (0.02, 0.5)])
def test_radiation_flow_scales_with_layer_emissivity(radiation_link, emissivity, scale):
    """Verify linear scaling of radiation-only flow with layer emissivity."""
    adjusted_link = replace(radiation_link, emissivity=emissivity)
    temperature = np.array([50.0, 300.0])

    original = radiation_link.calc_energy_balance_terms(temperature)
    adjusted = adjusted_link.calc_energy_balance_terms(temperature)

    np.testing.assert_allclose(adjusted.heat_flow, scale * original.heat_flow, rtol=1e-12, atol=0.0)


@pytest.mark.parametrize("layer_density_exponent, scale", [(1.0, 2.0), (2.0, 4.0), (2.56, 2.0**2.56)])
def test_doubling_density_at_fixed_layer_count_follows_spacer_power_law(spacer_link, layer_density_exponent, scale):
    """Verify spacer-flow scaling with density for the configured exponent at fixed layer count."""
    spacer_link = replace(spacer_link, layer_density_exponent=layer_density_exponent)
    denser_link = replace(spacer_link, layer_density=2000.0)
    temperature = np.array([50.0, 300.0])

    original = spacer_link.calc_energy_balance_terms(temperature)
    denser = denser_link.calc_energy_balance_terms(temperature)

    np.testing.assert_allclose(denser.heat_flow, scale * original.heat_flow, rtol=1e-12, atol=0.0)


@pytest.mark.parametrize("layer_density", [500.0, 2000.0])
def test_radiation_flow_is_independent_of_density_at_fixed_layer_count(radiation_link, layer_density):
    """Verify that changing layer density at fixed layer count preserves radiation-only flow."""
    adjusted_link = replace(radiation_link, layer_density=layer_density)
    temperature = np.array([50.0, 300.0])

    original = radiation_link.calc_energy_balance_terms(temperature)
    adjusted = adjusted_link.calc_energy_balance_terms(temperature)

    np.testing.assert_allclose(adjusted.heat_flow, original.heat_flow, rtol=1e-12, atol=0.0)


@pytest.mark.parametrize("temperature", [(50.0, 300.0), (300.0, 50.0)])
def test_four_to_nine_layers_at_fixed_density_halves_spacer_flow(spacer_link, temperature):
    """Verify that increasing layers from four to nine halves spacer-only flow at fixed density."""
    thicker_link = replace(spacer_link, layers=9)
    temperature = np.array(temperature)

    original = spacer_link.calc_energy_balance_terms(temperature)
    thicker = thicker_link.calc_energy_balance_terms(temperature)

    np.testing.assert_allclose(thicker.heat_flow, 0.5 * original.heat_flow, rtol=1e-12, atol=0.0)


@pytest.mark.parametrize("temperature", [(50.0, 300.0), (300.0, 50.0)])
def test_doubling_layers_at_fixed_density_halves_radiation_flow(radiation_link, temperature):
    """Verify that increasing layers from four to eight halves radiation-only flow at fixed density."""
    thicker_link = replace(radiation_link, layers=8)
    temperature = np.array(temperature)

    original = radiation_link.calc_energy_balance_terms(temperature)
    thicker = thicker_link.calc_energy_balance_terms(temperature)

    np.testing.assert_allclose(thicker.heat_flow, 0.5 * original.heat_flow, rtol=1e-12, atol=0.0)


@pytest.mark.parametrize("temperature", [(40.0, 100.0), (100.0, 40.0)])
def test_doubling_temperatures_quadruples_spacer_flow(spacer_link, temperature):
    """Verify quadratic absolute-temperature scaling of spacer-only flow."""
    temperature = np.array(temperature)

    original = spacer_link.calc_energy_balance_terms(temperature)
    doubled = spacer_link.calc_energy_balance_terms(2.0 * temperature)

    np.testing.assert_allclose(doubled.heat_flow, 4.0 * original.heat_flow, rtol=1e-12, atol=0.0)


@pytest.mark.parametrize("temperature", [(40.0, 100.0), (100.0, 40.0)])
def test_doubling_temperatures_scales_radiation_flow_with_empirical_exponent(radiation_link, temperature):
    """Verify absolute-temperature scaling of radiation-only flow with the empirical exponent 4.67."""
    temperature = np.array(temperature)

    original = radiation_link.calc_energy_balance_terms(temperature)
    doubled = radiation_link.calc_energy_balance_terms(2.0 * temperature)

    np.testing.assert_allclose(doubled.heat_flow, 2.0**4.67 * original.heat_flow, rtol=1e-12, atol=0.0)
