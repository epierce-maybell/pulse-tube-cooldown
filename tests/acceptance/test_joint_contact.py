"""Test physical invariants and known solutions for joint-contact links."""

from dataclasses import replace

import numpy as np
import pytest

from pulse_tube_cooldown.library.links.joint_contact import JointContact


@pytest.fixture
def link():
    """Return a joint with force 1000 N and quadratic conductance of 2 W/K at 10 K."""
    return JointContact(
        name="Flange",
        nodes=("Plate", "Bath"),
        force=1000.0,
        conductance_per_force=0.002,
        reference_temperature=10.0,
        temperature_exponent=2.0,
        source="Synthetic contact",
    )


@pytest.mark.parametrize("temperature", [(4.0, 300.0), (300.0, 4.0), (50.0, 77.0), (77.0, 77.0)])
def test_joint_contact_conserves_energy(link, temperature):
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
def test_joint_contact_adds_no_heat_capacity(link, temperature):
    """Verify that the joint-contact link contributes zero heat capacity."""
    terms = link.calc_energy_balance_terms(np.array(temperature))

    assert np.all(terms.heat_capacity == 0.0)


@pytest.mark.parametrize("force, scale", [(2000.0, 2.0), (500.0, 0.5)])
def test_heat_flow_scales_with_clamping_force(link, force, scale):
    """Verify that scaling clamping force scales heat flow by the same factor."""
    adjusted_link = replace(link, force=force)
    temperature = np.array([15.0, 25.0])

    original = link.calc_energy_balance_terms(temperature)
    adjusted = adjusted_link.calc_energy_balance_terms(temperature)

    np.testing.assert_allclose(adjusted.heat_flow, scale * original.heat_flow, rtol=1e-12, atol=0.0)


@pytest.mark.parametrize("conductance_per_force, scale", [(0.004, 2.0), (0.001, 0.5)])
def test_heat_flow_scales_with_conductance_per_force(link, conductance_per_force, scale):
    """Verify that scaling conductance per unit force scales heat flow by the same factor."""
    adjusted_link = replace(link, conductance_per_force=conductance_per_force)
    temperature = np.array([15.0, 25.0])

    original = link.calc_energy_balance_terms(temperature)
    adjusted = adjusted_link.calc_energy_balance_terms(temperature)

    np.testing.assert_allclose(adjusted.heat_flow, scale * original.heat_flow, rtol=1e-12, atol=0.0)


@pytest.mark.parametrize(
    "temperature, doubled_difference",
    [((18.0, 22.0), (16.0, 24.0)), ((22.0, 18.0), (24.0, 16.0))],
)
def test_doubling_temperature_difference_at_fixed_mean_doubles_heat_flow(link, temperature, doubled_difference):
    """Verify linear difference scaling at fixed mean under the mean-temperature approximation."""
    original = link.calc_energy_balance_terms(np.array(temperature))
    doubled = link.calc_energy_balance_terms(np.array(doubled_difference))

    np.testing.assert_allclose(doubled.heat_flow, 2.0 * original.heat_flow, rtol=1e-12, atol=0.0)


@pytest.mark.parametrize("temperature_exponent, scale", [(1.0, 2.0), (2.0, 4.0)])
@pytest.mark.parametrize(
    "temperature, doubled_mean",
    [((8.0, 12.0), (18.0, 22.0)), ((12.0, 8.0), (22.0, 18.0))],
)
def test_doubling_mean_at_fixed_difference_follows_temperature_power_law(
    link, temperature_exponent, scale, temperature, doubled_mean
):
    """Verify that doubling mean temperature multiplies flow by two or four for exponents one or two."""
    link = replace(link, temperature_exponent=temperature_exponent)

    original = link.calc_energy_balance_terms(np.array(temperature))
    doubled = link.calc_energy_balance_terms(np.array(doubled_mean))

    np.testing.assert_allclose(doubled.heat_flow, scale * original.heat_flow, rtol=1e-12, atol=0.0)


@pytest.mark.parametrize("temperature", [(4.0, 300.0), (300.0, 4.0), (8.0, 12.0), (18.0, 22.0)])
def test_changing_reference_for_same_conductance_curve_preserves_heat_flow(link, temperature):
    """Verify equivalent quadratic curves when reference temperature doubles and its coefficient quadruples."""
    rereferenced_link = replace(link, reference_temperature=20.0, conductance_per_force=0.008)
    temperature = np.array(temperature)

    original = link.calc_energy_balance_terms(temperature)
    rereferenced = rereferenced_link.calc_energy_balance_terms(temperature)

    np.testing.assert_allclose(rereferenced.heat_flow, original.heat_flow, rtol=1e-12, atol=0.0)


@pytest.mark.parametrize("temperature_exponent", [1.0, 1.5, 2.0])
@pytest.mark.parametrize("temperature, expected", [((8.5, 11.5), [6.0, -6.0]), ((11.5, 8.5), [-6.0, 6.0])])
def test_reference_temperature_conductance_sets_absolute_heat_flow(link, temperature_exponent, temperature, expected):
    """Verify that 2 W/K conductance at a mean of 10 K transfers 6 W across a difference of 3 K."""
    link = replace(link, temperature_exponent=temperature_exponent)

    terms = link.calc_energy_balance_terms(np.array(temperature))

    np.testing.assert_allclose(terms.heat_flow, expected, rtol=1e-12, atol=0.0, strict=True)
