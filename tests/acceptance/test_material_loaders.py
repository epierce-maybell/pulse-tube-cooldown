"""Test interpolation shape, continuity, and consistency of loaded material integrals."""

import numpy as np
import pytest
from scipy.integrate import quad

from pulse_tube_cooldown.library.materials.load import load_material
from pulse_tube_cooldown.library.property_fits import make_fit_from_table, make_measured_with_lattice_fit


@pytest.fixture
def table():
    """Return unevenly spaced points with a plateau, a narrow peak, and a trough."""
    return np.array([[1.0, 1.0], [2.0, 1.0], [3.0, 2.0], [4.0, 6.0], [4.001, 3.0], [6.0, 2.9], [7.0, 4.0]])


@pytest.fixture
def interpolant(table):
    """Return the interpolant through the synthetic table."""
    return make_fit_from_table(table)


@pytest.fixture
def material(tmp_path):
    """Return a material with a narrow specific-heat peak and a steep conductivity rise."""
    (tmp_path / "synthetic.toml").write_text(
        """
        name = "Synthetic material"
        density = 1000.0

        [specific_heat]
        form = "table"
        range = [1.0, 7.0]
        table = [[1.0, 1.0], [2.0, 1.0], [3.0, 2.0], [4.0, 6.0], [4.001, 3.0], [6.0, 2.9], [7.0, 4.0]]

        [thermal_conductivity]
        form = "table"
        range = [1.0, 7.0]
        table = [[1.0, 2.0], [2.0, 2.0], [3.0, 1.5], [4.0, 1.0], [4.001, 5.0], [6.0, 3.0], [7.0, 3.0]]
        """,
        encoding="utf-8",
    )
    return load_material("Synthetic material", directory=tmp_path)


@pytest.fixture(params=["specific_heat", "thermal_conductivity"])
def property_functions(request, material):
    """Return a loaded property and its temperature integral."""
    return {
        "specific_heat": (material.specific_heat, material.specific_heat_integral),
        "thermal_conductivity": (material.thermal_conductivity, material.conductivity_integral),
    }[request.param]


def test_interpolation_recovers_supplied_points(interpolant, table):
    """Verify that interpolation passes through the supplied data."""
    value = interpolant(table[:, 0])

    np.testing.assert_allclose(value, table[:, 1], rtol=1e-12, atol=0.0)


def test_interpolation_introduces_no_overshoot_or_undershoot(interpolant, table):
    """Verify that each interval stays within its endpoint values around sharp features."""
    for (start, first), (end, last) in zip(table[:-1], table[1:], strict=True):
        value = interpolant(np.linspace(start, end, 101))

        assert np.all(value >= min(first, last) - 1e-12), (start, end)
        assert np.all(value <= max(first, last) + 1e-12), (start, end)


def test_interpolation_preserves_each_intervals_direction(interpolant, table):
    """Verify that interpolation introduces no reversals between increasing or decreasing points."""
    for (start, first), (end, last) in zip(table[:-1], table[1:], strict=True):
        value = interpolant(np.linspace(start, end, 101))

        assert np.all(np.sign(last - first) * np.diff(value) >= -1e-12), (start, end)


def test_interpolation_preserves_flat_intervals(interpolant):
    """Verify that a plateau remains flat beside an increasing interval."""
    value = interpolant(np.linspace(1.0, 2.0, 101))

    np.testing.assert_allclose(value, 1.0, rtol=1e-12, atol=0.0)


def test_interpolation_has_continuous_slopes_at_interior_points(interpolant, table):
    """Verify matching one-sided slopes at knots, including the narrow peak and trough."""
    temperature, value = table.T
    spacing = np.diff(temperature)
    secant = np.abs(np.diff(value) / spacing)
    step = 1e-6 * np.minimum(spacing[:-1], spacing[1:])
    knots = temperature[1:-1]

    below = (interpolant(knots) - interpolant(knots - step)) / step
    above = (interpolant(knots + step) - interpolant(knots)) / step

    assert np.all(np.abs(below - above) <= 1e-4 * np.maximum(secant[:-1], secant[1:]))


@pytest.mark.parametrize("start, end", [(3.0, 4.0), (4.0, 4.001), (4.001, 6.0)])
def test_integral_resolves_sharp_and_broad_features(property_functions, start, end):
    """Verify integral differences against independent quadrature across each feature."""
    property_function, integral = property_functions
    expected, _ = quad(lambda t: property_function(np.array(t)).item(), start, end, epsabs=1e-10, epsrel=1e-10)

    value = integral(np.array(end)) - integral(np.array(start))

    assert value == pytest.approx(expected, rel=1e-7, abs=1e-10)


def test_positive_properties_give_increasing_integrals(property_functions):
    """Verify increasing integrals through the positive broad and narrow features."""
    _, integral = property_functions
    temperature = np.union1d(np.linspace(1.0, 7.0, 301), np.linspace(4.0, 4.001, 101))

    value = integral(temperature)

    assert np.all(np.diff(value) > 0.0)


@pytest.mark.parametrize("temperature, endpoint", [(0.5, 1.0), (8.0, 7.0)])
def test_property_holds_its_validity_endpoint_value(property_functions, temperature, endpoint):
    """Verify that evaluation outside validity holds the nearest endpoint value and warns."""
    property_function, _ = property_functions

    with pytest.warns(RuntimeWarning, match="Synthetic material"):
        value = property_function(np.array(temperature))

    assert value == pytest.approx(property_function(np.array(endpoint)), rel=1e-12)


@pytest.mark.parametrize("temperature, endpoint", [(0.5, 1.0), (8.0, 7.0)])
def test_integral_continues_linearly_with_the_endpoint_property(property_functions, temperature, endpoint):
    """Verify that integral continuation joins the endpoint with the endpoint property's slope."""
    property_function, integral = property_functions
    expected = integral(np.array(endpoint)) + property_function(np.array(endpoint)) * (temperature - endpoint)

    with pytest.warns(RuntimeWarning, match="Synthetic material"):
        value = integral(np.array(temperature))

    assert value == pytest.approx(expected, rel=1e-12, abs=1e-12)


def test_validity_endpoints_are_included_without_warning(property_functions, recwarn):
    """Verify that property and integral evaluation include both validity endpoints."""
    for function in property_functions:
        function(np.array([1.0, 4.0, 7.0]))

    assert not recwarn.list


def test_measured_specific_heat_joins_its_continuation_continuously():
    """Verify that the measured-to-lattice join introduces no jump in specific heat."""
    fit = make_measured_with_lattice_fit(
        measured=[[2.0, 1.0], [4.0, 5.0], [10.0, 2.0]],
        molar_mass=0.1,
        lattice={"debye": {"temperature": 200.0}},
    )

    value = fit(np.array([10.0 - 1e-7, 10.0, 10.0 + 1e-7]))

    np.testing.assert_allclose(value, value[1], rtol=1e-6, atol=0.0)
