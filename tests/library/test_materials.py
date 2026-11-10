"""Tests for loading materials from TOML files."""

import textwrap
import warnings

import numpy as np
import pytest
from scipy.integrate import quad

from pulse_tube_cooldown.core.errors import ConfigurationError
from pulse_tube_cooldown.core.types import Material
from pulse_tube_cooldown.library.materials.load import load_material

# Specific heat is log10(y) = log10(T), so y = T. Thermal conductivity is log10(y) = 1 + sqrt(T), so 1e3 at T = 4 and
# 1e4 at T = 9.
_POLY_MATERIAL = """
    name = "Poly Metal"
    density = 8000.0

    [specific_heat]
    form = "log_poly"
    coefficients = [0.0, 1.0]
    range = [4.0, 300.0]
    source = "Specific heat reference"

    [thermal_conductivity]
    form = "rational_poly"
    coefficients = [1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    range = [4.0, 9.0]
    source = "Conductivity reference"
"""

_TABLE_MATERIAL = """
    name = "Table Metal"
    density = 2700.0

    [specific_heat]
    form = "table"
    table = [[10.0, 1.0], [20.0, 3.0], [40.0, 7.0]]
    range = [10.0, 40.0]
    source = "Specific heat table"

    [thermal_conductivity]
    form = "table"
    table = [[10.0, 100.0], [40.0, 400.0]]
    range = [10.0, 40.0]
    source = "Conductivity table"
"""


def _write(directory, filename, text):
    (directory / filename).write_text(textwrap.dedent(text), encoding="utf-8")


def test_loads_analytic_fits_by_name(tmp_path):
    _write(tmp_path, "poly.toml", _POLY_MATERIAL)

    material = load_material("Poly Metal", directory=tmp_path)

    assert isinstance(material, Material)
    assert material.name == "Poly Metal"
    np.testing.assert_allclose(material.specific_heat(np.array([10.0, 100.0])), [10.0, 100.0], rtol=1e-12)
    np.testing.assert_allclose(material.thermal_conductivity(np.array([4.0, 9.0])), [1.0e3, 1.0e4], rtol=1e-12)


def test_loads_table_fits_by_name(tmp_path):
    _write(tmp_path, "table.toml", _TABLE_MATERIAL)

    material = load_material("Table Metal", directory=tmp_path)

    np.testing.assert_allclose(material.specific_heat(np.array([15.0, 30.0])), [2.0, 5.0])
    np.testing.assert_allclose(material.thermal_conductivity(np.array([25.0])), [250.0])


def test_loads_density(tmp_path):
    _write(tmp_path, "poly.toml", _POLY_MATERIAL)
    _write(tmp_path, "table.toml", _TABLE_MATERIAL)

    assert load_material("Poly Metal", directory=tmp_path).density == 8000.0
    assert load_material("Table Metal", directory=tmp_path).density == 2700.0


def test_loads_material_without_thermal_conductivity(tmp_path):
    text = _POLY_MATERIAL.split("[thermal_conductivity]")[0]
    _write(tmp_path, "poly.toml", text)

    material = load_material("Poly Metal", directory=tmp_path)

    np.testing.assert_allclose(material.specific_heat(np.array([10.0, 100.0])), [10.0, 100.0], rtol=1e-12)
    assert material.thermal_conductivity is None
    assert material.conductivity_integral is None


def test_evaluates_at_nearest_range_end_outside_range(tmp_path):
    _write(tmp_path, "poly.toml", _POLY_MATERIAL)

    material = load_material("Poly Metal", directory=tmp_path)

    # Specific heat equals T inside [4, 300] K; conductivity is 1e3 at 4 K and 1e4 at 9 K.
    with pytest.warns(RuntimeWarning):
        specific_heat = material.specific_heat(np.array([1.0, 4.0, 300.0, 1000.0]))
        thermal_conductivity = material.thermal_conductivity(np.array([2.0, 16.0]))

    np.testing.assert_allclose(specific_heat, [4.0, 4.0, 300.0, 300.0])
    np.testing.assert_allclose(thermal_conductivity, [1.0e3, 1.0e4], rtol=1e-12)


@pytest.mark.parametrize(
    ("temperature", "side"),
    [([1.0, 10.0], "below"), ([10.0, 1000.0], "above")],
)
def test_warns_naming_material_property_and_side_outside_range(tmp_path, temperature, side):
    _write(tmp_path, "poly.toml", _POLY_MATERIAL)

    material = load_material("Poly Metal", directory=tmp_path)

    with pytest.warns(RuntimeWarning, match=rf"Poly Metal.*specific_heat.*{side}.*4.*300"):
        material.specific_heat(np.array(temperature))


def test_does_not_warn_within_range_including_its_ends(tmp_path, recwarn):
    _write(tmp_path, "poly.toml", _POLY_MATERIAL)

    material = load_material("Poly Metal", directory=tmp_path)
    material.specific_heat(np.array([4.0, 100.0, 300.0]))
    material.thermal_conductivity(np.array([4.0, 9.0]))
    material.conductivity_integral(np.array([4.0, 9.0]))

    assert not recwarn.list


def _poly_conductivity_antiderivative(temperature):
    # k = 10^(1 + sqrt(T)). With u = sqrt(T) and a = ln(10), the integral of 20 u e^(a u) du is
    # 20 e^(a u) (u / a - 1 / a^2).
    u = np.sqrt(temperature)
    a = np.log(10.0)
    return 20.0 * np.exp(a * u) * (u / a - 1.0 / a**2)


def test_conductivity_integral_of_constant_conductivity_grows_linearly_from_range_start(tmp_path):
    text = _POLY_MATERIAL.replace(
        "coefficients = [1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]",
        "coefficients = [2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]",
    )
    _write(tmp_path, "poly.toml", text)

    material = load_material("Poly Metal", directory=tmp_path)

    # k = 10^2 = 100 W/(m K) over [4, 9] K.
    np.testing.assert_allclose(material.conductivity_integral(np.array([4.0, 6.5, 9.0])), [0.0, 250.0, 500.0])


def test_conductivity_integral_of_analytic_conductivity_matches_antiderivative(tmp_path):
    _write(tmp_path, "poly.toml", _POLY_MATERIAL)

    material = load_material("Poly Metal", directory=tmp_path)

    temperature = np.array([4.0, 5.0, 6.25, 8.0, 9.0])
    expected = _poly_conductivity_antiderivative(temperature) - _poly_conductivity_antiderivative(4.0)
    np.testing.assert_allclose(material.conductivity_integral(temperature), expected, rtol=1e-8, atol=1e-8)


def test_conductivity_integral_extends_linearly_outside_range(tmp_path):
    _write(tmp_path, "poly.toml", _POLY_MATERIAL)

    material = load_material("Poly Metal", directory=tmp_path)

    # k is 1e3 at 4 K and 1e4 at 9 K. Below, -1e3 * (4 - 2); above, the integral to 9 K plus 1e4 * (16 - 9).
    at_upper = _poly_conductivity_antiderivative(9.0) - _poly_conductivity_antiderivative(4.0)
    with pytest.warns(RuntimeWarning):
        integral = material.conductivity_integral(np.array([2.0, 16.0]))

    np.testing.assert_allclose(integral, [-2.0e3, at_upper + 7.0e4], rtol=1e-8)


@pytest.mark.parametrize(
    ("temperature", "side"),
    [([1.0, 5.0], "below"), ([5.0, 100.0], "above")],
)
def test_conductivity_integral_warns_naming_material_and_side_outside_range(tmp_path, temperature, side):
    _write(tmp_path, "poly.toml", _POLY_MATERIAL)

    material = load_material("Poly Metal", directory=tmp_path)

    with pytest.warns(RuntimeWarning, match=rf"Poly Metal.*conductivity_integral.*{side}.*4.*9"):
        material.conductivity_integral(np.array(temperature))


def test_conductivity_integral_warns_once_per_side_under_default_filter(tmp_path):
    _write(tmp_path, "poly.toml", _POLY_MATERIAL)
    material = load_material("Poly Metal", directory=tmp_path)

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("default")
        for temperature in [1.0, 2.0, 3.0, 40.0, 50.0]:
            material.conductivity_integral(np.array([temperature]))

    assert len(caught) == 2


def test_repeated_evaluation_outside_range_warns_once_per_side_under_default_filter(tmp_path):
    _write(tmp_path, "poly.toml", _POLY_MATERIAL)
    material = load_material("Poly Metal", directory=tmp_path)

    # A solver evaluates the same property from one line at many temperatures.
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("default")
        for temperature in [1.0, 2.0, 3.0, 400.0, 500.0]:
            material.specific_heat(np.array([temperature]))

    assert len(caught) == 2


def test_accepts_directory_as_string(tmp_path):
    _write(tmp_path, "poly.toml", _POLY_MATERIAL)

    assert load_material("Poly Metal", directory=str(tmp_path)).name == "Poly Metal"


def test_accepts_integer_range_bounds(tmp_path):
    _write(tmp_path, "poly.toml", _POLY_MATERIAL.replace("range = [4.0, 300.0]", "range = [4, 300]"))

    material = load_material("Poly Metal", directory=tmp_path)

    with pytest.warns(RuntimeWarning):
        np.testing.assert_allclose(material.specific_heat(np.array([1.0, 1000.0])), [4.0, 300.0])


def test_preserves_temperature_array_shape(tmp_path):
    _write(tmp_path, "poly.toml", _POLY_MATERIAL)

    material = load_material("Poly Metal", directory=tmp_path)

    np.testing.assert_allclose(
        material.specific_heat(np.array([[4.0, 10.0], [100.0, 300.0]])), [[4.0, 10.0], [100.0, 300.0]]
    )


def test_finds_material_by_name_field_rather_than_filename(tmp_path):
    _write(tmp_path, "a.toml", _TABLE_MATERIAL)
    _write(tmp_path, "b.toml", _POLY_MATERIAL)

    assert load_material("Poly Metal", directory=tmp_path).name == "Poly Metal"
    assert load_material("Table Metal", directory=tmp_path).name == "Table Metal"


def test_unknown_name_is_reported(tmp_path):
    _write(tmp_path, "poly.toml", _POLY_MATERIAL)

    with pytest.raises(ConfigurationError, match=r"0 materials named 'Copper'"):
        load_material("Copper", directory=tmp_path)


def test_duplicate_names_are_reported(tmp_path):
    _write(tmp_path, "first.toml", _POLY_MATERIAL)
    _write(tmp_path, "second.toml", _POLY_MATERIAL)

    with pytest.raises(ConfigurationError, match=r"2 materials named 'Poly Metal'"):
        load_material("Poly Metal", directory=tmp_path)


def test_loads_material_without_sources(tmp_path):
    _write(tmp_path, "poly.toml", _POLY_MATERIAL.replace('source = "Specific heat reference"\n', ""))

    np.testing.assert_allclose(load_material("Poly Metal", directory=tmp_path).specific_heat(np.array([10.0])), [10.0])


def test_conductivity_integral_has_conductivity_as_its_slope(tmp_path):
    _write(tmp_path, "poly.toml", _POLY_MATERIAL)

    material = load_material("Poly Metal", directory=tmp_path)

    temperature = np.array([4.5, 5.0, 7.3, 8.5])
    step = 1e-4
    slope = (
        material.conductivity_integral(temperature + step) - material.conductivity_integral(temperature - step)
    ) / (2.0 * step)
    np.testing.assert_allclose(slope, material.thermal_conductivity(temperature), rtol=1e-6)


_RATIONAL_POLY_SPECIFIC_HEAT_MATERIAL = """
    name = "Rational Metal"
    density = 8000.0

    [specific_heat]
    form = "rational_poly"
    coefficients = [0.0, 0.5, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    range = [4.0, 300.0]
"""


def _quadrature(material, temperature, lower, breakpoints=()):
    def specific_heat(t):
        return material.specific_heat(np.array([t])).item()

    return [
        quad(
            specific_heat,
            lower,
            t,
            points=[p for p in breakpoints if lower < p < t] or None,
            limit=1000,
            epsabs=1e-12,
            epsrel=1e-12,
        )[0]
        for t in temperature
    ]


@pytest.mark.parametrize(
    ("filename", "text", "name", "temperature", "breakpoints"),
    [
        ("poly.toml", _POLY_MATERIAL, "Poly Metal", [4.0, 7.0, 55.0, 300.0], ()),
        (
            "rational.toml",
            _RATIONAL_POLY_SPECIFIC_HEAT_MATERIAL,
            "Rational Metal",
            [4.0, 9.0, 120.0, 300.0],
            (),
        ),
        ("table.toml", _TABLE_MATERIAL, "Table Metal", [10.0, 15.0, 25.0, 40.0], (20.0,)),
    ],
    ids=["log_poly", "rational_poly", "table"],
)
def test_specific_heat_integral_agrees_with_quadrature(tmp_path, filename, text, name, temperature, breakpoints):
    _write(tmp_path, filename, text)

    material = load_material(name, directory=tmp_path)

    np.testing.assert_allclose(
        material.specific_heat_integral(np.array(temperature)),
        _quadrature(material, temperature, temperature[0], breakpoints),
        rtol=1e-7,
        atol=1e-9,
    )


def test_specific_heat_integral_of_measured_with_lattice_agrees_with_quadrature_across_its_peak():
    material = load_material("Gd2O2S")

    # Range [2, 300] K; the measured table, on a 0.05 K grid, peaks at 5.24 K and joins the lattice continuation at
    # 10 K. Between knots on the peak, the integral's spline is within 1e-5 of the quadrature.
    temperature = [2.0, 4.0, 5.24, 6.0, 10.0, 50.0, 300.0]
    measured_grid = np.round(np.arange(2.05, 10.01, 0.05), 2)
    np.testing.assert_allclose(
        material.specific_heat_integral(np.array(temperature)),
        _quadrature(material, temperature, 2.0, breakpoints=measured_grid),
        rtol=1e-5,
    )


def test_specific_heat_integral_of_table_matches_piecewise_integral(tmp_path):
    _write(tmp_path, "table.toml", _TABLE_MATERIAL)

    material = load_material("Table Metal", directory=tmp_path)

    # Points (10, 1), (20, 3), (40, 7) are collinear in slope 0.2, so PCHIP is the line c = 1 + 0.2 (T - 10),
    # integrating to (T - 10) + 0.1 (T - 10)^2.
    np.testing.assert_allclose(
        material.specific_heat_integral(np.array([10.0, 15.0, 25.0, 40.0])), [0.0, 7.5, 37.5, 120.0], rtol=1e-8
    )


def test_specific_heat_integral_extends_linearly_outside_range(tmp_path):
    _write(tmp_path, "poly.toml", _POLY_MATERIAL)

    material = load_material("Poly Metal", directory=tmp_path)

    # c = T on [4, 300] K. Below, -4 * (4 - 2); above, (300^2 - 4^2) / 2 plus 300 * (400 - 300).
    with pytest.warns(RuntimeWarning):
        integral = material.specific_heat_integral(np.array([2.0, 400.0]))

    np.testing.assert_allclose(integral, [-8.0, 44992.0 + 30000.0], rtol=1e-8)


@pytest.mark.parametrize(
    ("temperature", "side"),
    [([1.0, 5.0], "below"), ([5.0, 1000.0], "above")],
)
def test_specific_heat_integral_warns_naming_material_and_side_outside_range(tmp_path, temperature, side):
    _write(tmp_path, "poly.toml", _POLY_MATERIAL)

    material = load_material("Poly Metal", directory=tmp_path)

    with pytest.warns(RuntimeWarning, match=rf"Poly Metal.*specific_heat_integral.*{side}.*4.*300"):
        material.specific_heat_integral(np.array(temperature))


def test_conductivity_integral_of_table_extending_past_range_integrates_within_range(tmp_path):
    text = _TABLE_MATERIAL.replace(
        "table = [[10.0, 100.0], [40.0, 400.0]]", "table = [[5.0, 100.0], [25.0, 400.0], [50.0, 775.0]]"
    )
    _write(tmp_path, "table.toml", text)

    material = load_material("Table Metal", directory=tmp_path)

    # Over the range [10, 40] K, k = 100 + 15 (T - 5), integrating to 100 (T - 5) + 7.5 (T - 5)^2.
    # To 25 K: (2000 + 3000) - (500 + 187.5) = 4312.5. To 40 K: (3500 + 9187.5) - 687.5 = 12000.
    np.testing.assert_allclose(
        material.conductivity_integral(np.array([10.0, 25.0, 40.0])), [0.0, 4312.5, 12000.0], rtol=1e-8
    )
