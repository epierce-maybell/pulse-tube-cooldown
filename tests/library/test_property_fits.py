"""Tests for material property fits."""

import numpy as np
import pytest

from pulse_tube_cooldown.library.property_fits import (
    PROPERTY_FITS,
    make_fit_from_table,
    make_log_poly_fit,
    make_measured_with_lattice_fit,
    make_rational_poly_fit,
    make_wagner_fit,
)


def test_log_poly_constant_is_independent_of_temperature():
    fit = make_log_poly_fit([2.0])

    # log10(y) = 2 everywhere, so y = 100.
    np.testing.assert_allclose(fit(np.array([4.0, 77.0, 300.0])), [100.0, 100.0, 100.0])


@pytest.mark.parametrize(
    ("coefficients", "temperature", "expected"),
    [
        # log10(y) = 1 + 2 * log10(T): at T = 10, log10(y) = 3; at T = 100, log10(y) = 5.
        ([1.0, 2.0], [10.0, 100.0], [1.0e3, 1.0e5]),
        # log10(y) = log10(T)^2 - 1: at T = 10, log10(y) = 0; at T = 100, log10(y) = 3.
        ([-1.0, 0.0, 1.0], [10.0, 100.0], [1.0, 1.0e3]),
        # log10(y) = 0.5 * log10(T): y = sqrt(T).
        ([0.0, 0.5], [4.0, 25.0], [2.0, 5.0]),
    ],
)
def test_log_poly_nonconstant_matches_hand_calculation(coefficients, temperature, expected):
    fit = make_log_poly_fit(coefficients)

    np.testing.assert_allclose(fit(np.array(temperature)), expected, rtol=1e-12)


def test_log_poly_preserves_temperature_array_shape():
    temperature = np.array([[10.0, 100.0], [1000.0, 10.0]])

    result = make_log_poly_fit([0.0, 1.0])(temperature)

    np.testing.assert_allclose(result, temperature, rtol=1e-12)


def test_log_poly_is_registered():
    assert PROPERTY_FITS["log_poly"] is make_log_poly_fit


# Coefficients are ordered a through i in
# log10(y) = (a + c T^0.5 + e T + g T^1.5 + i T^2) / (1 + b T^0.5 + d T + f T^1.5 + h T^2).
@pytest.mark.parametrize(
    ("coefficients", "temperature", "expected"),
    [
        # Only a: log10(y) = 2 at every temperature.
        ([2.0, 0, 0, 0, 0, 0, 0, 0, 0], [4.0, 300.0], [1.0e2, 1.0e2]),
        # a and c: log10(y) = 1 + sqrt(T), so 3 at T = 4 and 4 at T = 9.
        ([1.0, 0, 1.0, 0, 0, 0, 0, 0, 0], [4.0, 9.0], [1.0e3, 1.0e4]),
        # a and b: log10(y) = 6 / (1 + sqrt(T)), so 2 at T = 4 and 1.5 at T = 9.
        ([6.0, 1.0, 0, 0, 0, 0, 0, 0, 0], [4.0, 9.0], [1.0e2, 10.0**1.5]),
        # d, e, f, g: log10(y) = (T + T^1.5) / (1 + T + T^1.5), so 2 / 3 at T = 1 and 12 / 13 at T = 4.
        ([0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0, 0], [1.0, 4.0], [10.0 ** (2 / 3), 10.0 ** (12 / 13)]),
        # h and i: log10(y) = T^2 / (1 + T^2), so 0.5 at T = 1 and 0.9 at T = 3.
        ([0, 0, 0, 0, 0, 0, 0, 1.0, 1.0], [1.0, 3.0], [10.0**0.5, 10.0**0.9]),
    ],
)
def test_rational_poly_matches_hand_calculation(coefficients, temperature, expected):
    fit = make_rational_poly_fit(coefficients)

    np.testing.assert_allclose(fit(np.array(temperature)), expected, rtol=1e-12)


def test_rational_poly_matches_nist_copper_conductivity():
    # OFHC copper RRR 100 thermal conductivity fit, valid from 4 K to 300 K, from NIST Cryogenic Material Properties:
    # https://trc.nist.gov/cryogenics/materials/OFHC%20Copper/OFHC_Copper_rev1.htm (revised 2010-02-03).
    coefficients = [2.2154, -0.47461, -0.88068, 0.13871, 0.29505, -0.02043, -0.04831, 0.001281, 0.003207]
    temperature = np.array([4.0, 20.0, 77.0, 300.0])

    # NIST publishes no tabulated values for this fit. Expected values in W/(m K) evaluate NIST's equation with the
    # legacy scalar `eval_rational_poly` in branch Advanced_Feature_SIM.
    expected = [642.2969607429585, 2422.510264536503, 547.1996980793672, 396.3239590150121]

    np.testing.assert_allclose(make_rational_poly_fit(coefficients)(temperature), expected, rtol=1e-12)


@pytest.mark.parametrize("count", [8, 10])
def test_rational_poly_rejects_wrong_coefficient_count_when_built(count):
    with pytest.raises(ValueError, match=rf"9 coefficients.*got {count}"):
        make_rational_poly_fit([0.0] * count)


def test_rational_poly_is_registered():
    assert PROPERTY_FITS["rational_poly"] is make_rational_poly_fit


_TABLE = [[10.0, 1.0], [20.0, 3.0], [40.0, 7.0]]


def test_table_returns_values_at_points():
    fit = make_fit_from_table(_TABLE)

    np.testing.assert_allclose(fit(np.array([10.0, 20.0, 40.0])), [1.0, 3.0, 7.0])


def test_table_interpolates_linearly_in_temperature():
    fit = make_fit_from_table(_TABLE)

    # Midpoints of each segment: (1 + 3) / 2 at 15 K and (3 + 7) / 2 at 30 K; a quarter of the way along the second
    # segment, 3 + 0.25 * 4 at 25 K.
    np.testing.assert_allclose(fit(np.array([15.0, 25.0, 30.0])), [2.0, 4.0, 5.0])


def test_table_holds_end_values_outside_points():
    fit = make_fit_from_table(_TABLE)

    np.testing.assert_allclose(fit(np.array([1.0, 9.9, 40.1, 300.0])), [1.0, 1.0, 7.0, 7.0])


def test_table_accepts_array_and_preserves_temperature_shape():
    fit = make_fit_from_table(np.array(_TABLE))

    result = fit(np.array([[10.0, 15.0], [30.0, 40.0]]))

    np.testing.assert_allclose(result, [[1.0, 2.0], [5.0, 7.0]])


def test_table_is_unaffected_by_later_changes_to_its_source():
    source = np.array(_TABLE)
    fit = make_fit_from_table(source)

    source[:, 1] = 0.0

    np.testing.assert_allclose(fit(np.array([15.0])), [2.0])


@pytest.mark.parametrize(
    ("table", "message"),
    [
        ([], "at least 2 points"),
        ([[10.0, 1.0]], "at least 2 points"),
        ([[10.0, 1.0], [20.0]], "2 values"),
        ([[10.0, 1.0, 0.0], [20.0, 3.0, 0.0]], "2 values"),
        ([[20.0, 1.0], [10.0, 3.0]], "strictly increasing"),
        ([[10.0, 1.0], [10.0, 3.0]], "strictly increasing"),
        ([[10.0, 1.0], [20.0, np.nan]], "finite"),
        ([[10.0, 1.0], [np.inf, 3.0]], "finite"),
    ],
)
def test_table_rejects_invalid_tables_when_built(table, message):
    with pytest.raises(ValueError, match=message):
        make_fit_from_table(table)


def test_table_is_registered():
    assert PROPERTY_FITS["table"] is make_fit_from_table


def _wagner_example():
    # ln(p / 2e5) = (5 / T) * (-4 * theta + 2 * theta^1.5), theta = 1 - T / 5.
    return make_wagner_fit(
        critical_temperature=5.0, critical_pressure=2.0e5, coefficients=[-4.0, 2.0], exponents=[1.0, 1.5]
    )


def test_wagner_matches_hand_calculation():
    fit = _wagner_example()

    # At 2.5 K, theta = 0.5 and ln(p / pc) = 2 * (-2 + 2 * 0.5^1.5) = -4 + 2^0.5.
    # At 4 K, theta = 0.2 and ln(p / pc) = 1.25 * (-0.8 + 2 * 0.2^1.5).
    expected = [2.0e5 * np.exp(-4.0 + np.sqrt(2.0)), 2.0e5 * np.exp(1.25 * (-0.8 + 2.0 * 0.2**1.5))]

    np.testing.assert_allclose(fit(np.array([2.5, 4.0])), expected, rtol=1e-12)


def test_wagner_reaches_critical_pressure_at_critical_temperature():
    np.testing.assert_allclose(_wagner_example()(np.array([5.0])), [2.0e5], rtol=1e-12)


def test_wagner_is_infinite_above_critical_temperature():
    fit = _wagner_example()

    assert np.all(np.isposinf(fit(np.array([5.0001, 10.0, 300.0]))))


def test_wagner_preserves_temperature_array_shape():
    assert _wagner_example()(np.array([[2.5, 4.0], [5.0, 6.0]])).shape == (2, 2)


def test_wagner_rejects_mismatched_coefficients_and_exponents_when_built():
    with pytest.raises(ValueError, match="same number"):
        make_wagner_fit(critical_temperature=5.0, critical_pressure=2.0e5, coefficients=[-4.0, 2.0], exponents=[1.0])


def test_wagner_is_registered():
    assert PROPERTY_FITS["wagner"] is make_wagner_fit


_GAS_CONSTANT = 8.314462618

# Zeros ending at 0.02 K, leaving the lattice terms above.
_NO_MEASUREMENT = [[0.01, 0.0], [0.02, 0.0]]


def test_measured_with_lattice_follows_measured_table_to_its_last_point():
    fit = make_measured_with_lattice_fit(
        measured=[[2.0, 1.0], [4.0, 5.0], [10.0, 2.0]], molar_mass=0.1, lattice={"debye": {"temperature": 200.0}}
    )

    np.testing.assert_allclose(fit(np.array([2.0, 3.0, 4.0, 7.0, 10.0])), [1.0, 3.0, 5.0, 3.5, 2.0])


def test_measured_with_lattice_debye_term_follows_t_cubed_law_far_below_debye_temperature():
    fit = make_measured_with_lattice_fit(
        measured=_NO_MEASUREMENT, molar_mass=0.5, lattice={"debye": {"temperature": 1000.0, "atoms": 4}}
    )

    # C = (12 pi^4 / 5) n R (T / theta)^3 per mole for n atoms, here per kg.
    temperature = np.array([5.0, 20.0])
    expected = 12.0 * np.pi**4 / 5.0 * 4 * _GAS_CONSTANT * (temperature / 1000.0) ** 3 / 0.5

    np.testing.assert_allclose(fit(temperature), expected, rtol=1e-6)


def test_measured_with_lattice_debye_term_approaches_dulong_petit_far_above_debye_temperature():
    fit = make_measured_with_lattice_fit(
        measured=_NO_MEASUREMENT, molar_mass=0.5, lattice={"debye": {"temperature": 10.0, "atoms": 4}}
    )

    # C = 3 n R (1 - (theta / T)^2 / 20) per mole as T / theta grows.
    expected = 3 * 4 * _GAS_CONSTANT * (1.0 - (10.0 / 1000.0) ** 2 / 20.0) / 0.5

    np.testing.assert_allclose(fit(np.array([1000.0])), [expected], rtol=1e-6)


def test_measured_with_lattice_einstein_branch_matches_hand_calculation():
    fit = make_measured_with_lattice_fit(
        measured=_NO_MEASUREMENT,
        molar_mass=1.0,
        lattice={"debye": {"temperature": 1.0e6}, "einstein": [{"temperature": 150.0, "branches": 2}]},
    )

    # At T = theta_E, x = 1 and each branch carries R x^2 e^x / (e^x - 1)^2 = R e / (e - 1)^2.
    expected = 2 * _GAS_CONSTANT * np.e / (np.e - 1.0) ** 2

    np.testing.assert_allclose(fit(np.array([150.0])), [expected], rtol=1e-9)


def test_measured_with_lattice_adds_electronic_term():
    fit = make_measured_with_lattice_fit(
        measured=_NO_MEASUREMENT, molar_mass=0.25, lattice={"debye": {"temperature": 1.0e6}, "electronic": 0.015}
    )

    # gamma T per mole, here per kg.
    np.testing.assert_allclose(fit(np.array([100.0])), [0.015 * 100.0 / 0.25], rtol=1e-6)


def test_measured_with_lattice_anharmonic_coefficient_divides_by_one_minus_alpha_t():
    plain = make_measured_with_lattice_fit(
        measured=_NO_MEASUREMENT, molar_mass=1.0, lattice={"debye": {"temperature": 10.0}}
    )
    anharmonic = make_measured_with_lattice_fit(
        measured=_NO_MEASUREMENT, molar_mass=1.0, lattice={"debye": {"temperature": 10.0, "anharmonic": 2.0e-4}}
    )

    temperature = np.array([300.0])

    np.testing.assert_allclose(anharmonic(temperature), plain(temperature) / (1.0 - 2.0e-4 * 300.0), rtol=1e-12)


def test_measured_with_lattice_magnetic_remainder_falls_as_inverse_square_above_last_point():
    fit = make_measured_with_lattice_fit(
        measured=[[5.0, 40.0], [10.0, 8.0]], molar_mass=1.0, lattice={"debye": {"temperature": 1.0e6}}
    )

    # With no lattice heat capacity, 8 J/(kg K) at 10 K falls to 8 (10 / T)^2 above.
    np.testing.assert_allclose(fit(np.array([10.0, 20.0, 40.0])), [8.0, 2.0, 0.5], rtol=1e-9)


def test_measured_with_lattice_is_continuous_at_last_measured_point():
    fit = make_measured_with_lattice_fit(
        measured=[[5.0, 40.0], [20.0, 80.0]],
        molar_mass=0.3,
        lattice={"debye": {"temperature": 150.0, "atoms": 3}, "einstein": [{"temperature": 250.0, "branches": 2}]},
    )

    below, above = fit(np.array([20.0 - 1e-9, 20.0 + 1e-9]))

    assert above == pytest.approx(below, rel=1e-6)


def test_measured_with_lattice_preserves_temperature_array_shape():
    fit = make_measured_with_lattice_fit(
        measured=[[5.0, 40.0], [20.0, 80.0]], molar_mass=0.3, lattice={"debye": {"temperature": 150.0}}
    )

    assert fit(np.array([[5.0, 10.0], [100.0, 300.0]])).shape == (2, 2)


def test_measured_with_lattice_is_registered():
    assert PROPERTY_FITS["measured_with_lattice"] is make_measured_with_lattice_fit
