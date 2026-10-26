"""Tests for material property fits."""

import numpy as np
import pytest

from pulse_tube_cooldown.library.property_fits import PROPERTY_FITS, eval_log_poly, eval_rational_poly


def test_log_poly_constant_is_independent_of_temperature():
    temperature = np.array([4.0, 77.0, 300.0])

    # log10(y) = 2 everywhere, so y = 100.
    np.testing.assert_allclose(eval_log_poly(temperature, [2.0]), [100.0, 100.0, 100.0])


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
    np.testing.assert_allclose(eval_log_poly(np.array(temperature), coefficients), expected, rtol=1e-12)


def test_log_poly_preserves_temperature_array_shape():
    temperature = np.array([[10.0, 100.0], [1000.0, 10.0]])

    result = eval_log_poly(temperature, [0.0, 1.0])

    np.testing.assert_allclose(result, temperature, rtol=1e-12)


def test_log_poly_is_registered():
    assert PROPERTY_FITS["log_poly"] is eval_log_poly


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
    np.testing.assert_allclose(eval_rational_poly(np.array(temperature), coefficients), expected, rtol=1e-12)


def test_rational_poly_matches_nist_copper_conductivity():
    # OFHC copper RRR 100 thermal conductivity fit, valid from 4 K to 300 K, from NIST Cryogenic Material Properties:
    # https://trc.nist.gov/cryogenics/materials/OFHC%20Copper/OFHC_Copper_rev1.htm (revised 2010-02-03).
    coefficients = [2.2154, -0.47461, -0.88068, 0.13871, 0.29505, -0.02043, -0.04831, 0.001281, 0.003207]
    temperature = np.array([4.0, 20.0, 77.0, 300.0])

    # NIST publishes no tabulated values for this fit. Expected values in W/(m K) evaluate NIST's equation with the
    # legacy scalar `eval_rational_poly` in branch Advanced_Feature_SIM.
    expected = [642.2969607429585, 2422.510264536503, 547.1996980793672, 396.3239590150121]

    np.testing.assert_allclose(eval_rational_poly(temperature, coefficients), expected, rtol=1e-12)


@pytest.mark.parametrize("count", [8, 10])
def test_rational_poly_rejects_wrong_coefficient_count(count):
    with pytest.raises(ValueError, match=rf"9 coefficients.*got {count}"):
        eval_rational_poly(np.array([10.0]), [0.0] * count)


def test_rational_poly_is_registered():
    assert PROPERTY_FITS["rational_poly"] is eval_rational_poly
