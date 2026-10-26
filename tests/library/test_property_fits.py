"""Tests for material property fits."""

import numpy as np
import pytest

from pulse_tube_cooldown.library.property_fits import PROPERTY_FITS, eval_log_poly


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
