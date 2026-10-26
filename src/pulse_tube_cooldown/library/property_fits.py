"""Fit forms for temperature-dependent material properties.

Each fit evaluates a property elementwise over temperatures in K from the coefficients given in a material file.
Fits apply no validity range; the material loader is responsible for range handling.
"""

from collections.abc import Sequence

import numpy as np
from numpy.polynomial import polynomial
from numpy.typing import NDArray

__all__ = ["PROPERTY_FITS", "eval_log_poly", "eval_rational_poly"]


def eval_log_poly(temperature: NDArray[np.float64], coefficients: Sequence[float]) -> NDArray[np.float64]:
    """Evaluates a polynomial in log-log space, the standard NIST cryogenic property form.

    The property y satisfies log10(y) = sum(c_i * log10(T)^i), summed over the coefficients c_i in ascending
    order of power.

    Args:
        temperature: Positive temperatures in K, of any shape.
        coefficients: Polynomial coefficients, starting with the constant term.

    Returns:
        Property values in the units of the fit, with the same shape as `temperature`.
    """
    return 10.0 ** polynomial.polyval(np.log10(temperature), coefficients)


def eval_rational_poly(temperature: NDArray[np.float64], coefficients: Sequence[float]) -> NDArray[np.float64]:
    """Evaluates the NIST rational polynomial in sqrt(T), used for OFHC copper thermal conductivity.

    The property y satisfies log10(y) = (a + c T^0.5 + e T + g T^1.5 + i T^2) / (1 + b T^0.5 + d T + f T^1.5 + h T^2).

    Args:
        temperature: Positive temperatures in K, of any shape.
        coefficients: The nine coefficients a through i, in alphabetical order.

    Returns:
        Property values in the units of the fit, with the same shape as `temperature`.

    Raises:
        ValueError: If `coefficients` does not contain exactly nine values.
    """
    if len(coefficients) != 9:
        raise ValueError(f"rational_poly takes 9 coefficients, a through i; got {len(coefficients)}.")
    a, b, c, d, e, f, g, h, i = coefficients
    root = np.sqrt(temperature)
    numerator = polynomial.polyval(root, [a, c, e, g, i])
    denominator = polynomial.polyval(root, [1.0, b, d, f, h])
    return 10.0 ** (numerator / denominator)


PROPERTY_FITS = {
    "log_poly": eval_log_poly,
    "rational_poly": eval_rational_poly,
}
