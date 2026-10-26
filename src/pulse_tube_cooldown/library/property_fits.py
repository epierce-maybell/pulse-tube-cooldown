"""Fit forms for temperature-dependent material properties.

Each fit evaluates a property elementwise over temperatures in K from the coefficients given in a material file.
Fits apply no validity range; the material loader is responsible for range handling.
"""

from collections.abc import Sequence

import numpy as np
from numpy.polynomial import polynomial
from numpy.typing import NDArray

__all__ = ["PROPERTY_FITS", "eval_log_poly"]


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


PROPERTY_FITS = {
    "log_poly": eval_log_poly,
}
