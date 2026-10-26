"""Fit forms for temperature-dependent material properties.

Each fit form builds a property function from the coefficients given in a material file. Building checks the
coefficients once; the returned function maps temperatures in K to property values elementwise and does no checking.

Validity ranges belong to the material loader. Analytic fits evaluate at any temperature; tables hold their end
values beyond their first and last points.
"""

from collections.abc import Callable, Sequence

import numpy as np
from numpy.polynomial import polynomial
from numpy.typing import ArrayLike, NDArray

__all__ = ["PROPERTY_FITS", "PropertyFit", "make_fit_from_table", "make_log_poly_fit", "make_rational_poly_fit"]

type PropertyFit = Callable[[NDArray[np.float64]], NDArray[np.float64]]


def make_log_poly_fit(coefficients: Sequence[float]) -> PropertyFit:
    """Builds a polynomial in log-log space, the standard NIST cryogenic property form.

    The property y satisfies log10(y) = sum(c_i * log10(T)^i), summed over the coefficients c_i in ascending
    order of power.

    Args:
        coefficients: Polynomial coefficients, starting with the constant term.

    Returns:
        A function from positive temperatures in K, of any shape, to property values in the units of the fit.
    """
    c = np.array(coefficients, dtype=np.float64)

    def fit(temperature: NDArray[np.float64]) -> NDArray[np.float64]:
        return 10.0 ** polynomial.polyval(np.log10(temperature), c)

    return fit


def make_rational_poly_fit(coefficients: Sequence[float]) -> PropertyFit:
    """Builds the NIST rational polynomial in sqrt(T), used for OFHC copper thermal conductivity.

    The property y satisfies log10(y) = (a + c T^0.5 + e T + g T^1.5 + i T^2) / (1 + b T^0.5 + d T + f T^1.5 + h T^2).

    Args:
        coefficients: The nine coefficients a through i, in alphabetical order.

    Returns:
        A function from positive temperatures in K, of any shape, to property values in the units of the fit.

    Raises:
        ValueError: If `coefficients` does not contain exactly nine values.
    """
    if len(coefficients) != 9:
        raise ValueError(f"rational_poly takes 9 coefficients, a through i; got {len(coefficients)}.")

    a, b, c, d, e, f, g, h, i = coefficients
    numerator = np.array([a, c, e, g, i], dtype=np.float64)
    denominator = np.array([1.0, b, d, f, h], dtype=np.float64)

    def fit(temperature: NDArray[np.float64]) -> NDArray[np.float64]:
        root = np.sqrt(temperature)
        return 10.0 ** (polynomial.polyval(root, numerator) / polynomial.polyval(root, denominator))

    return fit


def make_fit_from_table(table: ArrayLike) -> PropertyFit:
    """Builds linear interpolation in temperature over tabulated property values.

    Temperatures below the first point or above the last take that point's value.

    Args:
        table: Points as rows of [temperature in K, property value], with strictly increasing temperatures.

    Returns:
        A function from temperatures in K, of any shape, to property values in the units of the table.

    Raises:
        ValueError: If the table has fewer than two points, a row without exactly two values, temperatures that are
            not strictly increasing, or non-finite entries.
    """
    point_temperature, point_value = _parse_table(table)

    def fit(temperature: NDArray[np.float64]) -> NDArray[np.float64]:
        return np.interp(temperature, point_temperature, point_value)

    return fit


def _parse_table(table: ArrayLike) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Converts table rows into point temperatures and values, raising ValueError for an invalid table."""
    if len(table) < 2:
        raise ValueError(f"A property table needs at least 2 points; got {len(table)}.")

    try:
        points = np.array(table, dtype=np.float64)
    except ValueError as error:
        raise ValueError("Each property table row must hold 2 values, a temperature and a property value.") from error

    if points.ndim != 2 or points.shape[1] != 2:
        raise ValueError(f"Each property table row must hold 2 values; got rows of shape {points.shape[1:]}.")

    if not np.all(np.isfinite(points)):
        raise ValueError("Property table entries must be finite.")

    point_temperature, point_value = points.T

    if np.any(np.diff(point_temperature) <= 0):
        raise ValueError(f"Property table temperatures must be strictly increasing; got {point_temperature.tolist()}.")

    return point_temperature, point_value


PROPERTY_FITS = {
    "log_poly": make_log_poly_fit,
    "rational_poly": make_rational_poly_fit,
    "table": make_fit_from_table,
}
