"""Fit forms and validity ranges for temperature-dependent material and gas properties.

Each fit form builds a property function from the coefficients given in a material or gas file. Building checks the
coefficients once; the returned function maps temperatures in K to property values elementwise and does no checking.
Analytic fits evaluate at any temperature; tables hold their end values beyond their first and last points.

A property's table in a file holds its fit's `form`, its validity `range` in K, an optional `source` citing the data,
and the form's arguments. An infinite end of the range leaves the fit unbounded on that side.
"""

import warnings
from collections.abc import Callable, Mapping, Sequence
from typing import Any

import numpy as np
from numpy.polynomial import polynomial
from numpy.typing import ArrayLike, NDArray
from scipy.constants import gas_constant

__all__ = [
    "PROPERTY_FITS",
    "PropertyFit",
    "make_fit",
    "make_fit_from_table",
    "make_log_poly_fit",
    "make_measured_with_lattice_fit",
    "make_rational_poly_fit",
    "make_wagner_fit",
    "warn_above",
    "warn_below",
]

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


def make_measured_with_lattice_fit(measured: ArrayLike, molar_mass: float, lattice: Mapping[str, Any]) -> PropertyFit:
    """Builds a measured specific heat continued above its last point by lattice, electronic, and magnetic terms.

    At and below the measured table's last temperature T_m, the specific heat interpolates the table linearly. Above,
    it is the molar lattice and electronic heat capacity divided by the molar mass, plus the table's excess over those
    terms at T_m, falling as (T_m / T)^2.

    The lattice table holds, with x = theta / T and R the molar gas constant:

    - `debye`: `temperature` theta, `atoms` n (default 1), and `anharmonic` alpha (default 0), giving three acoustic
      branches per atom, 9 n R x^-3 integral_0^x u^4 e^u / (e^u - 1)^2 du / (1 - alpha T).
    - `einstein`: a list of optical modes, each with `temperature` theta, `branches` b, and `anharmonic` alpha
      (default 0), giving b R x^2 e^x / (e^x - 1)^2 / (1 - alpha T). Defaults to none.
    - `electronic`: the Sommerfeld coefficient gamma in J/(mol K^2), giving gamma T. Defaults to 0.

    Args:
        measured: Points as rows of [temperature in K, specific heat in J/(kg K)], with strictly increasing
            temperatures.
        molar_mass: Mass of one formula unit's mole in kg/mol.
        lattice: The terms above.

    Returns:
        A function from positive temperatures in K, of any shape, to specific heats in J/(kg K).

    Raises:
        ValueError: If the measured table is invalid, as for `make_fit_from_table`.
    """
    point_temperature, point_value = _parse_table(measured)
    debye = lattice["debye"]
    modes = [
        (mode["temperature"], mode["branches"], mode.get("anharmonic", 0.0)) for mode in lattice.get("einstein", [])
    ]
    electronic = lattice.get("electronic", 0.0)

    # Gauss-Legendre rule on [-1, 1] for the Debye integral, truncated at u = 50.
    nodes, weights = np.polynomial.legendre.leggauss(64)

    def acoustic(temperature: NDArray[np.float64]) -> NDArray[np.float64]:
        x = debye["temperature"] / temperature
        upper = np.minimum(x, 50.0)[..., np.newaxis]
        u = upper * (nodes + 1.0) / 2.0
        integral = upper[..., 0] / 2.0 * np.sum(weights * u**4 * np.exp(-u) / np.expm1(-u) ** 2, axis=-1)
        return 9.0 * gas_constant * integral / x**3

    def optical(temperature: NDArray[np.float64], theta: float) -> NDArray[np.float64]:
        x = theta / temperature
        return gas_constant * x**2 * np.exp(-x) / np.expm1(-x) ** 2

    def molar_heat_capacity(temperature: NDArray[np.float64]) -> NDArray[np.float64]:
        total = debye.get("atoms", 1) * acoustic(temperature) / (1.0 - debye.get("anharmonic", 0.0) * temperature)
        total = total + electronic * temperature
        for theta, branches, anharmonic in modes:
            total = total + branches * optical(temperature, theta) / (1.0 - anharmonic * temperature)
        return total

    last = point_temperature[-1]
    remainder = point_value[-1] - molar_heat_capacity(np.array([last]))[0] / molar_mass

    def fit(temperature: NDArray[np.float64]) -> NDArray[np.float64]:
        above = np.maximum(temperature, last)
        continued = molar_heat_capacity(above) / molar_mass + remainder * (last / above) ** 2
        return np.where(temperature <= last, np.interp(temperature, point_temperature, point_value), continued)

    return fit


def make_wagner_fit(
    critical_temperature: float, critical_pressure: float, coefficients: Sequence[float], exponents: Sequence[float]
) -> PropertyFit:
    """Builds a Wagner-form vapor-pressure equation, the standard ancillary form for saturation pressure.

    The pressure p satisfies ln(p / p_c) = (T_c / T) * sum(N_i * theta^t_i), with theta = 1 - T / T_c. Above the
    critical temperature, the pressure is infinite, since no liquid can form.

    Args:
        critical_temperature: The critical temperature T_c in K.
        critical_pressure: The critical pressure p_c in Pa.
        coefficients: The coefficients N_i.
        exponents: The exponents t_i, one for each coefficient.

    Returns:
        A function from positive temperatures in K, of any shape, to pressures in Pa.

    Raises:
        ValueError: If the coefficients and exponents differ in number.
    """
    if len(coefficients) != len(exponents):
        raise ValueError(
            f"wagner takes the same number of coefficients and exponents; got {len(coefficients)} and {len(exponents)}."
        )

    n = np.array(coefficients, dtype=np.float64)
    t = np.array(exponents, dtype=np.float64)

    def fit(temperature: NDArray[np.float64]) -> NDArray[np.float64]:
        theta = np.maximum(1.0 - temperature / critical_temperature, 0.0)
        log_ratio = critical_temperature / temperature * np.sum(n * theta[..., np.newaxis] ** t, axis=-1)
        return np.where(temperature > critical_temperature, np.inf, critical_pressure * np.exp(log_ratio))

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
    "measured_with_lattice": make_measured_with_lattice_fit,
    "rational_poly": make_rational_poly_fit,
    "table": make_fit_from_table,
    "wagner": make_wagner_fit,
}


def make_fit(spec: Mapping[str, Any], label: str) -> PropertyFit:
    """Builds a property table's fit, held at the nearest end of its range, with a warning, outside that range.

    Args:
        spec: The property's table.
        label: The item and property names, such as "N2 thermal_conductivity", for warnings.
    """
    arguments = {key: value for key, value in spec.items() if key not in {"form", "range", "source"}}
    fit = PROPERTY_FITS[spec["form"]](**arguments)
    lower, upper = valid_range = tuple(spec["range"])

    def clamped_fit(temperature: NDArray[np.float64]) -> NDArray[np.float64]:
        warn_below(temperature, valid_range, label, f"using its value at {lower:g} K")
        warn_above(temperature, valid_range, label, f"using its value at {upper:g} K")
        return fit(np.clip(temperature, lower, upper))

    return clamped_fit


def warn_below(temperature: NDArray[np.float64], valid_range: Sequence[float], label: str, consequence: str) -> None:
    """Warns, on behalf of the property's caller, for temperatures below the range."""
    if np.any(temperature < valid_range[0]):
        warnings.warn(
            f"{label} evaluated below {_range_text(valid_range)}; {consequence}.", RuntimeWarning, stacklevel=3
        )


def warn_above(temperature: NDArray[np.float64], valid_range: Sequence[float], label: str, consequence: str) -> None:
    """Warns, on behalf of the property's caller, for temperatures above the range."""
    if np.any(temperature > valid_range[1]):
        warnings.warn(
            f"{label} evaluated above {_range_text(valid_range)}; {consequence}.", RuntimeWarning, stacklevel=3
        )


def _range_text(valid_range: Sequence[float]) -> str:
    return f"its range [{valid_range[0]:g}, {valid_range[1]:g}] K"
