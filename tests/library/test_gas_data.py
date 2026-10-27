"""Tests for the library's gas definitions against their references."""

import numpy as np
import pytest

from pulse_tube_cooldown.library.gases.load import load_gas


def test_nitrogen_constants_are_si():
    nitrogen = load_gas("N2")

    assert nitrogen.molar_mass == pytest.approx(28.0134e-3, rel=1e-5)
    assert nitrogen.heat_capacity_ratio == pytest.approx(1.4, rel=1e-3)


def test_nitrogen_thermal_conductivity_agrees_with_nist_webbook():
    nitrogen = load_gas("N2")

    # NIST Chemistry WebBook isobar at 1 kPa, in W/(m K).
    temperature = np.array([63.151, 73.151, 103.15, 153.15, 203.15, 273.15, 293.15, 303.15])
    expected = [0.0056026, 0.0066171, 0.0095838, 0.014220, 0.018491, 0.023966, 0.025440, 0.026164]

    np.testing.assert_allclose(nitrogen.thermal_conductivity(temperature), expected, rtol=2e-3)


def test_nitrogen_saturation_pressure_agrees_with_nist_webbook():
    nitrogen = load_gas("N2")

    # NIST Chemistry WebBook saturation table, in Pa, and the normal boiling point from Span et al. (2000).
    temperature = np.array([65.0, 70.0, 77.355, 80.0, 90.0, 100.0, 120.0])
    expected = [17404.4, 38544.8, 101325.0, 136872.0, 360458.0, 778275.0, 2510580.0]

    np.testing.assert_allclose(nitrogen.saturation_pressure(temperature), expected, rtol=1e-3)


def test_nitrogen_does_not_condense_above_its_critical_point():
    assert np.all(np.isposinf(load_gas("N2").saturation_pressure(np.array([130.0, 300.0]))))


def test_nitrogen_is_valid_from_triple_point_to_room_temperature(recwarn):
    nitrogen = load_gas("N2")
    temperature = np.linspace(63.151, 300.0, 50)

    nitrogen.thermal_conductivity(temperature)
    nitrogen.saturation_pressure(temperature)

    assert not recwarn.list


def test_nitrogen_accommodation_coefficient_matches_room_temperature_measurement():
    # Trott et al. (2006): 0.80 +/- 0.02 for N2 on 304 stainless steel near room temperature.
    np.testing.assert_allclose(load_gas("N2").accommodation_coeff(np.array([77.0, 295.0])), [0.80, 0.80])


_AWAITS_ADSORPTION = pytest.mark.xfail(
    raises=KeyError, strict=True, reason="He-4 charcoal adsorption constants are supplied by task 9.3."
)


@_AWAITS_ADSORPTION
def test_helium_constants_are_si():
    helium = load_gas("He-4")

    assert helium.molar_mass == pytest.approx(4.002602e-3, rel=1e-6)
    assert helium.heat_capacity_ratio == pytest.approx(5.0 / 3.0, rel=1e-6)


@_AWAITS_ADSORPTION
def test_helium_thermal_conductivity_agrees_with_nist_webbook():
    # NIST Chemistry WebBook isobar at 1 kPa, in W/(m K); about 1.3% above the ab initio value at 10 K.
    np.testing.assert_allclose(
        load_gas("He-4").thermal_conductivity(np.array([10.0, 300.0])), [0.0166343, 0.155896], rtol=1.5e-2
    )


@_AWAITS_ADSORPTION
def test_helium_saturation_pressure_agrees_with_its90():
    # ITS-90 He-4 vapor-pressure equation for 2.1768-5.0 K, inverted for pressure in Pa.
    temperature = np.array([2.1768, 2.2, 3.0, 4.2221, 5.0])
    expected = [5041.8, 5335.1, 24046.4, 101325.0, 196016.5]

    np.testing.assert_allclose(load_gas("He-4").saturation_pressure(temperature), expected, rtol=1e-3)


@_AWAITS_ADSORPTION
def test_helium_accommodation_coefficient_matches_cited_measurements():
    # Corruccini (1959) about 0.6 at 20 K; Trott et al. (2011) 0.43 at room temperature.
    np.testing.assert_allclose(load_gas("He-4").accommodation_coeff(np.array([20.0, 300.0])), [0.6, 0.43])


@_AWAITS_ADSORPTION
def test_helium_does_not_condense_above_its_critical_point():
    assert np.all(np.isposinf(load_gas("He-4").saturation_pressure(np.array([5.5, 300.0]))))
