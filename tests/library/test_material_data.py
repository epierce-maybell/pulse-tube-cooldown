"""Tests for the library's material definitions against their references."""

import numpy as np
import pytest
from scipy.integrate import quad

from pulse_tube_cooldown.library.materials.load import load_material

# Molar mass of copper used by White & Collocott (1984), in kg/mol.
_COPPER_MOLAR_MASS = 0.06354


def _copper_reference_equation(temperature):
    """Returns copper's specific heat in J/(kg K) below 25 K from the reference equation of White & Collocott (1984).

    G. K. White and S. J. Collocott, "Heat Capacity of Reference Materials: Cu and W," J. Phys. Chem. Ref. Data 13,
    1251 (1984), Eq. (1): Cp = sum(A_i * T^(2i - 1)) mJ/(mol K) for i = 1 to 6.
    """
    a = [6.9434e-1, 4.7548e-2, 1.6314e-6, 9.4786e-8, -1.3639e-10, 5.3898e-14]
    molar = sum(a_i * temperature ** (2 * i + 1) for i, a_i in enumerate(a))
    return molar * 1e-3 / _COPPER_MOLAR_MASS


def test_copper_specific_heat_matches_nist_fit():
    material = load_material("OFHC Cu RRR100")
    temperature = np.array([4.0, 20.0, 77.0, 300.0])

    # NIST Cryogenic Material Properties, OFHC Copper (revised 2010-02-03), specific heat log-polynomial with
    # coefficients a through i = -1.91844, -0.15973, 8.61013, -18.996, 21.9661, -12.7328, 3.54322, -0.3797, 0.
    # Expected values in J/(kg K) evaluate NIST's equation with the legacy scalar `eval_log_poly` in branch
    # Advanced_Feature_SIM.
    expected = [0.09944090395343075, 7.50607290450624, 195.9208752033289, 389.40151139668626]

    np.testing.assert_allclose(material.specific_heat(temperature), expected, rtol=1e-12)


def test_copper_thermal_conductivity_matches_nist_fit():
    material = load_material("OFHC Cu RRR100")
    temperature = np.array([4.0, 20.0, 77.0, 300.0])

    # NIST Cryogenic Material Properties, OFHC Copper (revised 2010-02-03), RRR = 100 rational polynomial. Expected
    # values in W/(m K) evaluate NIST's equation with the legacy scalar `eval_rational_poly` in branch
    # Advanced_Feature_SIM.
    expected = [642.2969607429585, 2422.510264536503, 547.1996980793672, 396.3239590150121]

    np.testing.assert_allclose(material.thermal_conductivity(temperature), expected, rtol=1e-12)


@pytest.mark.parametrize(
    ("temperature", "expected", "tolerance"),
    [
        # Below 25 K, White & Collocott's reference equation; NIST states a fit error of 10% below 15 K and 5% above.
        (4.0, _copper_reference_equation(4.0), 0.10),
        (20.0, _copper_reference_equation(20.0), 0.05),
        # White & Collocott (1984), Table 2, Cp in J/(mol K) converted to J/(kg K).
        (100.0, 16.01 / _COPPER_MOLAR_MASS, 0.05),
        (200.0, 22.63 / _COPPER_MOLAR_MASS, 0.05),
        (300.0, 24.45 / _COPPER_MOLAR_MASS, 0.05),
    ],
)
def test_copper_specific_heat_agrees_with_measured_reference_within_fit_error(temperature, expected, tolerance):
    material = load_material("OFHC Cu RRR100")

    np.testing.assert_allclose(material.specific_heat(np.array([temperature])), [expected], rtol=tolerance)


def test_copper_is_valid_from_4_to_300_kelvin(recwarn):
    material = load_material("OFHC Cu RRR100")

    # NIST gives both fits from 4 K to 300 K; evaluation at the range ends does not warn, and beyond them it does.
    material.specific_heat(np.array([4.0, 300.0]))
    material.thermal_conductivity(np.array([4.0, 300.0]))
    assert not recwarn.list

    with pytest.warns(RuntimeWarning, match=r"OFHC Cu RRR100.*specific_heat.*below.*4.*300"):
        material.specific_heat(np.array([3.9]))
    with pytest.warns(RuntimeWarning, match=r"OFHC Cu RRR100.*thermal_conductivity.*above.*4.*300"):
        material.thermal_conductivity(np.array([301.0]))


def test_stainless_steel_specific_heat_matches_nist_fit():
    material = load_material("304 SS")
    temperature = np.array([4.0, 20.0, 77.0, 300.0])

    # NIST Cryogenic Material Properties, 304 Stainless (UNS S30400), specific heat log-polynomial with coefficients
    # a through i = 22.0061, -127.5528, 303.647, -381.0098, 274.0328, -112.9212, 24.7593, -2.239153, 0. Expected
    # values in J/(kg K) evaluate NIST's equation with the legacy scalar `eval_log_poly` in branch Advanced_Feature_SIM.
    expected = [2.0659856848217872, 13.452478061309394, 204.49302744915207, 469.448840674266]

    # The large alternating coefficients cancel, so different evaluation orders agree to about 1e-11.
    np.testing.assert_allclose(material.specific_heat(temperature), expected, rtol=1e-9)


def test_stainless_steel_thermal_conductivity_matches_nist_fit():
    material = load_material("304 SS")
    temperature = np.array([1.0, 4.0, 20.0, 77.0, 300.0])

    # NIST Cryogenic Material Properties, 304 Stainless (UNS S30400), thermal conductivity log-polynomial with
    # coefficients a through i = -1.4087, 1.3982, 0.2543, -0.6260, 0.2334, 0.4256, -0.4658, 0.1650, -0.0199. Expected
    # values in W/(m K) evaluate NIST's equation with the legacy scalar `eval_log_poly` in branch Advanced_Feature_SIM.
    expected = [0.03902114421141359, 0.2723961889664811, 2.168622406294646, 7.920651602258382, 15.308653824348552]

    np.testing.assert_allclose(material.thermal_conductivity(temperature), expected, rtol=1e-12)


@pytest.mark.parametrize(("temperature", "expected"), [(20.0, 16.3), (80.0, 349.0), (290.0, 3060.0)])
def test_stainless_steel_conductivity_integral_agrees_with_independent_compilation(temperature, expected):
    material = load_material("304 SS")

    # Conductivity integrals in W/m from 4.2 K for AISI 304: V. Parma, "Cryostat Design," CERN Accelerator School
    # (2014), arXiv:1501.07154, Table 3. Compilations of 304 conductivity from different samples differ by more than
    # NIST's 2% fit error, up to about 15% at low temperature, so this checks for gross transcription or unit errors.
    integral, _ = quad(lambda t: material.thermal_conductivity(np.array([t]))[0], 4.2, temperature, limit=200)

    assert integral == pytest.approx(expected, rel=0.15)


def test_stainless_steel_specific_heat_is_valid_from_4_and_conductivity_from_1_to_300_kelvin(recwarn):
    material = load_material("304 SS")

    # NIST gives the specific heat equation from 4 K to 300 K and the conductivity equation from 1 K to 300 K.
    material.specific_heat(np.array([4.0, 300.0]))
    material.thermal_conductivity(np.array([1.0, 300.0]))
    assert not recwarn.list

    with pytest.warns(RuntimeWarning, match=r"304 SS.*specific_heat.*below.*4.*300"):
        material.specific_heat(np.array([3.9]))
    with pytest.warns(RuntimeWarning, match=r"304 SS.*thermal_conductivity.*below.*1.*300"):
        material.thermal_conductivity(np.array([0.9]))
    with pytest.warns(RuntimeWarning, match=r"304 SS.*thermal_conductivity.*above.*1.*300"):
        material.thermal_conductivity(np.array([301.0]))
