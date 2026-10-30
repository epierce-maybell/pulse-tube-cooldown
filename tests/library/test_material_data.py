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


# Molar gas constant in J/(mol K) and the Avogadro constant in 1/mol, CODATA 2018.
_GAS_CONSTANT = 8.314462618
_AVOGADRO = 6.02214076e23

# Molar masses in kg/mol from IUPAC standard atomic weights.
_LEAD_MOLAR_MASS = 0.2072
_ER3NI_MOLAR_MASS = (3 * 167.259 + 58.6934) * 1e-3
_HOCU2_MOLAR_MASS = (164.93033 + 2 * 63.546) * 1e-3
_GOS_MOLAR_MASS = (2 * 157.25 + 2 * 15.999 + 32.06) * 1e-3

_REGENERATOR_MATERIALS = ["Pb", "Er3Ni", "HoCu2", "Gd2O2S"]


def _volumetric_heat_capacity(material, temperature):
    """Returns the heat capacity per unit volume in J/(cm^3 K), the unit of the regenerator literature."""
    return material.specific_heat(np.asarray(temperature, dtype=np.float64)) * material.density * 1e-6


def _er3ni_kuriyama(temperature):
    """Er3Ni volumetric heat capacity in J/(cm^3 K) below 30 K, as fitted in NIST REGEN3.3 (compute_mod.f90, material 9)
    to the data of T. Kuriyama (Toshiba), letter of 1990-12-21."""
    t = temperature
    if t > 20.0:
        return 0.04257 + 0.02552 * t
    if t > 7.136:
        return 1.0832 - 0.1667 * t + 0.01117 * t**2 - 0.0002082 * t**3
    if t > 4.0:
        return -0.225 + 0.08574 * t
    return 10 ** (2.5 * np.log10(t) - 2.436)


def _hocu2_regen(temperature):
    """HoCu2 volumetric heat capacity in J/(cm^3 K) below 20 K, as fitted in NIST REGEN3.3 (compute_mod.f90, material
    23)."""
    t = temperature
    if t > 9.76:
        return -6.6756744 + 1.405063 * np.log(t) + 8.5155067 / np.log(t)
    if t > 9.13:
        return 3.0188643 - 0.28247374 * t
    if t > 6.98:
        return 22.685639 - 8.1329417 * t + 0.97869285 * t**2 - 0.038857012 * t**3
    if t > 6.44:
        return 1 / (8217.5447 - 5110.8339 * t + 1192.1197 * t**2 - 123.55526 * t**3 + 4.8007214 * t**4)
    return (0.062912185 + 0.13235092 * t - 0.0065811947 * t**2) ** 2


def _gos_numazawa(temperature):
    """Gd2O2S volumetric heat capacity in J/(cm^3 K) below 10 K, as fitted in NIST REGEN3.3 (compute_mod.f90, material
    33) to Numazawa et al., Cryocoolers 12 (2003)."""
    t = temperature
    if t > 5.245:
        return (-0.059716577 + 0.026807868 * t**0.5) / (1 - 1.0024168 * t**0.5 + 0.24731957 * t)
    return (-0.16495511 + 0.31185708 * t - 0.042708217 * t**2) / (1 + 0.33018724 * t - 0.089944754 * t**2)


@pytest.mark.parametrize("name", _REGENERATOR_MATERIALS)
def test_regenerator_materials_are_valid_from_2_to_300_kelvin(recwarn, name):
    material = load_material(name)

    material.specific_heat(np.geomspace(2.0, 300.0, 50))

    assert not recwarn.list


@pytest.mark.parametrize(
    ("temperature", "molar_heat_capacity"),
    [
        # J. W. Arblaster, Calphad 39, 47 (2012), in J/(mol K). Table 11, superconducting state.
        (2.0, 13.92e-3),
        (4.0, 139.3e-3),
        (7.0, 1049e-3),
        # Table 12.
        (8.0, 1.527),
        (20.0, 10.956),
        (50.0, 21.459),
        (100.0, 24.302),
        (200.0, 25.707),
        # Table 13.
        (300.0, 26.508),
    ],
)
def test_lead_specific_heat_matches_arblaster_selected_values(temperature, molar_heat_capacity):
    material = load_material("Pb")

    np.testing.assert_allclose(
        material.specific_heat(np.array([temperature])), [molar_heat_capacity / _LEAD_MOLAR_MASS], rtol=1e-3
    )


@pytest.mark.parametrize(
    ("name", "reference", "temperatures"),
    [
        ("Er3Ni", _er3ni_kuriyama, [2.0, 3.0, 5.0, 7.0, 10.0, 15.0, 25.0]),
        ("HoCu2", _hocu2_regen, [2.0, 4.0, 6.0, 7.5, 9.5, 12.0, 18.0]),
        ("Gd2O2S", _gos_numazawa, [2.0, 3.0, 4.0, 4.8, 6.0, 8.0]),
    ],
)
def test_rare_earth_heat_capacity_follows_measurements_where_they_exist(name, reference, temperatures):
    material = load_material(name)

    expected = [reference(t) for t in temperatures]

    np.testing.assert_allclose(_volumetric_heat_capacity(material, temperatures), expected, rtol=0.01)


def test_gos_peaks_above_1_j_per_cm3_k_near_5_2_kelvin():
    material = load_material("Gd2O2S")
    temperature = np.linspace(4.0, 7.0, 3001)

    # Numazawa et al., Cryocoolers 12 (2003): more than 1 J/(cm^3 K) at 5.2 K.
    heat_capacity = _volumetric_heat_capacity(material, temperature)

    assert temperature[np.argmax(heat_capacity)] == pytest.approx(5.2, abs=0.1)
    assert heat_capacity.max() > 1.0


def test_hocu2_keeps_both_transitions():
    material = load_material("HoCu2")

    # Bischof et al., Phys. Status Solidi A 114, K229 (1989): anomalies at T_M = 6.7 K and T_N = 9 K, with the
    # heat capacity falling between and above them.
    lower = np.linspace(6.2, 7.2, 1001)
    upper = np.linspace(8.5, 9.5, 1001)
    lower_peak = lower[np.argmax(_volumetric_heat_capacity(material, lower))]
    upper_peak = upper[np.argmax(_volumetric_heat_capacity(material, upper))]
    between, above = _volumetric_heat_capacity(material, [7.8, 10.5])

    assert lower_peak == pytest.approx(6.7, abs=0.2)
    assert upper_peak == pytest.approx(9.0, abs=0.2)
    assert between < _volumetric_heat_capacity(material, [lower_peak])[0]
    assert between < _volumetric_heat_capacity(material, [upper_peak])[0]
    assert above < _volumetric_heat_capacity(material, [upper_peak])[0]


def test_er3ni_peaks_near_7_kelvin():
    material = load_material("Er3Ni")
    temperature = np.linspace(4.0, 12.0, 8001)

    # Kuriyama's data as fitted in REGEN3.3 peak at 7.136 K.
    heat_capacity = _volumetric_heat_capacity(material, temperature)

    assert temperature[np.argmax(heat_capacity)] == pytest.approx(7.1, abs=0.2)


@pytest.mark.parametrize(
    ("name", "molar_mass", "atoms_per_formula_unit"),
    [("Er3Ni", _ER3NI_MOLAR_MASS, 4), ("HoCu2", _HOCU2_MOLAR_MASS, 3), ("Gd2O2S", _GOS_MOLAR_MASS, 5)],
)
def test_rare_earth_heat_capacity_approaches_dulong_petit_at_room_temperature(name, molar_mass, atoms_per_formula_unit):
    material = load_material(name)

    # The Dulong-Petit limit, 3 R per atom.
    dulong_petit = 3 * atoms_per_formula_unit * _GAS_CONSTANT / molar_mass

    assert material.specific_heat(np.array([300.0]))[0] / dulong_petit == pytest.approx(1.0, abs=0.07)


@pytest.mark.parametrize(
    ("name", "expected"),
    [
        # Handbook values at room temperature, in kg/m^3.
        ("OFHC Cu RRR100", 8960.0),
        ("304 SS", 8000.0),
        ("Pb", 11340.0),
        # X-ray density, Z M / (N_A V), from the orthorhombic Imma cell a = 4.28, b = 6.76, c = 7.27 Angstrom with
        # Z = 4 of Smetana, Sima & Lebech (1986), COD 1525314.
        ("HoCu2", 4 * _HOCU2_MOLAR_MASS / (_AVOGADRO * 4.28 * 6.76 * 7.27e-30)),
        # X-ray density from the trigonal P-3m1 cell a = 3.851, c = 6.664 Angstrom with Z = 1 of Onoda et al.,
        # J. Solid State Chem. (2000), COD 1521699.
        ("Gd2O2S", _GOS_MOLAR_MASS / (_AVOGADRO * np.sqrt(3) / 2 * 3.851**2 * 6.664e-30)),
        # The sum of the elements' molar volumes, from Er 9066 and Ni 8908 kg/m^3.
        ("Er3Ni", _ER3NI_MOLAR_MASS / (3 * 0.167259 / 9066.0 + 0.0586934 / 8908.0)),
    ],
)
def test_material_density(name, expected):
    assert load_material(name).density == pytest.approx(expected, rel=2e-3)
