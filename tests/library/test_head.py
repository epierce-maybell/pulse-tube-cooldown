"""Requirements for one head's gas inventory and mass swing."""

import math

import numpy as np
import pytest
from scipy.integrate import quad

from pulse_tube_cooldown.core.types import Material
from pulse_tube_cooldown.library.gases import helium
from pulse_tube_cooldown.library.heads import HeadModel, MatrixLayer, PulseTube, RegeneratorSection
from pulse_tube_cooldown.library.links.cooler.head import calc_gas_mass, calc_mass_swing

_MOLAR_MASS = 4.002602e-3  # kg/mol
_GAS_CONSTANT = 8.314462618  # J/(mol K)
_SECOND_VIRIAL_AT_300_K = 11.8e-6  # m^3/mol
_MATRIX = Material(
    name="Matrix", density=1.0, specific_heat=lambda t: t, thermal_conductivity=None, conductivity_integral=None
)


def _model(reservoir_volume=5.0e-3):
    """A head whose stage-2 section's second layer splits into two nodes, giving four nodes in all."""

    def section(*layers):
        return RegeneratorSection(
            length=0.2,
            wall_material=_MATRIX,
            wall_area=1e-4,
            layers=tuple(
                MatrixLayer(material=_MATRIX, mass=1.0, pore_volume=volume, length_fraction=fraction, node_count=count)
                for volume, fraction, count in layers
            ),
        )

    def tube(volume):
        return PulseTube(
            length=0.3, volume=volume, wall_material=_MATRIX, wall_area=1e-5, reservoir_volume=reservoir_volume
        )

    return HeadModel(
        name="Test",
        regenerators=(section((1.0e-4, 1.0, 1)), section((3.0e-5, 0.5, 1), (4.0e-5, 0.5, 2))),
        pulse_tubes=(tube(6.0e-5), tube(2.0e-4)),
        clamp_fraction=0.5,
        dead_volume=4.5e-5,
        waveform_factor=math.sqrt(2),
        buffer_fraction=0.5,
        stroke_spread=0.1,
    )


# Chain order: warm flange, stage-1 section node, stage 1, stage-2 section's three nodes, stage 2.
_COLD = np.array([300.0, 60.0, 40.0, 20.0, 10.0, 6.0, 4.0])
_PORE = {1: 1.0e-4, 3: 3.0e-5, 4: 2.0e-5, 5: 2.0e-5}  # m^3 at each node's chain index
_TUBE = {2: 6.0e-5, 6: 2.0e-4}  # m^3 at each stage's chain index
_FREQUENCY = 1.2


def _mean_over_tube(function, warm, cold):
    """Mean of `function(T)` over a linear temperature profile from `warm` to `cold`, by adaptive quadrature."""
    return quad(function, cold, warm, epsrel=1e-10, limit=200)[0] / (warm - cold)


def test_gas_mass_at_uniform_temperature_fills_every_volume():
    # HEPAK's density at 300 K and 1.7 MPa times the pores, tubes, reservoirs, and dead volume.
    volume = 1.7e-4 + 2.6e-4 + 1.0e-2 + 4.5e-5

    assert calc_gas_mass(_model(), np.full(7, 300.0), 1.7e6) == pytest.approx(2.7069548683 * volume, rel=1e-8)


def test_gas_mass_places_each_volume_at_its_temperature():
    pressure = 1.4e6
    pores = sum(volume * float(helium.density(_COLD[index], pressure)) for index, volume in _PORE.items())
    tubes = sum(
        volume * _mean_over_tube(lambda t: float(helium.density(t, pressure)), _COLD[0], _COLD[index])
        for index, volume in _TUBE.items()
    )
    warm = (1.0e-2 + 4.5e-5) * float(helium.density(_COLD[0], pressure))

    assert calc_gas_mass(_model(), _COLD, pressure) == pytest.approx(pores + tubes + warm, rel=1e-5)


def test_mass_swing_in_ideal_gas_limit():
    # Uniform 300 K and 0.4 MPa, where helium is ideal to 0.3%. Pores swing isothermally with the second virial
    # coefficient; the dead volume, tubes, and reservoirs swing adiabatically from the mean state as an ideal gas with
    # gamma = 5/3, each reservoir by the swing reduced by sqrt(1 + (omega R B)^2) for B = V_r / (gamma P).
    mean_pressure, swing = 0.4e6, 0.4e6
    resistance = (1.8e7, 3.6e7)  # Pa s/m^3, omega R B near 1 and 2
    mean_density = 0.6407114641  # HEPAK at 300 K and 0.4 MPa

    def isothermal(volume):
        def virial_density(pressure):
            return _MOLAR_MASS / (_GAS_CONSTANT * 300.0 / pressure + _SECOND_VIRIAL_AT_300_K)

        return volume * (virial_density(mean_pressure + swing / 2) - virial_density(mean_pressure - swing / 2))

    def adiabatic(volume, swing=swing):
        upper, lower = (1 + swing / (2 * mean_pressure)), (1 - swing / (2 * mean_pressure))
        return volume * mean_density * (upper**0.6 - lower**0.6)

    response = 2 * math.pi * _FREQUENCY * np.array(resistance) * 5.0e-3 / (5 / 3 * mean_pressure)
    reservoir = [adiabatic(5.0e-3, swing / math.sqrt(1 + r**2)) for r in response]
    expected = [
        adiabatic(4.5e-5),
        isothermal(1.0e-4),
        adiabatic(6.0e-5) + reservoir[0],
        isothermal(3.0e-5),
        isothermal(2.0e-5),
        isothermal(2.0e-5),
        adiabatic(2.0e-4) + reservoir[1],
    ]

    mass_swing = calc_mass_swing(
        _model(), np.full(7, 300.0), mean_pressure, swing, resistance=resistance, frequency=_FREQUENCY
    )

    np.testing.assert_allclose(mass_swing, expected, rtol=3e-3)


def test_mass_swing_with_cold_pores_and_tubes():
    # Without reservoirs: pores swing at their node temperatures, and each tube's gas adiabatically from the mean
    # state along a linear profile from the warm flange to its stage, averaged by the default quadrature to within
    # 0.1% of adaptive quadrature.
    mean_pressure, swing = 1.4e6, 1.2e6
    upper, lower = mean_pressure + swing / 2, mean_pressure - swing / 2

    def adiabatic(t):
        return float(
            helium.isentropic_density(t, mean_pressure, upper) - helium.isentropic_density(t, mean_pressure, lower)
        )

    expected = np.zeros(7)
    expected[0] = 4.5e-5 * adiabatic(300.0)
    for index, volume in _PORE.items():
        expected[index] = volume * float(helium.density(_COLD[index], upper) - helium.density(_COLD[index], lower))
    for index, volume in _TUBE.items():
        expected[index] = volume * _mean_over_tube(adiabatic, _COLD[0], _COLD[index])

    mass_swing = calc_mass_swing(
        _model(reservoir_volume=0.0), _COLD, mean_pressure, swing, resistance=(1e9, 1e9), frequency=_FREQUENCY
    )

    np.testing.assert_allclose(mass_swing, expected, rtol=1e-3)
