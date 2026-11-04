"""Tests for the drive of a gas manifold: mean pressure from the gas charge and swing from the drive constraint."""

import numpy as np
import pytest

from pulse_tube_cooldown.library.links.cooler.drive import calc_mean_pressure, calc_pressure_swing
from pulse_tube_cooldown.library.links.cooler.head import calc_gas_mass, calc_mass_flow

# A head with one stage-1 layer and two stage-2 layers, in chain order: warm flange, stage-1 column, stage 1,
# stage-2 column, stage 2. It holds 4.1e-4 m^3 of gas.
_TEMPERATURE = np.array([300.0, 60.0, 40.0, 20.0, 8.0, 4.0])
_PORE_VOLUME = (np.array([1.0e-4]), np.array([3.0e-5, 2.0e-5]))
_PULSE_TUBE_VOLUME = (6.0e-5, 2.0e-4)
_FREQUENCY = 1.2

# One F-70H, from data/compressors.md: the swept displacement rate and the middle of the relief swing's range.
_DISPLACEMENT_RATE = 5.72e-3
_RELIEF_SWING = 1.93e6
_WARM_TEMPERATURE = 300.0
_WARM_VOLUME = 0.03
# Filled to 1.7 MPa with the head and warm volume at 293.15 K, where HEPAK's density is 2.7696416 kg/m^3.
_CHARGE = 2.7696416 * (4.1e-4 + _WARM_VOLUME)


def _gas_mass(temperature, pore_volume=_PORE_VOLUME, pulse_tube_volume=_PULSE_TUBE_VOLUME):
    return lambda pressure: calc_gas_mass(temperature, pore_volume, pulse_tube_volume, pressure)


def _head_flow(*temperatures, pore_volume=_PORE_VOLUME, pulse_tube_volume=_PULSE_TUBE_VOLUME, frequency=_FREQUENCY):
    def head_flow(supply_pressure, return_pressure):
        return sum(
            calc_mass_flow(temperature, pore_volume, pulse_tube_volume, supply_pressure, return_pressure, frequency)[0]
            for temperature in temperatures
        )

    return head_flow


def _swing(head_flow, mean_pressure=1.5e6, displacement_rate=(_DISPLACEMENT_RATE,), relief_swing=(_RELIEF_SWING,)):
    return calc_pressure_swing(
        head_flow,
        mean_pressure=mean_pressure,
        warm_temperature=_WARM_TEMPERATURE,
        displacement_rate=displacement_rate,
        relief_swing=relief_swing,
    )


def test_mean_pressure_at_rest_is_fill_pressure():
    mean_pressure = calc_mean_pressure(
        _gas_mass(np.full(6, 293.15)), charge=_CHARGE, warm_volume=_WARM_VOLUME, warm_temperature=293.15
    )

    assert mean_pressure == pytest.approx(1.7e6, rel=1e-6)


def test_mean_pressure_falls_as_head_cools():
    # The head cold and the warm volume at 300 K, from a root find on HEPAK densities directly.
    mean_pressure = calc_mean_pressure(
        _gas_mass(_TEMPERATURE), charge=_CHARGE, warm_volume=_WARM_VOLUME, warm_temperature=_WARM_TEMPERATURE
    )

    assert mean_pressure == pytest.approx(1.6057806e6, rel=1e-6)


def test_swing_meets_drive_constraint_below_relief_cap():
    # From a root find of f (M(P + dP/2) - M(P - dP/2)) = rho(300 K, P - dP/2) V on HEPAK densities directly, for the
    # head's gas mass M at a pressure.
    swing = _swing(_head_flow(_TEMPERATURE))

    assert swing == pytest.approx(1.5125596e6, rel=1e-6)


def test_warm_head_holds_swing_at_smallest_relief_cap():
    # At 300 K and a 1.9 MPa swing the head takes 1.5 g/s against the two compressors' 13.7 g/s.
    swing = _swing(
        _head_flow(np.full(6, 300.0)),
        mean_pressure=1.7e6,
        displacement_rate=[_DISPLACEMENT_RATE] * 2,
        relief_swing=[2.0e6, 1.9e6],
    )

    assert swing == 1.9e6


def test_two_heads_on_one_compressor_share_its_flow():
    # Two identical heads on one compressor see the swing that one of them sees on half the displacement.
    shared = _swing(_head_flow(_TEMPERATURE, _TEMPERATURE))
    alone = _swing(_head_flow(_TEMPERATURE), displacement_rate=[_DISPLACEMENT_RATE / 2])

    assert shared < _RELIEF_SWING
    assert shared == pytest.approx(alone, rel=1e-6)


@pytest.mark.parametrize(
    ("temperature", "one_compressor_swing"),
    [
        # Stages at 60 K and 20 K: both drives at the cap.
        ([300.0, 213.6, 93.6, 60.0, 56.8, 47.2, 34.4, 24.0, 20.0], _RELIEF_SWING),
        # Stages at 50 K and 10 K: one compressor has left the cap, by root finds on HEPAK densities directly.
        ([300.0, 210.0, 85.0, 50.0, 46.8, 37.2, 24.4, 14.0, 10.0], 1.8523216e6),
    ],
)
def test_two_compressors_hold_relief_cap_colder_than_one(temperature, one_compressor_swing):
    # PT1.1's 6.8024e-4 m^3 of gas volume, one node per layer at the layer's midpoint on a linear profile, filled to
    # 1.7 MPa at 293.15 K with 0.03 m^3 of warm volume, which sits at 300 K while running.
    pore_volume = (np.array([2.848e-4, 6.012e-5]), np.array([1.196e-5, 2.392e-5, 2.392e-5, 1.495e-5]))
    pulse_tube_volume = (5.967e-5, 2.009e-4)
    temperature = np.array(temperature)
    mean_pressure = calc_mean_pressure(
        _gas_mass(temperature, pore_volume, pulse_tube_volume),
        charge=2.7696416 * (6.8024e-4 + _WARM_VOLUME),
        warm_volume=_WARM_VOLUME,
        warm_temperature=_WARM_TEMPERATURE,
    )
    head_flow = _head_flow(temperature, pore_volume=pore_volume, pulse_tube_volume=pulse_tube_volume, frequency=1.0)
    one_compressor = _swing(head_flow, mean_pressure)
    two_compressors = _swing(head_flow, mean_pressure, [_DISPLACEMENT_RATE] * 2, [_RELIEF_SWING] * 2)

    assert one_compressor == pytest.approx(one_compressor_swing, rel=1e-5)
    assert two_compressors == _RELIEF_SWING
