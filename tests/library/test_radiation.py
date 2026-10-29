"""Tests for gray-body radiation between two surfaces."""

import numpy as np
import pytest

from pulse_tube_cooldown.core.errors import ConfigurationError
from pulse_tube_cooldown.core.types import Surface
from pulse_tube_cooldown.library.links.radiation import Radiation

SIGMA = 5.670374419184429e-8  # W m^-2 K^-4, exact in the 2019 SI


def _constant(emissivity: float) -> Surface:
    return Surface(name=f"Constant {emissivity}", emissivity=lambda t: emissivity + 0.0 * t)


_PLATING = (_constant(0.05), _constant(0.05))


def _shield(area=(0.1, 0.1), surface=_PLATING) -> Radiation:
    return Radiation(name="Shield", nodes=("STP", "CFP"), area=area, surface=surface)


def _network_flow(t_inner, t_outer, a_inner, a_outer, e_inner, e_outer):
    """Heat into the inner surface from the radiosity network with a unit view factor from inner to outer."""
    resistance = (1 - e_inner) / (e_inner * a_inner) + 1 / a_inner + (1 - e_outer) / (e_outer * a_outer)
    return SIGMA * (t_outer**4 - t_inner**4) / resistance


def test_parallel_plates_match_hand_calculation():
    # Q = sigma A (T2^4 - T1^4) / (1/e1 + 1/e2 - 1) = 5.670374419e-8 * 0.1 * (300^4 - 4^4) / 39 = 1.17769 W.
    flow = _shield().calc_heat_flow(np.array([4.0, 300.0]))

    np.testing.assert_allclose(flow, [1.177693, -1.177693], rtol=1e-6)


def test_unequal_enclosed_areas_match_radiosity_network():
    area, emissivity = (0.119, 0.45), (0.03, 0.08)
    radiation = _shield(area=area, surface=(_constant(emissivity[0]), _constant(emissivity[1])))

    flow = radiation.calc_heat_flow(np.array([45.0, 70.0]))

    expected = _network_flow(45.0, 70.0, *area, *emissivity)
    np.testing.assert_allclose(flow, [expected, -expected], rtol=1e-12)


def test_emissivity_is_evaluated_at_each_surface_temperature():
    inner = Surface(name="Rising", emissivity=lambda t: t / 1000.0)
    outer = Surface(name="Falling", emissivity=lambda t: 0.3 - t / 1000.0)
    radiation = _shield(area=(0.2, 0.5), surface=(inner, outer))

    flow = radiation.calc_heat_flow(np.array([50.0, 200.0]))

    expected = _network_flow(50.0, 200.0, 0.2, 0.5, 0.05, 0.1)
    np.testing.assert_allclose(flow, [expected, -expected], rtol=1e-12)


def test_equal_temperatures_carry_no_heat():
    flow = _shield(area=(0.1, 0.3)).calc_heat_flow(np.array([77.0, 77.0]))

    np.testing.assert_array_equal(flow, [0.0, 0.0])


def test_warmer_inner_surface_loses_heat_to_outer_surface():
    radiation = _shield(area=(0.1, 0.3), surface=(_constant(0.05), _constant(0.1)))

    forward = radiation.calc_heat_flow(np.array([40.0, 250.0]))
    reverse = radiation.calc_heat_flow(np.array([250.0, 40.0]))

    assert reverse[0] < 0.0 < reverse[1]
    assert reverse.sum() == pytest.approx(0.0, abs=1e-12 * abs(reverse[0]))
    np.testing.assert_allclose(reverse, -forward, rtol=1e-12)


def test_rejects_inner_area_larger_than_outer_area():
    with pytest.raises(ConfigurationError, match="Shield"):
        _shield(area=(0.5, 0.1))
