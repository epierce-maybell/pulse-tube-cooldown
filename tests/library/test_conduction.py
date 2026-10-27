"""Tests for conduction through a uniform cross-section."""

import math

import numpy as np
import pytest
from scipy.integrate import quad

from pulse_tube_cooldown.core.errors import ConfigurationError
from pulse_tube_cooldown.core.types import Material
from pulse_tube_cooldown.library.links.conduction import Conduction
from pulse_tube_cooldown.library.materials.load import load_material

CONSTANT = Material(
    name="Constant",
    specific_heat=lambda t: 1.0 + 0.0 * t,
    thermal_conductivity=lambda t: 2.0 + 0.0 * t,
    conductivity_integral=lambda t: 2.0 * t,
)


def _rod(material: Material, area: float = 3.0e-4, length: float = 0.1) -> Conduction:
    return Conduction(name="Rod", nodes=("CFP", "STP"), material=material, area=area, length=length)


def test_constant_conductivity_matches_analytic_flow():
    # Q = k A / L (T_STP - T_CFP) = 2 * 3e-4 / 0.1 * (300 - 10) = 1.74 W into CFP.
    flow = _rod(CONSTANT).calc_heat_flow(np.array([10.0, 300.0]))

    np.testing.assert_allclose(flow, [1.74, -1.74], rtol=1e-12)


@pytest.mark.parametrize("temperature", [(4.0, 300.0), (4.0, 10.0), (50.0, 77.0), (120.0, 290.0)])
def test_temperature_dependent_conductivity_matches_independent_integration(temperature):
    material = load_material("304 SS")
    cold, hot = temperature
    area, length = 2.0e-5, 0.15

    integral, _ = quad(lambda t: material.thermal_conductivity(np.array([t]))[0], cold, hot, limit=200)
    flow = _rod(material, area=area, length=length).calc_heat_flow(np.array([cold, hot]))

    np.testing.assert_allclose(flow, [area / length * integral, -area / length * integral], rtol=1e-6)


def test_reversed_temperatures_reverse_flow_and_conserve_energy():
    rod = _rod(load_material("OFHC Cu RRR100"))

    forward = rod.calc_heat_flow(np.array([40.0, 250.0]))
    reverse = rod.calc_heat_flow(np.array([250.0, 40.0]))

    assert forward[0] > 0.0
    assert forward.sum() == pytest.approx(0.0, abs=1e-12 * abs(forward[0]))
    assert reverse.sum() == pytest.approx(0.0, abs=1e-12 * abs(reverse[0]))
    np.testing.assert_allclose(reverse, -forward, rtol=1e-12)


def test_equal_temperatures_carry_no_heat():
    flow = _rod(load_material("304 SS")).calc_heat_flow(np.array([77.0, 77.0]))

    np.testing.assert_array_equal(flow, [0.0, 0.0])


@pytest.mark.parametrize("nodes", [("CFP",), ("CFP", "STP", "ICP")])
def test_rejects_other_than_two_nodes(nodes):
    with pytest.raises(ConfigurationError, match="Rod"):
        Conduction(name="Rod", nodes=nodes, material=CONSTANT, area=1.0e-4, length=0.1)  # type: ignore[arg-type]


@pytest.mark.parametrize("field", ["area", "length"])
@pytest.mark.parametrize("value", [0.0, -1.0, math.inf, math.nan])
def test_rejects_nonpositive_or_nonfinite_geometry(field, value):
    geometry = {"area": 1.0e-4, "length": 0.1, field: value}

    with pytest.raises(ConfigurationError, match=field):
        Conduction(name="Rod", nodes=("CFP", "STP"), material=CONSTANT, **geometry)
