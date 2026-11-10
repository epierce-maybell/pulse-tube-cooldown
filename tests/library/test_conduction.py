"""Tests for conduction through a uniform cross-section."""

import numpy as np
import pytest
from scipy.integrate import quad

from pulse_tube_cooldown.core.types import Material
from pulse_tube_cooldown.library.links.conduction import Conduction
from pulse_tube_cooldown.library.materials.load import load_material

CONSTANT = Material(
    name="Constant",
    density=1.0,
    specific_heat=lambda t: 1.0 + 0.0 * t,
    specific_heat_integral=lambda t: 1.0 * t,
    thermal_conductivity=lambda t: 2.0 + 0.0 * t,
    conductivity_integral=lambda t: 2.0 * t,
)


def _rod(material: Material, area: float = 3.0e-4, length: float = 0.1) -> Conduction:
    return Conduction(name="Rod", nodes=("CFP", "STP"), material=material, area=area, length=length)


@pytest.mark.parametrize("temperature", [(4.0, 300.0), (4.0, 10.0), (50.0, 77.0), (120.0, 290.0)])
def test_temperature_dependent_conductivity_matches_independent_integration(temperature):
    material = load_material("304 SS")
    cold, hot = temperature
    area, length = 2.0e-5, 0.15

    integral, _ = quad(lambda t: material.thermal_conductivity(np.array([t]))[0], cold, hot, limit=200)
    flow = _rod(material, area=area, length=length).calc_energy_balance_terms(np.array([cold, hot])).heat_flow

    np.testing.assert_allclose(flow, [area / length * integral, -area / length * integral], rtol=1e-6)
