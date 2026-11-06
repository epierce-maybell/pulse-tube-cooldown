"""Analytic ideal-gas checks of the isolated tube termination."""

import numpy as np
import pytest

from pulse_tube_cooldown.library.heads import load_head_model
from simulations.validation.test_stand_steady_state.phase_feasibility import balanced


@pytest.mark.parametrize("cold_temperature", [3.4, 58.0, 295.0])
def test_adiabatic_tube_admittance_and_power(monkeypatch, cold_temperature):
    """Constant-pressure tubes preserve acoustic power across a temperature gradient."""
    gamma, gas_constant = 5 / 3, 2077.0
    pressure, omega, warm_temperature = 1565067.2928375078, 2 * np.pi, 295.0

    def ideal_gas(temperature, pressure):
        return {
            "temperature": temperature,
            "density": pressure / (gas_constant * temperature),
            "cp": gamma * gas_constant / (gamma - 1),
            "expansion_coefficient": 1.0,
            "isothermal_compressibility": 1 / pressure,
            "gamma": gamma,
        }

    monkeypatch.setattr(balanced, "properties_TP", ideal_gas)
    model = load_head_model("PT1.5")
    temperature = np.array([warm_temperature, 200.0, 90.0, cold_temperature, cold_temperature])
    terminations = balanced.tube_terminations(model, temperature)
    cold_density = pressure / (gas_constant * cold_temperature)
    warm_density = pressure / (gas_constant * warm_temperature)
    for tube, resistance, termination in zip(model.pulse_tubes, [1e10, 3e10], terminations, strict=True):
        cold_y, warm_y, transfer, _ = termination
        reservoir_compliance = tube.reservoir_volume / (gamma * pressure)
        valve_y = 1 / (resistance + 1 / (1j * omega * reservoir_compliance))
        expected = cold_density * (valve_y + 1j * omega * tube.volume / (gamma * pressure))
        np.testing.assert_allclose(cold_y, expected, rtol=1e-9, atol=0)
        np.testing.assert_allclose(cold_y.real / cold_density, warm_y.real / warm_density, rtol=1e-9)
        warm_state = transfer @ np.array([1.0, cold_y])
        np.testing.assert_allclose(warm_state, [1.0, warm_y], rtol=1e-9, atol=1e-18)
