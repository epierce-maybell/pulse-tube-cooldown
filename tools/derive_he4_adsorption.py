"""Prints Langmuir constants for He-4 on charcoal from the source isotherm recorded in he4.toml.

The Roubeau isotherm is ln P = (a - b ln(Q + c)) (d - 1 / T) + e + ln Q + ln 100, for pressure P in Pa, temperature T in
K, and amount adsorbed Q in cm^3/g. The Langmuir constant K = (n / (n_s - n)) / P, for amount adsorbed n and capacity
n_s, gives 1 / T = (R / E) (ln K - ln K_0), fitted with weights T^2 over a grid of amounts and pressures.
"""

import tomllib
from pathlib import Path

import numpy as np
from scipy.constants import R


def main() -> None:
    """Prints adsorption_energy in J/mol, adsorption_constant in 1/Pa, and adsorption_capacity in mol/kg."""
    source = Path(__file__).resolve().parents[1] / "src/pulse_tube_cooldown/library/gases/data/he4.toml"
    data = tomllib.loads(source.read_text(encoding="utf-8"))["sources"]["martins2011"]
    mmol_per_cm3 = data["amount_pressure"] / (R * data["amount_temperature"]) * 1e-3

    volume_per_gram, log_pressure = np.meshgrid(
        np.linspace(*data["amount_range"], 60), np.linspace(*np.log(data["pressure_range"]), 60)
    )
    heat_over_gas_constant = data["a"] - data["b"] * np.log(volume_per_gram + data["c"])
    log_pressure_over_amount = log_pressure - np.log(volume_per_gram) - np.log(100.0)
    temperature = 1.0 / (data["d"] - (log_pressure_over_amount - data["e"]) / heat_over_gas_constant)

    amount = volume_per_gram * mmol_per_cm3
    capacity = data["capacity"] * mmol_per_cm3
    log_langmuir_constant = np.log(amount / (capacity - amount)) - log_pressure
    slope, intercept = np.polyfit(
        log_langmuir_constant.ravel(), 1.0 / temperature.ravel(), 1, w=temperature.ravel() ** 2
    )
    residual = 1.0 / (slope * log_langmuir_constant + intercept) - temperature

    print(f"adsorption_energy = {R / slope:.4g}")
    print(f"adsorption_constant = {np.exp(-intercept / slope):.4g}")
    print(f"adsorption_capacity = {capacity:.4g}")
    print(f"temperature residual: {np.sqrt(np.mean(residual**2)):.2f} K rms, {np.abs(residual).max():.2f} K max")


if __name__ == "__main__":
    main()
