"""Screen ADR 04 relaxation closures using SN09 geometry and HEPAK properties.

Run from the repository root with .venv/Scripts/python followed by this path.
Results are printed as JSON. Local losses use a prescribed linear 50--300 K
profile and harmonic flow; they are screening estimates, not cold-head fits.
"""

import csv
import json
import tomllib
from pathlib import Path

import numpy as np
from hepak.properties import properties_TP
from scipy.special import iv, roots_legendre

from pulse_tube_cooldown.library.materials.load import load_material


def laminar_loss_factor(radius, omega, diffusivity, viscosity, order=128):
    """Return -mean(u*theta) for unit mean velocity amplitude and unit gradient.

    The local linear radial equations have no-slip velocity and a fixed wall
    temperature. Complex amplitudes use exp(i*omega*t). Diffusivity and viscosity
    are respectively thermal diffusivity and kinematic viscosity, in m^2/s.
    """
    nodes, weights = roots_legendre(order)
    radial = (nodes + 1) / 2
    weights = weights / 2
    thermal_argument = radius * np.sqrt(1j * omega / diffusivity)
    viscous_argument = radius * np.sqrt(1j * omega / viscosity)
    thermal_profile = iv(0, thermal_argument * radial) / iv(0, thermal_argument)
    viscous_profile = iv(0, viscous_argument * radial) / iv(0, viscous_argument)
    viscous_average = 2 * iv(1, viscous_argument) / (viscous_argument * iv(0, viscous_argument))
    velocity = (1 - viscous_profile) / (1 - viscous_average)
    temperature = -(1 - thermal_profile - (viscous_profile - thermal_profile) / (1 - diffusivity / viscosity))
    temperature /= 1j * omega * (1 - viscous_average)
    return float(-np.sum(weights * radial * np.real(velocity * temperature.conjugate())))


def main():
    """Print geometry, pressure sensitivity, local closure and flow screens."""
    root = Path(__file__).resolve().parents[3]
    inputs = tomllib.loads(Path(__file__).with_suffix(".toml").read_text())
    with (root / "data/sn09/pulse_tube/geometry.csv").open(newline="") as stream:
        geometry = {row["quantity"]: row for row in csv.DictReader(stream)}
    radius = (
        float(geometry["stage1_pulse_tube_outer_radius"]["value"]) - float(geometry["tube_wall_thickness"]["value"])
    ) * 1e-3
    length = float(geometry["stage1_pulse_tube_length"]["value"]) * 1e-3
    wall_thickness = float(geometry["tube_wall_thickness"]["value"]) * 1e-3
    frequency = inputs["frequency_Hz"]
    omega = 2 * np.pi * frequency
    pressure = inputs["pressure_Pa"]
    cold, warm = inputs["cold_K"], inputs["warm_K"]
    transition = inputs["transition_Re_delta"]
    renewal_length = inputs["renewal_length_diameters"] * 2 * radius
    gradient = (warm - cold) / length
    area = np.pi * radius**2
    material = load_material(inputs["wall_material"])
    output = {
        "assumptions": {
            "radius_m": radius,
            "length_m": length,
            "frequency_Hz": frequency,
            "pressure_Pa": pressure,
            "cold_K": cold,
            "warm_K": warm,
            "transition_Re_delta": transition,
            "renewal_length_m": renewal_length,
            "reference_cooling_W": inputs["reference_cooling_W"],
            "turbulent_surrogate": "Nu = 0.023 Re_D^0.8 Pr^0.4 (steady-flow screening only)",
            "renewal_surrogate": "alpha = 4/R * sqrt(chi*U/(pi*ell)), fixed contact length ell",
        },
        "sources": inputs["sources"],
        "thermal": [],
        "closures_at_cold_temperature": [],
        "tube_compliance": [],
    }
    for temperature in inputs["temperatures_K"]:
        p = properties_TP(temperature, pressure, basis="mass")
        rho, cp, k, mu = (p[name] for name in ["density", "cp", "thermal_conductivity", "viscosity"])
        chi, nu = k / (rho * cp), mu / rho
        delta = np.sqrt(2 * chi / omega)
        alpha = omega * delta / radius
        exact = laminar_loss_factor(radius, omega, chi, nu)
        fine = laminar_loss_factor(radius, omega, chi, nu, order=256)
        wall_cp = float(material.specific_heat(temperature))
        wall_k = float(material.thermal_conductivity(temperature))
        wall_chi = wall_k / (material.density * wall_cp)
        wall_depth = min(wall_thickness, np.sqrt(2 * wall_chi / omega))
        wall_capacity = material.density * wall_cp * wall_depth
        output["thermal"].append(
            {
                "temperature_K": temperature,
                "alpha_over_omega": alpha / omega,
                "exact_laminar_loss_per_U2_W_s2_m2": area * rho * cp * gradient * exact,
                "lumped_loss_per_U2_W_s2_m2": area * rho * cp * gradient * alpha / (2 * (alpha**2 + omega**2)),
                "quadrature_relative_difference": abs(fine / exact - 1),
                "transition_U_m_s": transition * np.sqrt(nu * omega / 2),
                "transition_xi_over_length": transition * np.sqrt(nu * omega / 2) / (omega * length),
                "gas_layer_to_wall_heat_capacity_ratio": rho * cp * delta / wall_capacity,
                "wall_delta_m": np.sqrt(2 * wall_chi / omega),
            }
        )
        if temperature != cold:
            continue
        for velocity in inputs["velocities_m_s"]:
            diameter = 2 * radius
            reynolds = velocity * diameter / nu
            prandtl = nu / chi
            alphas = {
                "laminar_lumped": (alpha, 0.0),
                "turbulent_surrogate": (4 * chi / diameter**2 * 0.023 * reynolds**0.8 * prandtl**0.4, 0.8),
                "fixed_length_renewal": (4 / radius * np.sqrt(chi * velocity / (np.pi * renewal_length)), 0.5),
            }
            row = {
                "U_m_s": velocity,
                "Re_delta": velocity * np.sqrt(2 / (nu * omega)),
                "xi_over_length": velocity / (omega * length),
                "gradient_excursion_over_T": gradient * velocity / (omega * temperature),
                "exact_laminar_loss_W": area * rho * cp * gradient * exact * velocity**2,
                "models": {},
            }
            for name, (rate, exponent) in alphas.items():
                loss = area * rho * cp * gradient * velocity**2 * rate / (2 * (rate**2 + omega**2))
                slope = 2 + exponent * (omega**2 - rate**2) / (omega**2 + rate**2)
                row["models"][name] = {
                    "alpha_over_omega": rate / omega,
                    "loss_W": loss,
                    "loss_exponent": slope,
                    "slope_margin_at_reference_load_W": (slope - 2) * loss - 2 * inputs["reference_cooling_W"],
                }
            output["closures_at_cold_temperature"].append(row)
    nodes, weights = roots_legendre(128)
    temperatures = cold + (nodes + 1) * (warm - cold) / 2
    compressibility_integral = 0.0
    for temperature, weight in zip(temperatures, weights, strict=True):
        p = properties_TP(float(temperature), pressure, basis="mass")
        compressibility_integral += weight / p["speed_of_sound"] ** 2
    compressibility_integral *= area * length / 2
    cold_density = properties_TP(cold, pressure, basis="mass")["density"]
    for amplitude in inputs["pressure_amplitudes_Pa"]:
        mass_amplitude = omega * amplitude * compressibility_integral
        output["tube_compliance"].append(
            {
                "pressure_amplitude_Pa": amplitude,
                "cold_mass_amplitude_kg_s": mass_amplitude,
                "cold_U_m_s": mass_amplitude / (area * cold_density),
            }
        )
    output["pressure_sensitivity_at_cold_temperature"] = []
    for mean_pressure in inputs["mean_pressures_Pa"]:
        p = properties_TP(cold, mean_pressure, basis="mass")
        chi = p["thermal_conductivity"] / (p["density"] * p["cp"])
        nu = p["viscosity"] / p["density"]
        output["pressure_sensitivity_at_cold_temperature"].append(
            {
                "pressure_Pa": mean_pressure,
                "alpha_over_omega": np.sqrt(2 * chi / omega) / radius,
                "transition_U_m_s": transition * np.sqrt(nu * omega / 2),
            }
        )
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
