"""One bounded fixed-drive harmonic/mean-equilibrium diagnostic.

Steady Ergun and Wakao correlations are amplitude-based approximations for
oscillatory flow. All geometry assumptions are inherited from check.py.
"""

import json
import time
from pathlib import Path

import numpy as np
from hepak.properties import properties_TP
from scipy.integrate import quad
from scipy.optimize import least_squares
from scipy.special import ndtr

from pulse_tube_cooldown.library.heads import load_head_model
from simulations.validation.test_stand_steady_state.phase_feasibility.check import cell, closure, coefficients, profiles


def make_cells(model, temperature):
    """Creates shared faces at physical matrix-node boundaries."""
    sections = []
    for section, (point, ts) in zip(model.regenerators, profiles(model, temperature), strict=True):
        node_layers = [layer for layer in section.layers for _ in range(layer.node_count)]
        widths = section.length * np.array([layer.length_fraction / layer.node_count for layer in node_layers])
        faces = np.r_[0, np.cumsum(widths)]
        layer_edges = np.r_[0, np.cumsum([section.length * layer.length_fraction for layer in section.layers])]
        grid = np.unique(np.r_[faces, *[np.linspace(a, b, 5) for a, b in zip(point[:-1], point[1:], strict=True)]])
        entries = []
        for left, right in zip(grid[:-1], grid[1:], strict=True):
            if right - left < 1e-12:
                continue
            mid = (left + right) / 2
            interval = min(np.searchsorted(point, mid), len(point) - 1)
            layer = section.layers[min(np.searchsorted(layer_edges, mid) - 1, len(section.layers) - 1)]
            temp = np.interp(mid, point, ts)
            gradient = (ts[interval] - ts[interval - 1]) / (point[interval] - point[interval - 1])
            prop, area, r, h, cg, cm = closure(layer, section.length * layer.length_fraction, temp, 0.0003, 0.00004)
            phi = 0.7 if layer.material.name == "304 SS" else 0.38
            diameter = phi * 0.00004 / (1 - phi) if phi == 0.7 else 0.0003
            entries.append((left, right, prop, area, r, h, cg, cm, gradient, phi, diameter))
        sections.append((entries, faces, section, point, ts))
    return sections


def tube_terminations(model, temperature):
    """Ideal adiabatic pulse tubes, with existing warm valve/reservoir closure."""
    omega, pressure = 2 * np.pi, 1565067.2928375078
    results = []

    def compressibility(position, cold):
        prop = properties_TP(float(cold + position * (295 - cold)), pressure)
        return prop["isothermal_compressibility"] - prop["expansion_coefficient"] ** 2 / (
            prop["temperature"] * prop["density"] * prop["cp"]
        )

    for tube, cold, resistance in zip(model.pulse_tubes, [temperature[3], temperature[-1]], [1e10, 3e10], strict=True):
        tube_compliance = tube.volume * quad(compressibility, 0, 1, args=(cold,), epsabs=1e-15)[0]
        warm = properties_TP(295, pressure)
        cold_gas = properties_TP(float(cold), pressure)
        compliance = tube.reservoir_volume / (warm["gamma"] * pressure)
        warm_y = warm["density"] / (resistance + 1 / (1j * omega * compliance))
        cold_y = cold_gas["density"] * (warm_y / warm["density"] + 1j * omega * tube_compliance)
        total = np.array(
            [[1, 0], [-1j * omega * warm["density"] * tube_compliance, warm["density"] / cold_gas["density"]]],
            dtype=complex,
        )
        results.append((cold_y, warm_y, total, tube))
    return results


def solve_network(model, temperature, *, max_iterations=20, drive_amplitude=1e6):
    """Solves flow-dependent harmonic coefficients at prescribed mean temperatures."""
    sections = make_cells(model, temperature)
    tubes = tube_terminations(model, temperature)
    count = sum(len(s[0]) for s in sections)
    mass = np.full(count, 0.01)
    converged = False
    for _iteration in range(max_iterations):
        matrices, coeffs = [], []
        for entry, mass_amplitude in zip([e for s in sections for e in s[0]], mass, strict=True):
            left, right, prop, area, r, h, cg, cm, gradient, phi, diameter = entry
            bulk = area / phi
            re = mass_amplitude * diameter / (bulk * prop["viscosity"])
            conductivity = prop["thermal_conductivity"] or 0.024365244746865064
            pr = prop["cp"] * prop["viscosity"] / conductivity
            h *= (2 + 1.1 * re**0.6 * pr ** (1 / 3)) / 2
            r *= 1 + 8 / (3 * np.pi) * 1.75 / 150 / (1 - phi) * re
            a, b, _ = coefficients(prop, area, cg, cm, h, gradient)
            matrices.append(cell(prop, area, a, b, r, right - left, 1 / area))
            coeffs.append((a, b, h, r))
        split = len(sections[0][0])
        downstream = [None] * count
        admittance = tubes[1][0]
        for j in range(count - 1, -1, -1):
            if j == split - 1:
                admittance += tubes[0][0]
            downstream[j] = admittance
            t = matrices[j]
            admittance = (admittance * t[0, 0] - t[1, 0]) / (t[1, 1] - admittance * t[0, 1])
        state = np.array([drive_amplitude + 0j, drive_amplitude * admittance])
        inlet = state.copy()
        all_states, mass_new = [], []
        for j, t in enumerate(matrices):
            if j == split:
                state[1] -= tubes[0][0] * state[0]
            end_p = np.linalg.det(t) / (t[1, 1] - downstream[j] * t[0, 1]) * state[0]
            end = np.array([end_p, downstream[j] * end_p])
            all_states.append((state.copy(), end.copy()))
            mass_new.append(abs((state[1] + end[1]) / 2))
            state = end
        mass_new = np.array(mass_new)
        relative_change = float(np.max(abs(mass_new - mass) / np.maximum(mass_new, 1e-7)))
        if relative_change < 2e-5:
            converged = True
            break
        mass = 0.4 * mass + 0.6 * mass_new

    face_fluxes, sanity = [], []
    offset = 0
    for entries, faces, section, point, ts in sections:
        cell_fluxes, theta_ratios, theta_matrix_ratios, amplitude_locations = [], [], [], []
        for local, entry in enumerate(entries):
            left, right, prop, area, r, h, cg, cm, gradient, phi, diameter = entry
            index = offset + local
            a, b, h, r = coeffs[index]
            values = []
            for position, state in zip([left, right], all_states[index], strict=True):
                p, m = state
                theta = a * p + b * m
                theta_ratios.append(abs(theta) / prop["temperature"])
                theta_matrix_ratios.append(abs(h / (h + 1j * 2 * np.pi * cm) * theta) / prop["temperature"])
                amplitude_locations.append(
                    {
                        "x_from_section_warm_m": position,
                        "coefficient_mean_temperature_K": prop["temperature"],
                        "profile_temperature_at_position_K": float(np.interp(position, point, ts)),
                        "theta_g_peak_K": float(abs(theta)),
                        "pressure_temperature_peak_K": float(abs(a * p)),
                        "gradient_temperature_peak_K": float(abs(b * m)),
                        "mass_peak_kg_s": float(abs(m)),
                        "pressure_peak_Pa": float(abs(p)),
                        "gradient_K_per_m": float(gradient),
                        "gas_displacement_peak_m": float(abs(m) / (2 * np.pi * prop["density"] * area)),
                        "gas_capacity_J_per_m_K": float(cg),
                        "matrix_capacity_J_per_m_K": float(cm),
                    }
                )
                hhat = prop["cp"] * theta + (1 - prop["expansion_coefficient"]) * p / prop["density"]
                values.append(0.5 * np.real(hhat * np.conj(m)))
            cell_fluxes.append(values)
        # Each physical face has one shared flux, averaged from its two cells.
        f = []
        for face in faces:
            candidates = []
            for entry, flux in zip(entries, cell_fluxes, strict=True):
                if np.isclose(entry[0], face, rtol=0, atol=1e-11):
                    candidates.append(flux[0])
                if np.isclose(entry[1], face, rtol=0, atol=1e-11):
                    candidates.append(flux[1])
            interval = min(max(np.searchsorted(point, face), 1), len(point) - 1)
            conduct = (
                -section.wall_area
                * (
                    section.wall_material.conductivity_integral(ts[interval])
                    - section.wall_material.conductivity_integral(ts[interval - 1])
                )
                / (point[interval] - point[interval - 1])
            )
            f.append(float(np.mean(candidates) + conduct))
        face_fluxes.append(np.array(f))
        sanity.append(
            {
                "max_theta_g_over_T": float(max(theta_ratios)),
                "max_theta_m_over_T": float(max(theta_matrix_ratios)),
                "maximum_ratio_location": amplitude_locations[int(np.argmax(theta_ratios))],
                "largest_absolute_theta_g_peak_K": max(s["theta_g_peak_K"] for s in amplitude_locations),
            }
        )
        offset += len(entries)
    junction, cold = all_states[split - 1][1], all_states[-1][1]
    tube_power, losses, walls, warm_tube_mass = [], [], [], []
    for state, (y, ywarm, _transfer, tube), temp in zip(
        [junction, cold], tubes, [temperature[3], temperature[-1]], strict=True
    ):
        rho = properties_TP(float(temp), 1565067.2928375078)["density"]
        tube_power.append(float(0.5 * abs(state[0]) ** 2 * np.real(y) / rho))
        mwarm = abs(ywarm * state[0])
        warm_tube_mass.append(float(mwarm))
        mass_tube = tube.volume * np.mean(
            [properties_TP(float(t), 1565067.2928375078)["density"] for t in np.linspace(temp, 295, 32)]
        )
        stroke = 2 * mwarm / (2 * np.pi)
        width = 0.1 * mass_tube
        z = (stroke - 0.5 * mass_tube) / width
        positive = (stroke - 0.5 * mass_tube) * ndtr(z) + width * np.exp(-z * z / 2) / np.sqrt(2 * np.pi)
        delta_h = (
            properties_TP(295, 1565067.2928375078)["enthalpy"]
            - properties_TP(float(temp), 1565067.2928375078)["enthalpy"]
        )
        losses.append(float(delta_h * positive))
        walls.append(
            float(
                tube.wall_area
                * (tube.wall_material.conductivity_integral(295) - tube.wall_material.conductivity_integral(temp))
                / tube.length
            )
        )
    residual = np.r_[np.diff(face_fluxes[0]) * -1, np.diff(face_fluxes[1]) * -1]
    q = np.array(
        [
            face_fluxes[0][-1] - face_fluxes[1][0] - tube_power[0] + losses[0] + walls[0],
            face_fluxes[1][-1] - tube_power[1] + losses[1] + walls[1],
        ]
    )
    return {
        "residual": residual,
        "Q": q,
        "harmonic_converged": converged,
        "harmonic_iterations": _iteration + 1,
        "harmonic_relative_change": relative_change,
        "amplitude_sanity": sanity,
        "stage_pressure_over_drive": [abs(junction[0]) / drive_amplitude, abs(cold[0]) / drive_amplitude],
        "tube_power_W": tube_power,
        "breakthrough_W": losses,
        "wall_W": walls,
        "face_flux_W": [f.tolist() for f in face_fluxes],
        "inlet_mass_amplitude_kg_s": float(abs(inlet[1])),
    }


class BudgetExpiredError(Exception):
    """The single operating-point attempt reached its wall-time budget."""


def serializable(result):
    """Converts top-level arrays in a network result to JSON-compatible lists."""
    return {key: value.tolist() if isinstance(value, np.ndarray) else value for key, value in result.items()}


def main():
    """Runs one bounded mean-equilibrium solve from the archived baseline profile."""
    start = time.perf_counter()
    model = load_head_model("PT1.5")
    temperature = np.array(
        json.loads((Path(__file__).parent.parent / "coupling_cases.json").read_text())[0]["temperature"]
    )
    matrix_indices = np.r_[1:3, 4:15]
    initial = solve_network(model, temperature)
    print(json.dumps({"initial": serializable(initial)}, indent=2), flush=True)
    best = {"norm": np.inf, "temperature": temperature.copy(), "result": initial, "evaluations": 0}

    def residual(values):
        if time.perf_counter() - start > 55:
            raise BudgetExpiredError
        trial = temperature.copy()
        trial[matrix_indices] = values
        answer = solve_network(model, trial)
        best["evaluations"] += 1
        norm = float(np.max(abs(answer["residual"])))
        if norm < best["norm"]:
            best.update(norm=norm, temperature=trial, result=answer)
        return answer["residual"]

    status = "not started"
    try:
        solution = least_squares(
            residual,
            temperature[matrix_indices],
            bounds=(2.5, 295),
            max_nfev=30,
            ftol=1e-8,
            xtol=1e-8,
            gtol=1e-8,
            diff_step=1e-4,
        )
        status = solution.message
    except BudgetExpiredError:
        status = "55 s wall-time limit"
    result = best["result"]
    output = {
        "fixed_mean_pressure_Pa": 1565067.2928375078,
        "fixed_inlet_peak_Pa": 1e6,
        "status": status,
        "initial": serializable(initial),
        "best": serializable(result),
        "temperature_K": best["temperature"].tolist(),
        "mean_max_residual_W": best["norm"],
        "mean_evaluations": best["evaluations"],
        "runtime_s": time.perf_counter() - start,
        "stability": "not evaluated unless mean equilibrium and amplitude gates pass",
        "relaxed_slopes": "not evaluated unless mean equilibrium and amplitude gates pass",
    }
    Path(__file__).with_name("balanced_results.json").write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps(output, indent=2), flush=True)


if __name__ == "__main__":
    main()
