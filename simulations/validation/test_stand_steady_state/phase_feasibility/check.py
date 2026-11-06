"""Bounded fixed-drive, prescribed-profile harmonic transport check.

Run as python -m simulations.validation.test_stand_steady_state.phase_feasibility.check.
Sphere and wire sizes are assumptions, not measured PT1.5 dimensions.
"""

import json
from pathlib import Path
import time

import numpy as np
from scipy.linalg import expm

from hepak.properties import properties_TP
from pulse_tube_cooldown.library.heads import load_head_model


def closure(layer, length, temperature, diameter, wire):
    """Returns r [Pa s/(kg m)], H [W/(m K)], Cg,Cm [J/(m K)]."""
    prop = properties_TP(float(temperature), 1565067.2928375078)
    conductivity = prop["thermal_conductivity"]
    if temperature < 4:
        conductivity = properties_TP(4, 1565067.2928375078)["thermal_conductivity"]
    pore_area = layer.pore_volume / length
    screen = layer.material.name == "304 SS"
    porosity = 0.70 if screen else 0.38
    bulk_area = pore_area / porosity
    if screen:
        hydraulic_diameter = porosity * wire / (1 - porosity)
        permeability = porosity * hydraulic_diameter**2 / (32 * 2)
        surface_per_length = 4 * (1 - porosity) * bulk_area / wire
        h = 2 * conductivity / hydraulic_diameter
    else:
        permeability = porosity**3 * diameter**2 / (150 * (1 - porosity)**2)
        surface_per_length = 6 * (1 - porosity) * bulk_area / diameter
        h = 2 * conductivity / diameter
    r = prop["viscosity"] / (prop["density"] * bulk_area * permeability)
    conductance = h * surface_per_length
    gas_capacity = prop["density"] * pore_area * prop["cp"]
    matrix_capacity = layer.mass / length * float(layer.material.specific_heat(temperature))
    return prop, pore_area, r, conductance, gas_capacity, matrix_capacity


def coefficients(prop, area, gas_capacity, matrix_capacity, conductance, gradient):
    """Eliminates gas and matrix temperature amplitudes: theta_g = a*p + b*m."""
    omega = 2 * np.pi
    denominator = 1j * omega * gas_capacity
    if matrix_capacity:
        denominator += conductance * 1j * omega * matrix_capacity / (conductance + 1j * omega * matrix_capacity)
    a = 1j * omega * prop["expansion_coefficient"] * area / denominator
    b = -prop["cp"] * gradient / denominator
    return a, b, denominator


def cell(prop, area, a, b, resistance, length, inertance):
    """Transfers [p,m] from warm to cold for a regenerator cell."""
    omega = 2 * np.pi
    rho_p = prop["density"] * prop["isothermal_compressibility"]
    rho_t = -prop["density"] * prop["expansion_coefficient"] / prop["temperature"]
    matrix = np.array([[0, -resistance - 1j * omega * inertance],
                       [-1j * omega * area * (rho_p + rho_t * a), -1j * omega * area * rho_t * b]])
    return expm(matrix * length)


def profiles(model, temperature):
    sections = []
    offset = 1
    for section in model.regenerators:
        centers, layers, position = [], [], 0.0
        for layer in section.layers:
            length = section.length * layer.length_fraction
            centers.extend(position + (np.arange(layer.node_count) + 0.5) * length / layer.node_count)
            layers.extend([layer] * layer.node_count)
            position += length
        sections.append((np.r_[0, centers, section.length], temperature[offset - 1:offset + len(centers) + 1]))
        offset += len(centers) + 1
    return sections


def network(model, temperature, diameter=0.0003, wire=0.00004, physical=True):
    """Solves the complex linear network; all mean profiles and drive are prescribed."""
    omega, mean_pressure, amplitude = 2 * np.pi, 1565067.2928375078, 1e6
    transfers, endpoints, samples, section_cells = [], [], [], []
    for section, (x, ts) in zip(model.regenerators, profiles(model, temperature), strict=True):
        total = np.eye(2, dtype=complex)
        position = 0.0
        boundary_coefficients = []
        cells = []
        for layer in section.layers:
            length = section.length * layer.length_fraction
            for j in range(12):
                left = position + j * length / 12
                right = left + length / 12
                mid = (left + right) / 2
                temp = np.interp(mid, x, ts)
                interval = min(np.searchsorted(x, mid), len(x) - 1)
                gradient = (ts[interval] - ts[interval - 1]) / (x[interval] - x[interval - 1])
                prop, area, r, h, cg, cm = closure(layer, length, temp, diameter, wire)
                a, b, denominator = coefficients(prop, area, cg, cm, h if physical else 1e15, gradient)
                matrix = cell(prop, area, a, b, r if physical else 0, right - left, 1 / area)
                total = matrix @ total
                porosity = .70 if layer.material.name == "304 SS" else .38
                cells.append((matrix, prop, area / porosity, a, b, porosity))
                samples.append((r, h, cg, cm, float(temp), layer.material.name))
            position += length
        for end, layer in [(0, section.layers[0]), (-1, section.layers[-1])]:
            length = section.length * layer.length_fraction
            prop, area, r, h, cg, cm = closure(layer, length, ts[end], diameter, wire)
            neighbor = 1 if end == 0 else -2
            gradient = (ts[end] - ts[neighbor]) / (x[end] - x[neighbor])
            a, b, _ = coefficients(prop, area, cg, cm, h if physical else 1e15, gradient)
            boundary_coefficients.append((prop, a, b, gradient))
        transfers.append(total)
        endpoints.append(boundary_coefficients)
        section_cells.append(cells)

    tube_admittances = []
    for tube, cold, resistance in zip(model.pulse_tubes, [temperature[3], temperature[-1]], [1e10, 3e10], strict=True):
        total = np.eye(2, dtype=complex)
        area = tube.volume / tube.length
        for j in range(24):
            temp = cold + (j + .5) / 24 * (295 - cold)
            prop = properties_TP(float(temp), mean_pressure)
            a, b, _ = coefficients(prop, area, prop["density"] * area * prop["cp"], 0, 0, (295 - cold) / tube.length)
            total = cell(prop, area, a, b, 0, tube.length / 24, 0) @ total
        prop = properties_TP(295, mean_pressure)
        compliance = tube.reservoir_volume / (prop["gamma"] * mean_pressure)
        warm_y = prop["density"] / (resistance + 1 / (1j * omega * compliance))
        tube_admittances.append((warm_y * total[0, 0] - total[1, 0]) / (total[1, 1] - warm_y * total[0, 1]))

    def inlet_y(transfer, outlet_y):
        return (outlet_y * transfer[0, 0] - transfer[1, 0]) / (transfer[1, 1] - outlet_y * transfer[0, 1])

    y2 = inlet_y(transfers[1], tube_admittances[1])
    y1 = inlet_y(transfers[0], tube_admittances[0] + y2)
    inlet = np.array([amplitude, amplitude * y1])
    junction = transfers[0] @ inlet
    second_in = np.array([junction[0], junction[0] * y2])
    cold = transfers[1] @ second_in
    sanity = []
    for state, cells in zip([inlet, second_in], section_cells, strict=True):
        theta_ratios, reynolds, inertial_ratios = [], [], []
        for matrix, prop, bulk_area, a, b, porosity in cells:
            theta_ratios.append(float(abs(a * state[0] + b * state[1]) / prop["temperature"]))
            if porosity == .38:
                re = abs(state[1]) * diameter / (bulk_area * prop["viscosity"])
                reynolds.append(float(re))
                inertial_ratios.append(float(8 / (3*np.pi) * 1.75 / 150 / (1-porosity) * re))
            state = matrix @ state
        sanity.append({"max_theta_g_over_mean_T": max(theta_ratios),
                       "max_sphere_Re_superficial": max(reynolds),
                       "max_Forchheimer_fundamental_over_Darcy": max(inertial_ratios)})

    def flux(state, coeff):
        prop, a, b, gradient = coeff
        p, m = state
        theta = a * p + b * m
        enthalpy = prop["cp"] * theta + (1 - prop["expansion_coefficient"]) * p / prop["density"]
        return float(.5 * np.real(enthalpy * np.conj(m)))

    f1 = flux(junction, endpoints[0][1])
    f2warm = flux(second_in, endpoints[1][0])
    f2cold = flux(cold, endpoints[1][1])
    tube_power = []
    for state, y, temp in zip([junction, cold], tube_admittances, [temperature[3], temperature[-1]], strict=True):
        rho = properties_TP(float(temp), mean_pressure)["density"]
        tube_power.append(float(.5 * abs(state[0])**2 * np.real(y) / rho))
    return {"Q_transport_W": [f1 - f2warm - tube_power[0], f2cold - tube_power[1]],
            "stage_pressure_over_drive": [float(abs(junction[0]) / amplitude), float(abs(cold[0]) / amplitude)],
            "stage_mass_amplitude_kg_s": [float(abs(junction[1])), float(abs(cold[1]))],
            "amplitude_sanity": sanity,
            "inlet_mass_amplitude_kg_s": float(abs(inlet[1])), "samples": samples}


def main():
    start = time.perf_counter()
    model = load_head_model("PT1.5")
    case_file = Path(__file__).parent.parent / "coupling_cases.json"
    temperature = np.array(json.loads(case_file.read_text())[0]["temperature"])
    rows = []
    for layer, length, temp in [(model.regenerators[0].layers[0], .18*.72, 228.2),
                                (model.regenerators[0].layers[1], .18*.28, 96.8),
                                (model.regenerators[1].layers[0], .2445*.16, 50.0),
                                (model.regenerators[1].layers[1], .2445*.32, 31.6),
                                (model.regenerators[1].layers[2], .2445*.32, 14.8),
                                (model.regenerators[1].layers[3], .2445*.20, 7.5),
                                (model.regenerators[1].layers[3], .2445*.20, 4.0),
                                (model.regenerators[1].layers[3], .2445*.20, 3.4)]:
        prop, area, r, h, cg, cm = closure(layer, length, temp, .0003, .00004)
        _, _, d = coefficients(prop, area, cg, cm, h, -100)
        re = .01 * .0003 / ((area/.38) * prop["viscosity"])
        pr = prop["cp"] * prop["viscosity"] / max(prop["thermal_conductivity"], .0243652447)
        nu = 2 + 1.1 * re**.6 * pr**(1/3)
        rows.append({"material": layer.material.name, "temperature_K": temp, "r_Pa_s_per_kg_m": r,
                     "H_W_per_m_K": h, "omega_Cg_over_H": 2*np.pi*cg/h,
                     "omega_Cm_over_H": 2*np.pi*cm/h,
                     "sphere_H_Wakao_at_0.01_kg_s": h * nu/2 if layer.material.name != "304 SS" else None,
                     "sphere_r_eff_at_0.01_kg_s": r * (1 + 8/(3*np.pi)*1.75/150/.62*re) if layer.material.name != "304 SS" else None,
                     "a_equivalent_sinusoidal": float(.5 * cm * np.real(1/d))})
    results = []
    for physical, diameter, wire in [(False, .0003, .00004), (True, .0003, .00004), (True, .0002, .00003), (True, .0005, .00006)]:
        base = network(model, temperature, diameter, wire, physical)
        derivatives = {}
        for path in ["fixed_matrix", "smooth_profile"]:
            derivative = []
            for index in [3, len(temperature)-1]:
                direction = np.zeros_like(temperature)
                direction[index] = 1
                if path == "smooth_profile":
                    for section_index, (x, _) in enumerate(profiles(model, temperature)):
                        offset = 0 if section_index == 0 else 3
                        if index == 3:
                            direction[offset:offset+len(x)] = x / x[-1] if section_index == 0 else 1 - x / x[-1]
                        elif section_index == 1:
                            direction[offset:offset+len(x)] = x / x[-1]
                upper, lower = temperature + .005 * direction, temperature - .005 * direction
                upper_q = np.array(network(model, upper, diameter, wire, physical)["Q_transport_W"])
                lower_q = np.array(network(model, lower, diameter, wire, physical)["Q_transport_W"])
                derivative.append((upper_q - lower_q) / .01)
            derivatives[path] = np.array(derivative).T.tolist()
        base.pop("samples")
        results.append({"physical_r_H": physical, "sphere_m": diameter, "wire_m": wire,
                        **base, "unrelaxed_derivative_W_K": derivatives})
    output = {"scales": rows, "network": results, "runtime_s": time.perf_counter()-start}
    path = Path(__file__).with_name("results.json")
    path.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
