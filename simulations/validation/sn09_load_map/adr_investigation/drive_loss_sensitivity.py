"""Separate local drive and first-regenerator effects at the May middle loads.

Each ablation agrees with the saved equilibrium at its reference point. Linearized
matrix balances then give the response to T2 with T1 fixed. Clamped quantities
remain at their reference values; the other quantities respond self-consistently.
"""

import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import fit
from pulse_tube_cooldown.library.links.cooler.head import calc_acoustic_power


def main():
    inputs, _, warm, drive = fit.read_experiment()
    head = fit.read_head(inputs["head"])
    nodes = fit.column_indices(head)
    walls = fit.build_walls()
    saved = json.loads((fit.HERE / "results.json").read_text())
    parameters = fit.Parameters(
        packing_constant=tuple(saved["packing_constants"]),
        conductance=tuple(saved["conductances_W_per_Pa2"]),
        exchanger_effectiveness=saved["exchanger_effectiveness"], frequency=saved["frequency_Hz"],
    )
    rows = []
    for point in (1, 4):
        reference = saved["points"][point]
        template = np.array(reference["temperature_K"])
        mean0 = (reference["supply_pressure_Pa"] + reference["return_pressure_Pa"])/2
        swing0 = reference["supply_pressure_Pa"] - reference["return_pressure_Pa"]
        mass0 = fit.calc_mass_flow(template, head.pore_volume, head.pulse_tube_volume,
                                   mean0+swing0/2, mean0-swing0/2, parameters.frequency)
        gross0 = tuple(fit.calc_enthalpy_flow(c, swing0, drive.waveform_factor) for c in parameters.conductance)
        ratio0 = fit.calc_capacity_ratio(template, head.calc_heat_capacity(template), head.pore_volume, mean0)
        acoustic0 = calc_acoustic_power(template, ratio0, gross0[1], mean0)[0]
        in_phase0 = float(fit.helium.density(template[head.stage_1], mean0))*2*(gross0[0]+acoustic0)/(swing0/2)
        state0 = np.r_[np.log(template[nodes]), mean0/1e6, swing0/1e6]
        for mode in ("full", "fixed_swing", "fixed_mean", "fixed_pressures", "fixed_gross_power", "fixed_first_column_flow",
                     "first_flow_tracks_in_phase_component"):

            def evaluate(state, plate2):
                temperature = template.copy()
                temperature[nodes] = np.exp(state[:-2])
                temperature[-1] = plate2
                mean, swing = state[-2:]*1e6
                flow = fit.calc_mass_flow(temperature, head.pore_volume, head.pulse_tube_volume,
                                          mean+swing/2, mean-swing/2, parameters.frequency)
                thermal_flow = flow.copy()
                if mode == "fixed_first_column_flow":
                    thermal_flow[:head.stage_1] = mass0[:head.stage_1]
                capacity = head.calc_heat_capacity(temperature)
                ratio = fit.calc_capacity_ratio(temperature, capacity, head.pore_volume, mean)
                gross = (gross0 if mode == "fixed_gross_power" else
                         tuple(fit.calc_enthalpy_flow(c, swing, drive.waveform_factor) for c in parameters.conductance))
                acoustic = calc_acoustic_power(temperature, ratio, gross[1], mean)[0]
                in_phase = float(fit.helium.density(temperature[head.stage_1], mean))*2*(gross[0]+acoustic)/(swing/2)
                if mode == "first_flow_tracks_in_phase_component":
                    thermal_flow[:head.stage_1] = mass0[:head.stage_1]*in_phase/in_phase0
                heat = fit.calc_heat_flow(temperature, capacity, ratio, thermal_flow, gross,
                                         parameters.packing_constant, parameters.exchanger_effectiveness,
                                         mean, parameters.frequency)
                inventory = fit.calc_gas_mass(temperature, head.pore_volume, head.pulse_tube_volume, mean)
                inventory += drive.warm_volume*float(fit.helium.density(warm, mean))
                pressure_residual = (inventory-drive.charge)/drive.charge
                if mode in ("fixed_mean", "fixed_pressures"):
                    pressure_residual = (mean-mean0)/1e6
                compressor = float(fit.helium.density(warm, mean-swing/2))*sum(drive.displacement_rate)
                drive_residual = min((min(drive.relief_swing)-swing)/min(drive.relief_swing),
                                     (compressor-flow[0])/0.005)
                if mode in ("fixed_swing", "fixed_pressures"):
                    drive_residual = (swing-swing0)/1e6
                residual = np.r_[heat[nodes]/100, pressure_residual, drive_residual]
                wall = fit.calc_link_heat_flow(walls, {fit.WARM_FLANGE: warm,
                    fit.STAGE_1: temperature[head.stage_1], fit.STAGE_2: plate2})
                cooling = -heat[[head.stage_1, head.stage_2]] - [wall[fit.STAGE_1], wall[fit.STAGE_2]]
                k = head.stage_1-1
                cp = fit.helium.isobaric_specific_heat((temperature[k]+temperature[k+1])/2, mean)
                loss = 2*parameters.packing_constant[0]*(thermal_flow[k]*cp)**2 / (parameters.frequency*capacity[0][-1])
                loss *= temperature[k]-temperature[k+1]
                quantities = np.r_[cooling, sum(gross), loss, wall[fit.STAGE_1]+wall[fit.STAGE_2],
                                    mean, swing, flow[0], thermal_flow[k], acoustic, in_phase]
                return residual, quantities

            base_residual, base_quantities = evaluate(state0, template[-1])
            assert max(abs(base_residual)) < 1e-7
            estimates = []
            for step in (1e-5, 5e-6):
                jacobian, response = [], []
                for direction in np.eye(len(state0)):
                    plus = evaluate(state0+step*direction, template[-1])
                    minus = evaluate(state0-step*direction, template[-1])
                    jacobian.append((plus[0]-minus[0])/(2*step))
                    response.append((plus[1]-minus[1])/(2*step))
                jacobian = np.array(jacobian).T
                response = np.array(response).T
                plate_step = 0.001
                plus, minus = evaluate(state0, template[-1]+plate_step), evaluate(state0, template[-1]-plate_step)
                forcing = (plus[0]-minus[0])/(2*plate_step)
                explicit = (plus[1]-minus[1])/(2*plate_step)
                state_derivative = np.linalg.solve(jacobian, -forcing)
                derivative = response@state_derivative+explicit
                np.testing.assert_allclose(derivative[0]+derivative[1], derivative[2]-derivative[3]-derivative[4], atol=1e-5)
                estimates.append(dict(state_step=step, d_dT2=derivative.tolist(),
                                      condition_number=float(np.linalg.cond(jacobian))))
            row = dict(point=point, heater_W=reference["heater_power_W"], mode=mode,
                       reference=base_quantities.tolist(), estimates=estimates)
            rows.append(row)
            values = estimates[-1]["d_dT2"]
            print(point, mode, "dQsum/dT2", values[0]+values[1], "dH", values[2], "dL1", values[3], flush=True)
    output = dict(quantities=["Q1_W", "Q2_W", "H1_plus_H2_W", "L1_W", "walls_W", "mean_Pa", "swing_Pa",
                              "warm_mass_flow_kg_s", "first_column_cold_mass_flow_kg_s", "stage2_warm_acoustic_W",
                              "first_column_cold_in_phase_mass_flow_kg_s"], rows=rows,
                  in_phase_ablation="First-column loss flows track the required in-phase cold-end flow, normalized to retain each reference equilibrium. The boundary shape of the flow is held fixed. This is a local sensitivity experiment, not a replacement flow solution or refit.")
    (HERE/"drive_loss_sensitivity.json").write_text(json.dumps(output, indent=2)+"\n")


if __name__ == "__main__":
    main()
