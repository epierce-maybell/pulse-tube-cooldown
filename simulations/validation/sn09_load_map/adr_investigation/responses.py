"""Measure the May map's local response and decompose the fitted total cooling."""

import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import fit


def main():
    inputs, data, warm, drive = fit.read_experiment()
    head = fit.read_head(inputs["head"])
    walls = fit.build_walls()
    saved = json.loads((fit.HERE / "results.json").read_text())
    p = fit.Parameters(
        packing_constant=tuple(saved["packing_constants"]),
        conductance=tuple(saved["conductances_W_per_Pa2"]),
        exchanger_effectiveness=saved["exchanger_effectiveness"], frequency=saved["frequency_Hz"],
    )
    plane = np.column_stack([np.ones(6), data.T1.values, data.T2.values])
    observations = np.column_stack([data.Q1.values, data.Q2.values])
    coefficients = np.linalg.lstsq(plane, observations, rcond=None)[0]

    def quantities(stage, guess):
        flows, column = fit.calc_stage_heat_flow(
            head, p, drive, walls, warm_temperature=warm, stage_temperature=tuple(stage), guess=guess,
        )
        t = column.temperature
        mean = (column.supply_pressure + column.return_pressure) / 2
        swing = column.supply_pressure - column.return_pressure
        gross = sum(fit.calc_enthalpy_flow(c, swing, drive.waveform_factor) for c in p.conductance)
        wall = fit.calc_link_heat_flow(walls, {fit.WARM_FLANGE: warm, fit.STAGE_1: stage[0], fit.STAGE_2: stage[1]})
        wall_sum = wall[fit.STAGE_1] + wall[fit.STAGE_2]
        mass = fit.calc_mass_flow(t, head.pore_volume, head.pulse_tube_volume, column.supply_pressure,
                                  column.return_pressure, p.frequency)
        k = head.stage_1 - 1
        cp = fit.helium.isobaric_specific_heat((t[k] + t[k+1])/2, mean)
        capacity = head.calc_heat_capacity(t)[0][-1]
        first_loss = 2*p.packing_constant[0]*(mass[k]*cp)**2/(p.frequency*capacity)*(t[k]-t[k+1])
        np.testing.assert_allclose(-flows.sum(), gross-first_loss-wall_sum, atol=1e-7)
        return np.r_[-flows, gross, first_loss, wall_sum, swing]

    derivatives = []
    for index in (1, 4):
        row = saved["points"][index]
        stage = np.array(row["stage_temperature_K"])
        guess = np.array(row["temperature_K"])
        estimates = []
        for step in (0.002, 0.001):
            jac = []
            for j in (0, 1):
                delta = np.eye(2)[j]*step
                jac.append(((quantities(stage+delta, guess)-quantities(stage-delta, guess))/(2*step)).tolist())
            estimates.append(dict(step_K=step, d_dT1=jac[0], d_dT2=jac[1]))
        derivatives.append(dict(point=index, heater_W=row["heater_power_W"], estimates=estimates))
    output = dict(
        empirical_plane_rows=["constant", "T1", "T2"],
        empirical_plane_columns=["Q1", "Q2"],
        empirical_plane_coefficients=coefficients.tolist(),
        empirical_plane_rms_W=np.sqrt(np.mean((plane@coefficients-observations)**2, axis=0)).tolist(),
        model_derivative_order=["Q1", "Q2", "H1+H2", "L1", "W1+W2", "pressure_swing"],
        derivatives=derivatives,
        interpretation="Empirical planes summarize six operating points; model derivatives hold the other plate fixed. Neither is a measured differential response.",
    )
    (HERE / "responses.json").write_text(json.dumps(output, indent=2)+"\n")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
