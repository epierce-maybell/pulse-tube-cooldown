"""Decomposes the saved joint fit's stage balances without fitting parameters."""

import json

import matplotlib

matplotlib.use("Agg")
import fit
import matplotlib.pyplot as plt
import numpy as np

from pulse_tube_cooldown.library.gases import helium
from pulse_tube_cooldown.library.links.cooler.head import calc_enthalpy_flow, calc_mass_flow


def decompose(head, parameters, drive, temperature, mean, swing):
    """Returns the ADR's boundary flows and stage cooling powers, in W."""
    mass = calc_mass_flow(
        temperature,
        head.pore_volume,
        head.pulse_tube_volume,
        mean + swing / 2,
        mean - swing / 2,
        parameters.frequency,
    )
    capacity = np.r_[0, head.calc_heat_capacity(temperature)[0], 0, head.calc_heat_capacity(temperature)[1], 0]
    midpoint = (temperature[:-1] + temperature[1:]) / 2
    specific_heat = helium.isobaric_specific_heat(midpoint, mean)
    packing = np.repeat(parameters.packing_constant, [head.stage_1, head.stage_2 - head.stage_1])
    conductance = 2 * packing * (mass * specific_heat) ** 2 / (parameters.frequency * (capacity[:-1] + capacity[1:]))
    passive = conductance * (temperature[:-1] - temperature[1:])
    h1, h2 = [calc_enthalpy_flow(c, swing, drive.waveform_factor) for c in parameters.conductance]
    boundary_temperature = midpoint[head.stage_1 :].copy()
    boundary_temperature[[0, -1]] = temperature[[head.stage_1, -1]]
    beta_t = helium.isobaric_expansion_coefficient(boundary_temperature, mean) * boundary_temperature
    density = helium.density(boundary_temperature, mean)
    real_gas = (1 - beta_t) * density[-1] / density * h2
    total = passive[head.stage_1 :] + real_gas
    matrix_temperatures = temperature[head.stage_1 + 1 : -1]
    pore_capacity = (
        head.pore_volume[1]
        * helium.density(matrix_temperatures, mean)
        * helium.isobaric_specific_heat(matrix_temperatures, mean)
    )
    cold_pore_capacity = (
        head.pore_volume[1][-1] * density[-1] * float(helium.isobaric_specific_heat(temperature[-1], mean))
    )
    cold_matrix_capacity = head.mass[1][-1] * float(head.material[1][-1].specific_heat(np.array(temperature[-1])))
    cold_capacity_ratio = cold_pore_capacity / cold_matrix_capacity
    wall = fit.calc_link_heat_flow(
        fit.build_walls(),
        {
            fit.WARM_FLANGE: temperature[0],
            fit.STAGE_1: temperature[head.stage_1],
            fit.STAGE_2: temperature[-1],
        },
    )
    return dict(
        T1_K=float(temperature[head.stage_1]),
        T2_K=float(temperature[-1]),
        H1_W=h1,
        H2_W=h2,
        beta_T2=float(beta_t[-1]),
        corrected_gross_W=float(beta_t[-1] * h2),
        stage1_loss_W=float(passive[head.stage_1 - 1]),
        stage2_passive_boundary_W=passive[head.stage_1 :].tolist(),
        stage2_real_gas_boundary_W=real_gas.tolist(),
        stage2_total_boundary_W=total.tolist(),
        stage2_pore_to_matrix_heat_capacity=(pore_capacity / head.calc_heat_capacity(temperature)[1]).tolist(),
        cold_endpoint_pore_to_matrix_heat_capacity=cold_capacity_ratio,
        omitted_pressure_capacity_power_scale_W=float(
            cold_capacity_ratio / (1 + cold_capacity_ratio) * beta_t[-1] * h2
        ),
        cold_mass_flow_kg_per_s=float(mass[-1]),
        cold_specific_heat_J_per_kg_K=float(specific_heat[-1]),
        cold_matrix_capacity_J_per_K=float(capacity[-2]),
        cold_boundary_conductance_W_per_K=float(conductance[-1]),
        wall1_W=wall[fit.STAGE_1],
        wall2_W=wall[fit.STAGE_2],
        Q1_W=float(h1 - passive[head.stage_1 - 1] + total[0] - wall[fit.STAGE_1]),
        Q2_W=float(h2 - total[-1] - wall[fit.STAGE_2]),
    )


def main():
    """Checks the saved stage powers and writes their component flows and plot."""
    saved = json.loads((fit.HERE / "results.json").read_text())
    inputs, _, _, drive = fit.read_experiment()
    head = fit.read_head(inputs["head"])
    parameters = fit.Parameters(
        packing_constant=tuple(saved["packing_constants"]),
        conductance=tuple(saved["conductances_W_per_Pa2"]),
        frequency=saved["frequency_Hz"],
    )
    rows = []
    for point in saved["points"]:
        row = decompose(
            head,
            parameters,
            drive,
            np.array(point["temperature_K"]),
            (point["supply_pressure_Pa"] + point["return_pressure_Pa"]) / 2,
            point["supply_pressure_Pa"] - point["return_pressure_Pa"],
        )
        row["heater_W"] = point["heater_power_W"]
        np.testing.assert_allclose(
            [row["Q1_W"], row["Q2_W"]],
            np.array(point["heater_power_W"]) - point["residual_W"],
            atol=1e-8,
        )
        rows.append(row)
    (fit.HERE / "decomposition.json").write_text(json.dumps(rows, indent=2) + "\n")
    print("T1/K T2/K H2/W beta*T2 beta*T2*H2/W passive2/W total2/W wall2/W Q2/W")
    for row in rows:
        values = [row[key] for key in ("T1_K", "T2_K", "H2_W", "beta_T2", "corrected_gross_W")]
        values += [
            row["stage2_passive_boundary_W"][-1],
            row["stage2_total_boundary_W"][-1],
            row["wall2_W"],
            row["Q2_W"],
        ]
        print(" ".join(f"{value:.6g}" for value in values))
    for index in (1, 4):
        print(f"Point {index} pore/matrix capacity ratios: {rows[index]['stage2_pore_to_matrix_heat_capacity']}")
    print("Cold-end pressure/capacity term scale, W:", [r["omitted_pressure_capacity_power_scale_W"] for r in rows])
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), layout="constrained")
    for indices, label in (([2, 1, 0], "Q1 heater = 20 W"), ([5, 4, 3], "Q1 heater = 0 W")):
        group = [rows[i] for i in indices]
        t = [r["T2_K"] for r in group]
        axes[0].plot(t, [r["corrected_gross_W"] for r in group], "o-", label=label + ": beta*T2*H2")
        axes[0].plot(t, [r["stage2_passive_boundary_W"][-1] for r in group], "x--", label=label + ": passive loss")
        axes[1].plot(t, [r["stage2_total_boundary_W"][-1] for r in group], "o-", label=label)
        axes[2].plot(t, [r["Q2_W"] for r in group], "x-", label=label + ": model")
        axes[2].plot(t, [r["heater_W"][1] for r in group], "o:", label=label + ": measured")
    for ax, title in zip(
        axes, ("Stage-2 cooling terms", "Net column flow from stage 1 to 2", "Stage-2 available cooling"), strict=True
    ):
        ax.set(xlabel="Measured T2 (K)", ylabel="Power (W)", title=title)
        ax.legend(fontsize=8)
        ax.grid(alpha=0.2)
    fig.savefig(fit.HERE / "decomposition.png", dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    main()
