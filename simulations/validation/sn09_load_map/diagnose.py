"""Decomposes a saved May fit into packing, real-fluid, finite-capacity, and exchanger terms."""

import argparse
import json
from dataclasses import replace

import fit
import matplotlib
import numpy as np

from pulse_tube_cooldown.library.gases import helium
from pulse_tube_cooldown.library.links.cooler.head import calc_enthalpy_flow, calc_mass_flow

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def decompose(head, parameters, drive, temperature, mean, swing):
    """Returns independently reconstructed stage powers and boundary transports in W."""
    mass = calc_mass_flow(
        temperature, head.pore_volume, head.pulse_tube_volume, mean + swing / 2, mean - swing / 2, parameters.frequency
    )
    capacities = head.calc_heat_capacity(temperature)
    capacity = np.r_[0, capacities[0], 0, capacities[1], 0]
    midpoint = (temperature[:-1] + temperature[1:]) / 2
    cp = helium.isobaric_specific_heat(midpoint, mean)
    packing = np.repeat(parameters.packing_constant, [head.stage_1, head.stage_2 - head.stage_1])
    passive = 2 * packing * (mass * cp) ** 2 / (parameters.frequency * (capacity[:-1] + capacity[1:]))
    passive *= temperature[:-1] - temperature[1:]
    h1, h2 = [calc_enthalpy_flow(c, swing, drive.waveform_factor) for c in parameters.conductance]
    eps = fit.calc_capacity_ratio(temperature, capacities, head.pore_volume, mean)[1]
    boundary = midpoint[head.stage_1 :].copy()
    boundary[[0, -1]] = temperature[[head.stage_1, -2 if head.cold_exchanger_coupling else -1]]
    density = helium.density(boundary, mean)
    attenuation = np.log(density[1:] / density[:-1]) / (1 + eps)
    acoustic = h2 * np.exp(np.r_[np.cumsum(attenuation[::-1])[::-1], 0.0])
    boundary_eps = np.r_[eps[0], (eps[:-1] + eps[1:]) / 2, eps[-1]]
    beta_t = helium.isobaric_expansion_coefficient(boundary, mean) * boundary
    real_fluid = (1 - beta_t) * acoustic
    finite_capacity = boundary_eps / (1 + boundary_eps) * beta_t * acoustic
    pressure_flow = real_fluid + finite_capacity
    count = -1 if head.cold_exchanger_coupling else None
    total = passive[head.stage_1 : count] + pressure_flow[:count]
    gross = beta_t[-1] * h2 / (1 + eps[-1])
    loss = passive[-2 if head.cold_exchanger_coupling else -1]
    remainder = pressure_flow[-2 if head.cold_exchanger_coupling else -1] - (h2 - gross)
    exchanger = (
        parameters.exchanger_effectiveness * mass[-1] * float(helium.isobaric_specific_heat(temperature[-1], mean))
    )
    wall = fit.calc_link_heat_flow(
        fit.build_walls(),
        {fit.WARM_FLANGE: temperature[0], fit.STAGE_1: temperature[head.stage_1], fit.STAGE_2: temperature[-1]},
    )
    return dict(
        T1_K=float(temperature[head.stage_1]),
        T2_K=float(temperature[-1]),
        cold_matrix_temperature_K=float(temperature[-2]),
        H1_W=h1,
        H2_W=h2,
        beta_cold_T=float(beta_t[-1]),
        corrected_gross_W=float(gross),
        cold_packing_loss_W=float(loss),
        cold_pressure_span_W=float(remainder),
        stage2_acoustic_boundary_W=acoustic.tolist(),
        stage2_passive_boundary_W=passive[head.stage_1 : count].tolist(),
        stage2_real_fluid_boundary_W=real_fluid[:count].tolist(),
        stage2_finite_capacity_boundary_W=finite_capacity[:count].tolist(),
        stage2_total_boundary_W=total.tolist(),
        stage2_pore_to_matrix_heat_capacity=eps.tolist(),
        exchanger_conductance_W_per_K=exchanger,
        wall1_W=wall[fit.STAGE_1],
        wall2_W=wall[fit.STAGE_2],
        Q1_W=float(h1 - passive[head.stage_1 - 1] + total[0] - wall[fit.STAGE_1]),
        Q2_W=float(h2 - total[-1] - wall[fit.STAGE_2]),
    )


def main(folder=fit.HERE):
    """Checks both stage balances and plots the saved fit's cold-end decomposition."""
    saved = json.loads((folder / "results.json").read_text())
    inputs, _, _, drive = fit.read_experiment()
    head = replace(fit.read_head(inputs["head"]), cold_exchanger_coupling=saved["cold_exchanger_coupling"])
    parameters = fit.Parameters(
        packing_constant=tuple(saved["packing_constants"]),
        conductance=tuple(saved["conductances_W_per_Pa2"]),
        exchanger_effectiveness=saved["exchanger_effectiveness"],
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
            [row["Q1_W"], row["Q2_W"]], np.array(point["heater_power_W"]) - point["residual_W"], atol=1e-7, rtol=0
        )
        np.testing.assert_allclose(
            row["Q2_W"],
            row["corrected_gross_W"] - row["cold_packing_loss_W"] - row["cold_pressure_span_W"] - row["wall2_W"],
            atol=1e-10,
        )
        rows.append(row)
    (folder / "decomposition.json").write_text(json.dumps(rows, indent=2) + "\n")
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), layout="constrained")
    for indices, label in (([2, 1, 0], "20 W stage 1"), ([5, 4, 3], "0 W stage 1")):
        group = [rows[i] for i in indices]
        t = [r["T2_K"] for r in group]
        for key, style in [("corrected_gross_W", "o-"), ("cold_packing_loss_W", "x--"), ("cold_pressure_span_W", "+:")]:
            axes[0].plot(t, [r[key] for r in group], style, label=label + ": " + key)
        axes[1].plot(t, [r["cold_matrix_temperature_K"] for r in group], "o-", label=label)
        axes[2].plot(t, [r["Q2_W"] for r in group], "x-", label=label + ": model")
        axes[2].plot(t, [r["heater_W"][1] for r in group], "o:", label=label + ": measured")
    for ax, title in zip(axes, ("Cold-end power terms", "Matrix minimum", "Available stage-2 cooling"), strict=True):
        ax.set(xlabel="Measured T2 (K)", ylabel="Temperature (K)" if ax is axes[1] else "Power (W)", title=title)
        ax.legend(fontsize=6)
        ax.grid(alpha=0.2)
    fig.savefig(folder / "decomposition.png", dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folder", nargs="?", type=fit.Path, default=fit.HERE)
    main(parser.parse_args().folder)
