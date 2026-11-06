"""One compressor-supply check at the audit's fixed mean pressure."""

import json
import time
import tomllib
from pathlib import Path

import numpy as np
from hepak.properties import properties_TP
from scipy.optimize import least_squares

from pulse_tube_cooldown.library.heads import load_head_model
from simulations.validation.test_stand_steady_state.phase_feasibility.balanced import (
    BudgetExpiredError,
    serializable,
    solve_network,
)


def main():
    """Balances matrix heat and average positive inlet mass flow together."""
    start = time.perf_counter()
    directory = Path(__file__).parent
    saved = json.loads((directory / "audit_results.json").read_text())
    temperature = np.array(saved["temperature_K"])
    drive = tomllib.loads(Path("data/test_stand/pulse_tube/drive_one_compressor.toml").read_text())
    compressor = drive["compressors"][0]
    model = load_head_model("PT1.5")
    internal = np.r_[1:3, 4:15]
    mean_pressure = 1565067.2928375078
    best = {"norm": np.inf, "evaluations": 0}

    def residual(values):
        if time.perf_counter() - start > 55:
            raise BudgetExpiredError
        trial = temperature.copy()
        trial[internal] = values[:-1]
        amplitude = values[-1] * 1e6
        answer = solve_network(model, trial, drive_amplitude=amplitude)
        supply = properties_TP(295, mean_pressure - amplitude)["density"] * compressor["displacement_rate"]
        demand = answer["inlet_mass_amplitude_kg_s"] / np.pi
        vector = np.r_[answer["residual"], 1000 * (demand - supply)]
        best["evaluations"] += 1
        norm = float(np.max(abs(vector)))
        if norm < best["norm"] and answer["harmonic_converged"]:
            best.update(
                norm=norm,
                temperature_K=trial.tolist(),
                drive_peak_Pa=amplitude,
                supply_kg_s=supply,
                demand_kg_s=demand,
                result=serializable(answer),
            )
        return vector

    try:
        solution = least_squares(
            residual,
            np.r_[temperature[internal], 0.7],
            bounds=(
                np.r_[np.full(len(internal), 2.5), 0.01],
                np.r_[np.full(len(internal), 295), compressor["relief_swing"] / 2e6],
            ),
            max_nfev=20,
            ftol=1e-8,
            xtol=1e-8,
            gtol=1e-8,
            diff_step=1e-4,
        )
        status = solution.message
    except BudgetExpiredError:
        status = "55 s wall-time limit"
    output = {
        "fixed_mean_pressure_Pa": mean_pressure,
        "status": status,
        "runtime_s": time.perf_counter() - start,
        **best,
    }
    (directory / "supply_audit_results.json").write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
