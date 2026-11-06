"""Replays one saved equilibrium after correcting ideal-tube transport."""

import json
import time
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

from pulse_tube_cooldown.library.heads import load_head_model
from simulations.validation.test_stand_steady_state.phase_feasibility.balanced import (
    BudgetExpiredError,
    serializable,
    solve_network,
)


def main():
    """Checks one corrected equilibrium within a fixed runtime budget."""
    start = time.perf_counter()
    directory = Path(__file__).parent
    saved = json.loads((directory / "balanced_results.json").read_text())
    model = load_head_model("PT1.5")
    temperature = np.array(saved["temperature_K"])
    internal = np.r_[1:3, 4:15]
    initial = solve_network(model, temperature)
    best = {"norm": np.inf, "temperature": temperature, "result": initial, "evaluations": 0}

    def residual(values):
        if time.perf_counter() - start > 55:
            raise BudgetExpiredError
        trial = temperature.copy()
        trial[internal] = values
        answer = solve_network(model, trial)
        best["evaluations"] += 1
        norm = float(np.max(abs(answer["residual"])))
        if norm < best["norm"] and answer["harmonic_converged"]:
            best.update(norm=norm, temperature=trial, result=answer)
        return answer["residual"]

    try:
        solution = least_squares(
            residual,
            temperature[internal],
            bounds=(2.5, 295),
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
        "status": status,
        "fixed_drive_peak_Pa": 1e6,
        "corrected_at_saved_profile": serializable(initial),
        "best": serializable(best["result"]),
        "temperature_K": best["temperature"].tolist(),
        "max_matrix_residual_W": best["norm"],
        "mean_evaluations": best["evaluations"],
        "runtime_s": time.perf_counter() - start,
    }
    (directory / "audit_results.json").write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
