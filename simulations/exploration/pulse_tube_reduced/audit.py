"""Compare the measured stage-2 cooldown budget with a nearby loaded point."""

import json
import tomllib
import warnings
from pathlib import Path

import numpy as np
import xarray as xr

from pulse_tube_cooldown.library.heads.load import load_head_model
from simulations.exploration.pulse_tube_candidate.run import Stand
from simulations.validation.head_comparison.comparison import load_population


def main() -> None:
    """Records local measured rates and the consequences of stage-lumped storage."""
    here = Path(__file__).resolve().parent
    repository = here.parents[2]
    config = tomllib.loads((here / "settings.toml").read_text())
    audit = tomllib.loads((here / "five_coefficient.toml").read_text())["audit"]
    head = next(
        h for h in load_population(repository / config["records"], root=repository) if h.serial == config["serial"]
    )
    point = head.load_map.data.isel(point=audit["map_point"])
    stand = Stand.load(repository / config["stand"])
    design = load_head_model(config["head_design"])
    data = head.cooldown.data
    nearest = int(np.nanargmin(abs(data.T2.values - point.T2.item())))
    center_time = data.time[nearest].item()
    rows = []
    for half_window in audit["rate_half_windows"]:
        selected = data.where(abs(data.time - center_time) <= half_window, drop=True).dropna("time")
        coordinates = selected.time.values - center_time
        # A local quadratic resolves the curvature across the measurement window.
        fits = np.array([np.polynomial.polynomial.polyfit(coordinates, selected[name], 2) for name in ("T1", "T2")])
        temperature, rate = fits[:, 0], fits[:, 1]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            capacity = stand.stages[1].calc_heat_capacity(float(temperature[1]))
            capacity += sum(
                float(layer.material.specific_heat(temperature[1])) * layer.mass
                for layer in design.regenerators[1].layers
            )
        rows.append(
            {
                "half_window_s": half_window,
                "samples": selected.sizes["time"],
                "time_s": center_time,
                "temperature_K": temperature.tolist(),
                "rate_K_per_s": rate.tolist(),
                "stage_2_lumped_capacity_J_per_K": capacity,
                "implied_net_cooling_W": float(-capacity * rate[1]),
                "capacity_for_map_power_J_per_K": float(-point.Q2.item() / rate[1]),
                "map_power_to_implied_power_ratio": float(point.Q2.item() / (-capacity * rate[1])),
            }
        )
    destination = here / "results/five_coefficient/lumped_regenerators"
    with xr.open_dataset(destination / "map_calibration/load_map.nc") as prediction:
        map_rmse = {
            name: float(np.sqrt(np.mean((prediction[name].values - head.load_map.data[name].values) ** 2)))
            for name in ("T1", "T2")
        }
    report = {
        "map_point": {name: point[name].item() for name in ("Q1", "Q2", "T1", "T2")},
        "cooldown_budget": rows,
        "five_coefficient_map_calibration_rmse_K": map_rmse,
        "assumptions": [
            "The test-stand stage-2 external heater and parasitic heat inputs during cooldown are zero.",
            "All stage-2 regenerator material shares the stage-2 temperature in this storage approximation.",
            "Comparing cooling powers assumes unchanged head operation between the selected runs.",
            "Positive heat-flow derivative dQ2/dT1 would make cooling stronger at the lower cooldown T1.",
        ],
        "recorded_run_differences": {"cooldown": head.cooldown.note, "load_map": head.load_map.note},
    }
    (destination / "budget_audit.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
