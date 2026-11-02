"""First-pass May cooldown and August predictions using a frozen candidate May fit.

Run from the repository root with
`python simulations/validation/sn09_bare_head/run.py`.
Inputs, temperature histories, comparisons, and figures live beside this file.
"""

import hashlib
import json
import subprocess
import time
import tomllib
from collections.abc import Mapping
from pathlib import Path

import numpy as np
import xarray as xr
from numpy.typing import NDArray
from scipy.optimize import least_squares

from pulse_tube_cooldown.core.types import Model, Node, Phase, Procedure
from pulse_tube_cooldown.engine.equations import calc_net_heat_flow
from pulse_tube_cooldown.engine.solver import simulate_cooldown
from pulse_tube_cooldown.library.links.conduction import Conduction
from pulse_tube_cooldown.library.links.cooler import Compressor, PulseTubeCooler, calc_charge, load_pulse_tube_head
from pulse_tube_cooldown.library.links.cooler.drive import calc_mean_pressure, calc_pressure_swing
from pulse_tube_cooldown.library.links.cooler.head import calc_gas_mass, calc_mass_flow
from pulse_tube_cooldown.library.materials.load import load_material
from pulse_tube_cooldown.measurements.test_stand import read_cooldown, read_load_map

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
DATA = ROOT / "data/sn09/pulse_tube"


def read_inputs() -> dict:
    """Reads the experiment's inputs and their sources."""
    return tomllib.loads((HERE / "inputs.toml").read_text(encoding="utf-8"))


def read_measurements(run: str) -> xr.Dataset:
    """Reads one cooldown, with seconds since its first timestamp as the clock."""
    config = read_inputs()["runs"][run]
    data = read_cooldown(DATA / config["cooldown"])
    origin = data.time.values[0]
    if "end_exclusive" in config:
        data = data.isel(time=data.time.values < np.datetime64(config["end_exclusive"]))
    seconds = np.asarray((data.time.values - origin) / np.timedelta64(1, "s"), dtype=np.float64)
    return data.assign_coords(time=seconds).assign_attrs(
        origin=str(origin), end_exclusive=config.get("end_exclusive", "end of file"), time_units="s"
    )


def build_walls(geometry: Path, material: str) -> tuple[Conduction, ...]:
    """Builds the same four axial wall links used for the candidate May fit."""
    rows = np.genfromtxt(geometry, delimiter=",", names=True, dtype=None, encoding="utf-8", usecols=(0, 1))
    values = {row["quantity"]: row["value"] for row in rows}
    thickness = float(values["tube_wall_thickness"]) * 1e-3
    endpoints = {
        "stage1_regenerator": ("Warm flange", "Stage 1"),
        "stage1_pulse_tube": ("Warm flange", "Stage 1"),
        "stage2_regenerator": ("Stage 1", "Stage 2"),
        "stage2_pulse_tube": ("Warm flange", "Stage 2"),
    }
    walls = []
    for tube, nodes in endpoints.items():
        radius = float(values[f"{tube}_outer_radius"]) * 1e-3
        walls.append(
            Conduction(
                name=f"{tube} wall",
                nodes=nodes,
                material=load_material(material),
                area=np.pi * (radius**2 - (radius - thickness) ** 2),
                length=float(values[f"{tube}_length"]) * 1e-3,
            )
        )
    return tuple(walls)


def build_model(run: str, max_time: float | None = None) -> Model:
    """Builds the stand network with the run's compressor count and frozen fit.

    Args:
        run: Run key in inputs.toml, either may or august.
        max_time: Integration duration in s, defaulting to the measurement window.

    Returns:
        The bare head with a fixed warm flange, stage inventories, matrix Nodes,
        and four tube-wall links. Initial temperatures equal the first TPR reading.
    """
    inputs = read_inputs()
    measured = read_measurements(run)
    warm = float(measured.TPR[0])
    head = load_pulse_tube_head(
        name="SN09", stages=("Stage 1", "Stage 2"), model=inputs["head"], parameters=HERE / inputs["parameters"]
    )
    inventory = inputs["inventory"]
    plates = tuple(
        Node(name=name, composition=tuple((load_material(material), mass) for material, mass in composition))
        for name, composition in (("Stage 1", inventory["stage_1"]), ("Stage 2", inventory["stage_2"]))
    )
    # The warm flange is a fixed boundary; its composition supplies the Node's thermal mass.
    warm_node = Node(name="Warm flange", composition=plates[0].composition)
    nodes = (warm_node, *plates, *head.matrix_nodes)
    compressor = Compressor(displacement_rate=inputs["displacement_rate"], relief_swing=inputs["relief_swing"])
    cooler = PulseTubeCooler(
        name="Cooler",
        warm_flange=warm_node.name,
        heads=(head,),
        compressors=(compressor,) * inputs["runs"][run]["compressors"],
        charge=calc_charge((head,), inputs["fill_pressure"], warm, inputs["warm_volume"]),
        warm_volume=inputs["warm_volume"],
        waveform_factor=inputs["waveform_factor"],
    )
    links = (cooler, *build_walls(ROOT / inputs["geometry"], inputs["wall_material"]))
    return Model(
        meta={"experiment": "SN09 bare head, provisional May fit", "run": run},
        nodes=nodes,
        procedure=Procedure(phases=(Phase(name="cooldown", links=links, until=None),), stop_condition=None),
        initial_temperature=np.full(len(nodes), warm),
        fixed_boundary_nodes=(0,),
        solver_options={**inputs["solver"], "max_time": float(measured.time[-1]) if max_time is None else max_time},
    )


def compare_cooldown(measured: xr.Dataset, temperature: xr.DataArray) -> dict:
    """Scores temperatures interpolated to valid measured timestamps in their shared window.

    Each channel reports sample count, sum of squared natural-log residuals,
    RMS log residual, and RMS temperature error in K.
    """
    comparison = {}
    for channel, node in (("T1", "Stage 1"), ("T2", "Stage 2")):
        observed = measured[channel].values
        valid = np.isfinite(observed) & (observed > 0)
        valid &= measured.time.values >= temperature.time.values[0]
        valid &= measured.time.values <= temperature.time.values[-1]
        if channel == "T1":
            valid &= observed < 219.0
        prediction = np.interp(measured.time.values[valid], temperature.time, temperature.sel(node=node))
        log_residual = np.log(prediction / observed[valid])
        comparison[channel] = {
            "samples": int(valid.sum()),
            "sum_squared_log_residual": float(np.sum(log_residual**2)),
            "rms_log_residual": float(np.sqrt(np.mean(log_residual**2))) if valid.any() else None,
            "rms_temperature_error_K": float(np.sqrt(np.mean((prediction - observed[valid]) ** 2)))
            if valid.any()
            else None,
        }
    return comparison


def predict_load_map(model: Model, cooldown: xr.DataArray, load_map: xr.Dataset) -> list[dict]:
    """Predicts available stage cooling at August's measured stage temperatures.

    Solves only matrix-node balances, starting from the cooldown profile nearest
    each measured pair of stage temperatures. Reports stage power residuals,
    the model clamp temperature, and the independent T3 reading (logged as TPR).
    An unconverged point has no power prediction.
    """
    graph = model.graphs[0]
    head = graph.links[0].heads[0]
    stage_indices = [graph.get_node_index(node) for node in head.stages]
    matrix_indices = np.array([graph.get_node_index(node.name) for node in head.matrix_nodes])
    history = cooldown.values
    predictions = []

    def balance(log_temperature: NDArray[np.float64], temperature: NDArray[np.float64]) -> NDArray[np.float64]:
        temperature[matrix_indices] = np.exp(log_temperature)
        return calc_net_heat_flow(graph, temperature)[matrix_indices]

    for point in range(load_map.sizes["point"]):
        observed = np.array([float(load_map[name][point]) for name in ("T1", "T2")])
        heater = np.array([float(load_map[name][point]) for name in ("Q1", "Q2")])
        nearest = np.sum(np.log(history[:, stage_indices] / observed) ** 2, axis=1).argmin()
        temperature = history[nearest].copy()
        temperature[stage_indices] = observed

        solution = least_squares(
            balance,
            np.log(np.clip(temperature[matrix_indices], 2.001, 299.999)),
            args=(temperature,),
            bounds=(np.log(2.0), np.log(300.0)),
            diff_step=1e-5,
            x_scale="jac",
            ftol=1e-10,
            xtol=1e-10,
            gtol=1e-10,
            max_nfev=300,
        )
        balance(solution.x, temperature)
        flow = calc_net_heat_flow(graph, temperature)
        max_balance = float(np.max(np.abs(flow[matrix_indices])))
        converged = bool(solution.success and max_balance < 1e-6)
        predictions.append(
            {
                "point": point,
                "converged": converged,
                "solver_message": solution.message,
                "stage_temperature_K": observed.tolist(),
                "heater_power_W": heater.tolist(),
                "predicted_cooling_W": (-flow[stage_indices]).tolist() if converged else None,
                "residual_W": (flow[stage_indices] + heater).tolist() if converged else None,
                "max_matrix_balance_W": max_balance,
                "clamp_temperature_K": float(temperature[graph.get_node_index(head.clamp_node)]) if converged else None,
                "TPR_K": float(load_map.TPR[point]),
                "temperature_K": temperature.tolist(),
            }
        )
    return predictions


def describe_state(model: Model, temperature: NDArray[np.float64]) -> dict:
    """Reports the drive state and largest free-node temperature rate at a sample."""
    graph = model.graphs[0]
    cooler = graph.links[0]
    head = cooler.heads[0]
    chain = temperature[[graph.get_node_index(name) for name in head.chain_nodes(cooler.warm_flange)]]
    mean = calc_mean_pressure(
        lambda pressure: calc_gas_mass(chain, head.pore_volume, head.pulse_tube_volumes, pressure),
        cooler.charge,
        cooler.warm_volume,
        float(chain[0]),
    )
    swing = calc_pressure_swing(
        lambda supply, back: calc_mass_flow(
            chain, head.pore_volume, head.pulse_tube_volumes, supply, back, head.frequency
        )[0],
        mean,
        float(chain[0]),
        [compressor.displacement_rate for compressor in cooler.compressors],
        [compressor.relief_swing for compressor in cooler.compressors],
    )
    flow = calc_net_heat_flow(graph, temperature)
    free = [index for index in range(len(model.nodes)) if index not in model.fixed_boundary_nodes]
    rates = [flow[index] / model.nodes[index].calc_heat_capacity(float(temperature[index])) for index in free]
    return {
        "mean_pressure_Pa": mean,
        "swing_Pa": swing,
        "relief_cap_Pa": min(compressor.relief_swing for compressor in cooler.compressors),
        "max_free_node_rate_K_per_s": float(np.max(np.abs(rates))),
    }


def plot_cooldowns(histories: Mapping[str, xr.DataArray], destination: Path) -> None:
    """Plots T1 and T2 overlays and the measured pulse-tube T3 diagnostic."""
    import matplotlib.pyplot as plt

    figure, axes = plt.subplots(3, len(histories), figsize=(12, 10), layout="constrained", squeeze=False)
    for column, (run, temperature) in enumerate(histories.items()):
        measured = read_measurements(run)
        for row, (channel, node) in enumerate((("T1", "Stage 1"), ("T2", "Stage 2"), ("TPR", None))):
            axis = axes[row, column]
            valid = np.isfinite(measured[channel].values)
            if channel == "T1":
                valid &= measured[channel].values < 219.0
            label = "Measured T3 (logged TPR)" if row == 2 else f"Measured {channel}"
            axis.plot(measured.time.values[valid] / 60, measured[channel].values[valid], label=label)
            if node is not None:
                axis.plot(temperature.time / 60, temperature.sel(node=node), label="Model")
            axis.set(
                xlabel="Time since first record (min)",
                ylabel="Temperature (K)",
                title=f"{run.title()}: " + ("T3 on pulse tube (unscored)" if row == 2 else channel),
            )
            axis.grid(alpha=0.25)
            axis.legend()
    figure.suptitle("SN09 bare head: provisional May fit, tube walls only")
    figure.savefig(destination, dpi=160)
    plt.close(figure)


def plot_load_map(points: list[dict], destination: Path) -> None:
    """Plots predicted August cooling and its residual at the measured temperatures."""
    import matplotlib.pyplot as plt

    converged = [point for point in points if point["converged"]]
    figure, axes = plt.subplots(2, 2, figsize=(10, 7), layout="constrained")
    for stage in range(2):
        temperature = [point["stage_temperature_K"][stage] for point in converged]
        axes[0, stage].scatter(
            temperature, [point["heater_power_W"][stage] for point in converged], label="Measured heater"
        )
        axes[0, stage].scatter(
            temperature, [point["predicted_cooling_W"][stage] for point in converged], label="Prediction", marker="x"
        )
        axes[0, stage].set(title=f"Stage {stage + 1}", xlabel="Measured temperature (K)", ylabel="Power (W)")
        axes[0, stage].legend()
        axes[1, stage].axhline(0, color="0.5", linewidth=0.8)
        axes[1, stage].scatter(temperature, [point["residual_W"][stage] for point in converged])
        axes[1, stage].set(xlabel="Measured temperature (K)", ylabel="Net stage imbalance (W)")
    figure.suptitle("August prediction: unchanged May fit, two compressors assumed")
    figure.savefig(destination, dpi=160)
    plt.close(figure)


def main() -> None:
    """Runs both cooldowns and the August map prediction, archiving their evidence."""
    import matplotlib

    matplotlib.use("Agg")
    inputs = read_inputs()
    report = {"provisional": True, "runs": {}, "sources": inputs["sources"]}
    histories = {}
    for run in inputs["runs"]:
        print(f"Simulating {run}...", flush=True)
        start = time.perf_counter()
        model = build_model(run)
        result = simulate_cooldown(model)
        measured = read_measurements(run)
        histories[run] = result.data.temperature
        result.data.to_netcdf(HERE / f"{run}_cooldown.nc")
        report["runs"][run] = {
            "status": result.status,
            "elapsed_s": time.perf_counter() - start,
            "origin": measured.attrs["origin"],
            "end_exclusive": measured.attrs["end_exclusive"],
            "end_time_s": float(result.data.time[-1]),
            "initial_temperature_K": float(model.initial_temperature[0]),
            "compressors": len(model.graphs[0].links[0].compressors),
            "helium_charge_kg": model.graphs[0].links[0].charge,
            "final_stage_temperature_K": result.data.temperature.sel(node=["Stage 1", "Stage 2"]).values[-1].tolist(),
            "measured_final_stage_temperature_K": [float(measured[name][-1]) for name in ("T1", "T2")],
            "comparison": compare_cooldown(measured, result.data.temperature),
            "final_state": describe_state(model, result.data.temperature.values[-1]),
        }
        print(json.dumps(report["runs"][run]), flush=True)
    print("Predicting August load map...", flush=True)
    load_map = read_load_map(DATA / inputs["runs"]["august"]["load_map"])
    report["august_load_map"] = predict_load_map(build_model("august"), histories["august"], load_map)
    tracked = [
        HERE / "inputs.toml",
        HERE / "parameters.toml",
        HERE / "may_fit_snapshot.json",
        HERE / "run.py",
        ROOT / inputs["geometry"],
        ROOT / "src/pulse_tube_cooldown/library/heads/data/pt1_1.toml",
    ]
    tracked.extend(DATA / config["cooldown"] for config in inputs["runs"].values())
    tracked.append(DATA / inputs["runs"]["august"]["load_map"])
    tracked.extend((ROOT / "src/pulse_tube_cooldown").rglob("*.py"))
    tracked.extend((ROOT / "src/pulse_tube_cooldown/library/materials/data").glob("*.toml"))
    report["sha256"] = {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest() for path in tracked}
    report["git_revision"] = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    (HERE / "results.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    plot_cooldowns(histories, HERE / "cooldown_comparison.png")
    plot_load_map(report["august_load_map"], HERE / "august_load_map.png")
    print("Wrote results.json, two netCDF histories, and comparison figures.", flush=True)


if __name__ == "__main__":
    main()
