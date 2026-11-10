"""Derives the retained SN09 and FP1 legacy baselines.

Run from the repository root:
    python -m simulations.validation.legacy_baseline.run

inputs.toml selects the legacy commit, source inputs, and the two baseline runs.
Results, complete run configurations, scores, and figures are written beside it.
"""

import argparse
import copy
import csv
import importlib
import io
import json
import subprocess
import sys
import tempfile
import tomllib
import zipfile
from pathlib import Path

import numpy as np
import xarray as xr
from matplotlib.figure import Figure

from pulse_tube_cooldown.measurements import log_temperature_residual, read_refrigerator_cooldown


def _read_measurements(system: str) -> xr.Dataset:
    return read_refrigerator_cooldown(Path(__file__).parents[3] / "data" / system / "cooldown.toml")


def legacy_system(base: dict, system: str) -> dict:
    """Returns a legacy saved system for `system`.

    It is `base`, a saved system in the legacy `ptc_systems.json` format, with the plate and DHX materials of
    `data/<system>/bom_by_node.csv` and the system's coolers from `inputs.toml`.
    """
    with (Path(__file__).with_name("inputs.toml")).open("rb") as file:
        inputs = tomllib.load(file)
    with (Path(__file__).parents[3] / "data" / system / "bom_by_node.csv").open(encoding="utf-8") as file:
        materials = {
            row["node"]: [
                {"material": material, "mass": float(row[column])}
                for column, material in inputs["materials"].items()
                if float(row[column]) > 0.0
            ]
            for row in csv.DictReader(file)
        }

    config = copy.deepcopy(base)
    for node, stage in config["stages"].items():
        stage["materials"] = materials[node]
    for linked_mass in config["linked_masses"]:
        linked_mass["materials"] = materials[linked_mass["name"]]
    config["coolers"] = inputs["coolers"][system]
    config["cooler"] = config["coolers"][0]
    return config


def with_run_options(config: dict, run: dict) -> dict:
    """Returns `config` with a run's options from `inputs.toml` applied.

    Args:
        config: A saved system in the legacy `ptc_systems.json` format.
        run: The run's entry. `gas_gap_tables = false` removes the tabulated conductance (`k_table`) from every
            gas gap, leaving the legacy internal gas model. `linked_mass_gas_links` replaces the gas links of each
            named linked mass.
    """
    config = copy.deepcopy(config)
    if not run.get("gas_gap_tables", True):
        for key in ["gas_links_pr", "gas_links_rc", "gas_links_cs", "gas_links_si", "gas_links_im"]:
            for link in config[key]:
                link.pop("k_table", None)
    for linked_mass in config["linked_masses"]:
        if linked_mass["name"] in run.get("linked_mass_gas_links", {}):
            linked_mass["gas_links"] = copy.deepcopy(run["linked_mass_gas_links"][linked_mass["name"]])
    return config


def export_commit(repository: Path, commit: str, directory: Path) -> Path:
    """Writes `repository`'s tree at `commit` into `directory` and returns `directory`."""
    archive = subprocess.run(
        ["git", "archive", "--format=zip", commit], cwd=repository, capture_output=True, check=True
    ).stdout
    zipfile.ZipFile(io.BytesIO(archive)).extractall(directory)
    return directory


def simulate_legacy(
    config: dict,
    initial_temperature: float,
    max_time: float,
    legacy: Path = Path(__file__).parents[4] / "pulse_tube_cooldown_simulation",
) -> xr.DataArray:
    """Runs the legacy engine on a saved system, as the legacy GUI does, from `initial_temperature` to `max_time`.

    Every node's target is 0 K, so the run ends at `max_time` and no linked mass is held at its target.

    Args:
        config: A saved system in the legacy `ptc_systems.json` format.
        initial_temperature: Every node's starting temperature and the ambient temperature, in K.
        max_time: The run's length, in s.
        legacy: The legacy repository.

    Returns:
        `temperature(time, node)` in K, time in s, for the six plates and the linked masses.
    """
    sys.path.insert(0, str(legacy))
    engine = importlib.import_module("ptc_cooldown_simulator")

    def gas_links(links):
        pascals_per_torr = 133.322
        return [
            {
                "gas": link["gas"],
                "pressure": float(link["pressure_torr"]) * pascals_per_torr,
                "area": float(link["area_cm2"]) * 1e-4,
                "length": float(link["length_cm"]) * 1e-2,
                "shutoff_T": float(link["shutoff_K"]),
                "p_mode": link["p_mode"],
                "series_links": None,
            }
            | ({"k_table": sorted([float(t), float(k)] for t, k in link["k_table"])} if "k_table" in link else {})
            for link in links
        ]

    def material_list(materials):
        return [(material["material"], float(material["mass"])) for material in materials]

    linked_masses = [
        {
            "name": linked_mass["name"],
            "stage": linked_mass["stage"],
            "target_T": 0.0,
            "fixed_heat_W": float(linked_mass["fixed_heat_W"]),
            "mat_list": material_list(linked_mass["materials"]),
            "links": [
                {
                    "material": link["material"],
                    "area": float(link["area_cm2"]) * 1e-4,
                    "length": float(link["length_cm"]) * 1e-2,
                    "qty": int(float(link["qty"])),
                }
                for link in linked_mass["links"]
            ],
            "lib_links": [],
            "gas_links": gas_links(linked_mass.get("gas_links", [])),
        }
        for linked_mass in config["linked_masses"]
    ]
    radiative_loads = [
        {
            "name": load["name"],
            "stage": load["stage"],
            "hot_mode": load["hot_mode"],
            "hot_T": float(load["hot_T"]),
            "area": float(load["area"]),
            "e_hot": float(load["e_hot"]),
            "e_cold": float(load["e_cold"]),
            "geometry": load["geometry"],
        }
        for load in config["rad_loads"]
    ]
    conductive_loads = [
        {
            "name": load["name"],
            "stage": load["stage"],
            "warm_mode": load["warm_mode"],
            "warm_T": float(load["warm_T"]),
            "material": load["material"],
            "area": float(load["area_cm2"]) * 1e-4,
            "length": float(load["length_cm"]) * 1e-2,
            "qty": int(float(load["qty"])),
        }
        for load in config["cond_loads"]
    ]
    fixed_heat_loads = [
        {"name": load["name"], "stage": load["stage"], "power_W": float(load["power_W"])}
        for load in config["fixed_heat_loads"]
    ]
    radiation = config["radiation"]
    multi_layer_insulation = (
        {
            "n_layers": int(float(radiation["mli_layers"])),
            "layer_density": float(radiation["mli_density"]),
            "e_mli": float(radiation["mli_emissivity"]),
            "prp_area": float(radiation["prp_area"]),
            "cfp_area": float(radiation["cfp_area"]),
            "prp_inner_emissivity": float(radiation["prp_inner_emissivity"]),
            "cfp_emissivity": float(radiation["cfp_emissivity"]),
            "mli_degradation": float(radiation["mli_degradation"]),
            "cfp_geometry": radiation["cfp_geometry"],
        }
        if radiation["enabled"]
        else None
    )
    exchange_gas = config["precool_xg"]

    time, stage_temperatures, linked_mass_temperatures, _ = engine.simulate_cooldown(
        config["coolers"],
        {stage: material_list(config["stages"][stage]["materials"]) for stage in engine.STAGES},
        {stage: 0.0 for stage in engine.STAGES},
        initial_temperature,
        include_par=config["parasitic"],
        gas_links_prp_rgp=gas_links(config["gas_links_pr"]),
        gas_links_rgp_cfp=gas_links(config["gas_links_rc"]),
        gas_links_cfp_stp=gas_links(config["gas_links_cs"]),
        gas_links_stp_icp=gas_links(config["gas_links_si"]),
        gas_links_icp_mxp=gas_links(config["gas_links_im"]),
        linked_masses=linked_masses,
        mli_cfg=multi_layer_insulation,
        rad_loads=radiative_loads,
        cond_loads=conductive_loads,
        fixed_heat_loads=fixed_heat_loads,
        max_time_hrs=max_time / 3600.0,
        contact_factor=float(config["contact_factor"]),
        exchange_gas=exchange_gas["enabled"],
        xg_gas=exchange_gas["gas"],
        xg_area=float(exchange_gas["area"]),
        xg_cfp_cutoff=float(exchange_gas["cfp_cutoff"]),
        xg_scale=float(exchange_gas["scale"]),
    )

    nodes = list(engine.STAGES) + [linked_mass["name"] for linked_mass in linked_masses]
    temperature = np.column_stack([stage_temperatures[stage] for stage in engine.STAGES] + linked_mass_temperatures)
    return xr.DataArray(
        temperature,
        dims=("time", "node"),
        coords={"time": ("time", time, {"units": "s"}), "node": nodes},
        attrs={"units": "K"},
    )


def plot_baseline(
    measured: xr.Dataset,
    temperature: xr.DataArray,
    rms_by_node: dict[str, float],
    *,
    label: str,
) -> Figure:
    """Plots one baseline over the measured cooldown, with a panel per thermometer."""
    nodes = measured.node.values
    figure = Figure(figsize=(9.0, 2.0 * len(nodes)), layout="constrained")
    axes = figure.subplots(len(nodes), 1, sharex=True, squeeze=False)[:, 0]
    for ax, node, scored in zip(axes, nodes, measured.scored.values, strict=True):
        ax.plot(measured.time / 3600.0, measured.temperature.sel(node=node), color="#2b2b2b", lw=1.5, label="Measured")
        curve_label = f"{label}, RMS {rms_by_node[node]:.2f}" if scored else label
        ax.plot(temperature.time / 3600.0, temperature.sel(node=node), color="#2a78d6", lw=1.5, label=curve_label)
        ax.set_yscale("log")
        ax.set_ylabel(f"{node} (K)" if scored else f"{node} (K), not scored")
        ax.grid(which="major", color="0.85", lw=0.5)
        ax.legend(loc="best", fontsize="small", frameon=False)
    axes[-1].set_xlim(0.0, float(measured.time[-1]) / 3600.0)
    axes[-1].set_xlabel("Time since compressor start (h)")
    return figure


def main(
    repository: Path = Path(__file__).parents[4] / "pulse_tube_cooldown_simulation",
    *,
    output: Path | None = None,
) -> None:
    """Derives the selected baselines and saves their configurations, traces, scores, and plots.

    Args:
        repository: The legacy repository. Every run uses its tree at `inputs.toml`'s `legacy_commit`, whatever
            is checked out.
        output: Destination directory. Defaults to the directory containing this script and its inputs.
    """
    source = Path(__file__).with_name("inputs.toml")
    archive = Path(__file__).parent if output is None else output
    archive.mkdir(parents=True, exist_ok=True)
    with source.open("rb") as file:
        inputs = tomllib.load(file)
    if archive.resolve() != source.parent.resolve():
        (archive / source.name).write_bytes(source.read_bytes())
    commit = inputs["legacy_commit"]

    scores = [
        f'legacy_commit = "{commit}"',
        "# Each run sets every legacy target to 0 K, so it ends at max_time and no linked mass is held at its target.",
        "# initial_temperature is the mean of the scored nodes' first valid readings, in K; max_time is in s.",
    ]
    with tempfile.TemporaryDirectory() as directory:
        legacy = export_commit(repository, commit, Path(directory))
        with (legacy / "ptc_systems.json").open(encoding="utf-8") as file:
            base = json.load(file)[inputs["base_system"]]
        for name, run in inputs["runs"].items():
            system = run["system"]
            measured = _read_measurements(system)
            scored = measured.temperature.sel(node=measured.node[measured.scored])
            initial_temperature = float(
                np.mean([readings.dropna("time")[0] for readings in scored.transpose("node", ...)])
            )
            max_time = float(measured.time[-1]) + 3600.0

            config = with_run_options(legacy_system(base, system) if run["system_inputs"] else base, run)
            config["T_init"] = initial_temperature
            config["max_time_hrs"] = max_time / 3600.0
            temperature = simulate_legacy(config, initial_temperature, max_time, legacy)
            rms, rms_by_node = log_temperature_residual(temperature, measured)

            (archive / f"{name}.json").write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
            temperature.to_dataset(name="temperature").to_netcdf(archive / f"{name}.nc")
            scores += [
                f"\n[{name}]",
                f'system = "{system}"',
                f"rms = {rms}",
                f"initial_temperature = {initial_temperature}",
                f"max_time = {max_time}",
                f"\n[{name}.rms_by_node]",
                *(f"{node} = {value}" for node, value in rms_by_node.items()),
            ]
            print(f"{name}: rms {rms:.4f}", {node: round(value, 4) for node, value in rms_by_node.items()})
            figure = plot_baseline(measured, temperature, rms_by_node, label=name)
            figure.suptitle(f"{system.upper()}: legacy baseline and measured cooldown")
            figure.savefig(archive / f"{system}.png", dpi=150)

    (archive / "scores.toml").write_text("\n".join(scores) + "\n", encoding="utf-8")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--legacy-repository", type=Path, default=Path(__file__).parents[4] / "pulse_tube_cooldown_simulation"
    )
    parser.add_argument("--output", type=Path, help="Directory for derived baseline evidence")
    args = parser.parse_args()
    main(args.legacy_repository, output=args.output)
