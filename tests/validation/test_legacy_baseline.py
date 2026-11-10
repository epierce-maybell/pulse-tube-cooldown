"""The legacy configurations for SN09 and FP1, and the archived legacy baseline."""

import tomllib
from pathlib import Path

import pytest
import xarray as xr

from pulse_tube_cooldown.measurements import log_temperature_residual, read_refrigerator_cooldown
from simulations.validation.legacy_baseline.run import legacy_system, with_run_options

ARCHIVE = Path(__file__).parents[2] / "simulations" / "validation" / "legacy_baseline"


def _fridge():
    """A legacy saved system shaped like `Fridge`."""
    copper = [{"material": "OFHC Copper (RRR=100)", "mass": "1.0"}]
    return {
        "coolers": ["PT1.1 (1X F-70)"],
        "cooler": "PT1.1 (1X F-70)",
        "T_init": "295.0",
        "stages": {
            node: {"materials": copper, "target": "14.1"} for node in ["PRP", "RGP", "CFP", "STP", "ICP", "MXP"]
        },
        "linked_masses": [
            {
                "name": "DHX",
                "stage": "MXP",
                "target_T": "20.0",
                "materials": [{"material": "304 Stainless Steel", "mass": "0.4"}],
                "links": [{"material": "304 Stainless Steel", "area_cm2": "0.4", "length_cm": "9.0", "qty": "1"}],
            }
        ],
        "radiation": {"enabled": True, "mli_layers": "15"},
    }


def test_sn09_system():
    """SN09's BOM masses by node from `data/sn09/bom_by_node.csv` and one PT1.1 on two compressors."""
    config = legacy_system(_fridge(), "sn09")

    assert config["stages"]["CFP"] == {
        "materials": [
            {"material": "OFHC Copper (RRR=100)", "mass": 24.286},
            {"material": "304 Stainless Steel", "mass": 0.06},
            {"material": "6061-T6 Aluminum", "mass": 0.176},
            {"material": "G-10 Fiberglass Epoxy (norm)", "mass": 0.046},
        ],
        "target": "14.1",
    }
    [dhx] = config["linked_masses"]
    assert dhx["materials"] == [{"material": "304 Stainless Steel", "mass": 0.758}]
    assert dhx["links"] == _fridge()["linked_masses"][0]["links"]
    assert dhx["target_T"] == "20.0"
    assert config["coolers"] == ["PT1.1"]
    assert config["radiation"] == _fridge()["radiation"]
    assert config["T_init"] == "295.0"


def test_fp1_system():
    """FP1's two PT1.5 heads, each on its own compressor, and its STP copper."""
    config = legacy_system(_fridge(), "fp1")

    assert config["coolers"] == ["PT1.5 (1X F-70)", "PT1.5 (1X F-70)"]
    assert {"material": "OFHC Copper (RRR=100)", "mass": 31.556} in config["stages"]["STP"]["materials"]


def test_run_options():
    """`gas_gap_tables = false` removes every gas gap's table; `linked_mass_gas_links` replaces DHX's gas links."""
    tabled = {"gas": "N2", "k_table": [[5.0, 0.0002], [300.0, 0.03]]}
    config = _fridge() | {
        key: [dict(tabled)] for key in ["gas_links_pr", "gas_links_rc", "gas_links_cs", "gas_links_si", "gas_links_im"]
    }
    dhx_links = [{"gas": "He-4", "pressure_torr": 750.0, "p_mode": "T_avg"}]

    options = with_run_options(config, {"gas_gap_tables": False, "linked_mass_gas_links": {"DHX": dhx_links}})

    assert options["gas_links_cs"] == [{"gas": "N2"}]
    assert options["linked_masses"][0]["gas_links"] == dhx_links
    assert config["gas_links_cs"] == [tabled]
    assert "gas_links" not in config["linked_masses"][0]
    assert with_run_options(config, {}) == config


@pytest.mark.parametrize("run", ["sn09_tuned", "fp1_generous"])
def test_archived_scores(run):
    """Each archived run covers its window and scores as recorded."""
    with (ARCHIVE / "scores.toml").open("rb") as file:
        scores = tomllib.load(file)[run]
    with xr.open_dataset(ARCHIVE / f"{run}.nc") as trace:
        measured = read_refrigerator_cooldown(Path(__file__).parents[2] / "data" / scores["system"] / "cooldown.toml")
        rms, rms_by_node = log_temperature_residual(trace.temperature.load(), measured)

    assert scores["rms"] == pytest.approx(rms)
    assert scores["rms_by_node"] == pytest.approx(rms_by_node)
