"""Tests for the public functions and the `run-ptc` command."""

import shutil
import subprocess
import sys
from pathlib import Path

import xarray as xr

import pulse_tube_cooldown as ptc
from pulse_tube_cooldown.__main__ import main
from pulse_tube_cooldown.core.types import Result

_CONFIG = """
[meta]
experiment = "Toy CLI run"
version = "abc1234"

[[fridge.nodes]]
name = "Plate"
composition = [["OFHC Cu RRR100", 1.0]]

[[fridge.nodes]]
name = "Bath"
composition = [["OFHC Cu RRR100", 1.0]]

[[fridge.links]]
name = "Rod A"
type = "conduction"
nodes = ["Plate", "Bath"]
material = "304 SS"
area = 1.0e-4
length = 0.1

[[fridge.links]]
name = "Rod B"
type = "conduction"
nodes = ["Plate", "Bath"]
material = "304 SS"
area = 1.0e-4
length = 0.1

[[procedure.phases]]
name = "single rod"
disable = ["Rod B"]
until = { condition = "gas_pushing_begins", node = "Plate", threshold = 290.0 }

[[procedure.phases]]
name = "both rods"
enable = ["Rod B"]

[simulation]
max_time = 3600.0
save_every = 600.0
initial_temperature = [300.0, 4.0]
fixed_boundary_nodes = ["Bath"]
"""

_RUN_PTC = shutil.which("run-ptc", path=Path(sys.executable).parent)


def _config(tmp_path, text=_CONFIG):
    path = tmp_path / "toy.toml"
    path.write_text(text)
    return path


def _run_ptc(*args):
    return subprocess.run([_RUN_PTC, *map(str, args)], capture_output=True, text=True)


def test_public_functions_run_the_canonical_calls(tmp_path):
    output = tmp_path / "run.nc"

    model = ptc.configure_model(_config(tmp_path))
    result = ptc.simulate_cooldown(model)
    ptc.save_output(result, output)

    with xr.open_dataset(output) as saved:
        assert saved.attrs["experiment"] == "Toy CLI run"
        assert list(dict.fromkeys(saved["phase"].values)) == ["single rod", "both rods"]


def test_installed_cli_configures_simulates_and_saves_a_run(tmp_path):
    output = tmp_path / "run.nc"

    completed = _run_ptc("--config", _config(tmp_path), "--output", output)

    assert completed.returncode == 0, completed.stderr
    with xr.open_dataset(output) as saved:
        assert saved.attrs["experiment"] == "Toy CLI run"
        assert saved.attrs["status"] == "max_time"
        assert list(dict.fromkeys(saved["phase"].values)) == ["single rod", "both rods"]
        assert saved["time"].values[-1] == 3600.0


def test_bad_configuration_fails_with_its_error_and_saves_nothing(tmp_path):
    output = tmp_path / "run.nc"
    config = _config(tmp_path, _CONFIG.replace('material = "304 SS"', 'material = "Unobtanium"', 1))

    completed = _run_ptc("--config", config, "--output", output)

    assert completed.returncode != 0
    assert "Unobtanium" in completed.stderr
    assert "Traceback" not in completed.stderr
    assert not output.exists()


def test_missing_configuration_file_fails_naming_the_path(tmp_path):
    completed = _run_ptc("--config", tmp_path / "absent.toml", "--output", tmp_path / "run.nc")

    assert completed.returncode != 0
    assert "absent.toml" in completed.stderr
    assert "Traceback" not in completed.stderr


def test_malformed_toml_fails_with_the_parse_error(tmp_path):
    completed = _run_ptc("--config", _config(tmp_path, "[fridge\n"), "--output", tmp_path / "run.nc")

    assert completed.returncode != 0
    assert "line 1" in completed.stderr
    assert "Traceback" not in completed.stderr


def test_unwritable_output_fails_naming_the_path(tmp_path):
    output = tmp_path / "absent" / "run.nc"

    completed = _run_ptc("--config", _config(tmp_path), "--output", output)

    assert completed.returncode != 0
    assert "Traceback" not in completed.stderr
    assert "run.nc" in completed.stderr


def test_solver_failure_saves_the_run_and_fails(tmp_path, monkeypatch, capsys):
    output = tmp_path / "run.nc"
    simulate_cooldown = ptc.simulate_cooldown

    def failing(model):
        result = simulate_cooldown(model)
        data = result.data.assign_attrs(status="solver_failure", message="Required step size is too small.")
        return Result(status="solver_failure", model=model, data=data)

    monkeypatch.setattr(ptc, "simulate_cooldown", failing)

    code = main(["--config", str(_config(tmp_path)), "--output", str(output)])

    assert code != 0
    assert "Required step size is too small." in capsys.readouterr().err
    with xr.open_dataset(output) as saved:
        assert saved.attrs["status"] == "solver_failure"
