import subprocess
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pytest
import xarray as xr

import pulse_tube_cooldown as ptc
from pulse_tube_cooldown.core.types import EnergyBalanceTerms, Material, Model, Node, Phase, Procedure


@pytest.fixture
def config_path():
    """Return the configuration path for interface tests."""
    return Path(__file__).parent / "fixtures" / "config_for_interface_tests.toml"

@pytest.fixture
def result(config_path):
    """Return the simulation result for the interface test configuration."""
    model = ptc.configure_model(config_path)
    return ptc.simulate_cooldown(model)


@dataclass(frozen=True)
class CustomLink:
    name: str
    nodes: tuple[str, str]
    conductance: float

    def calc_energy_balance_terms(self, temperature):
        """Return conductive heat flows into the two linked nodes."""
        flow = self.conductance * (temperature[1] - temperature[0])
        return EnergyBalanceTerms(heat_flow=np.array([flow, -flow]))


def test_run_from_command_line_interface(config_path, tmp_path):
    """Verify that the CLI runs a configured model and writes an output file."""
    output_path = tmp_path / "result.nc"

    result = subprocess.run(
        ["run-ptc", "--config", config_path, "--output", output_path], capture_output=True, text=True
    )

    assert result.returncode == 0, result.stderr
    assert output_path.is_file()


def test_run_from_config_file(config_path, tmp_path):
    """Verify that the public API runs a configured model and saves its result."""
    output_path = tmp_path / "result.nc"

    model = ptc.configure_model(config_path)
    result = ptc.simulate_cooldown(model)
    ptc.save_output(result, output_path)

    assert output_path.is_file()


def test_run_from_user_defined_objects(tmp_path):
    """Verify that a model built from user-defined objects can run and save."""
    custom_material = Material(
        name="Constant heat capacity",
        density=1000.0,
        specific_heat=lambda temperature: np.full_like(temperature, 100.0),
        specific_heat_integral=lambda temperature: 100.0 * temperature,
        thermal_conductivity=None,
        conductivity_integral=None,
    )
    nodes = (
        Node(name="Plate", composition=((custom_material, 1.0),)),
        Node(name="Bath", composition=((custom_material, 1.0),)),
    )
    link = CustomLink(
        name="Strap",
        nodes=("Plate", "Bath"),
        conductance=1.0,
    )
    phase = Phase(name="cooling", links=(link,), until=None)
    model = Model(
        meta={"experiment": "Custom model"},
        nodes=nodes,
        procedure=Procedure(phases=(phase,), stop_condition=None),
        initial_temperature=np.array([300.0, 250.0]),
        fixed_boundary_nodes=(1,),
        solver_options={"max_time": 10.0, "save_every": 1.0},
    )
    output_path = tmp_path / "result.nc"

    result = ptc.simulate_cooldown(model)
    ptc.save_output(result, output_path)

    assert output_path.is_file()


def test_result_object_contract(result):
    """Verify result data, units, phase-dependent heat flows, and metadata."""
    data = result.data

    plate_temperature = data["temperature"].sel(node="Plate")
    rod_a_heat_flow = data["heat_flow"].sel(link="Rod A", node="Plate")

    assert plate_temperature.sizes["time"] > 1
    assert np.isfinite(plate_temperature).all()
    assert np.isfinite(rod_a_heat_flow).all()

    assert data["time"].attrs["units"] == "s"
    assert plate_temperature.attrs["units"] == "K"
    assert rod_a_heat_flow.attrs["units"] == "W"

    phases = list(dict.fromkeys(data["phase"].values))
    assert phases == ["single rod", "both rods"]

    rod_b_heat_flow = data["heat_flow"].sel(link="Rod B", node="Plate")
    before = rod_b_heat_flow.where(data["phase"] == "single rod", drop=True)
    after = rod_b_heat_flow.where(data["phase"] == "both rods", drop=True)

    assert before.isnull().all()
    assert np.isfinite(after).all()

    for name, value in result.model.meta.items():
        assert data.attrs[name] == value

def test_saved_data_roundtrips(result, tmp_path):
    """Verify that saving and loading preserves the result dataset exactly."""
    output_path = tmp_path / "result.nc"

    ptc.save_output(result, output_path)

    with xr.open_dataset(output_path) as saved:
        restored = saved.load()

    xr.testing.assert_identical(restored, result.data)
