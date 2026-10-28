"""The canonical toy run through the public calls and the `run-ptc` command, checked against a reference."""

import shutil
import subprocess
import sys
import tomllib
from functools import cache
from pathlib import Path

import numpy as np
import pytest
import xarray as xr
from scipy.integrate import quad
from scipy.optimize import brentq

import pulse_tube_cooldown as ptc
from pulse_tube_cooldown.library.materials.load import load_material

_CONFIG = Path(__file__).parents[2] / "simulations" / "validation" / "toy" / "toy.toml"
_RUN_PTC = shutil.which("run-ptc", path=Path(sys.executable).parent)

# The Plate cools toward the Bath, held fixed, through Rod A, and through Rods A and B once it reaches the switch
# temperature. With one free node, the time to cool from T0 to T is
# t(T) = integral from T to T0 of m * c(T') / (G * (K(T') - K(T_bath))) dT', where K is the rods' conductivity
# integral and G the sum of A / L over the active rods.
with _CONFIG.open("rb") as _file:
    _TOML = tomllib.load(_file)

_NODES = {node["name"]: node for node in _TOML["fridge"]["nodes"]}
_LINKS = {link["name"]: link for link in _TOML["fridge"]["links"]}
((_PLATE_MATERIAL, _MASS),) = _NODES["Plate"]["composition"]
(_ROD_MATERIAL,) = {link["material"] for link in _LINKS.values()}
_COPPER, _STEEL = load_material(_PLATE_MATERIAL), load_material(_ROD_MATERIAL)
_CONDUCTANCE = {name: link["area"] / link["length"] for name, link in _LINKS.items()}

_SIMULATION = _TOML["simulation"]
_INITIAL, _BATH = _SIMULATION["initial_temperature"]
_MAX_TIME, _SAVE_EVERY = _SIMULATION["max_time"], _SIMULATION["save_every"]
_SWITCH = _TOML["procedure"]["phases"][0]["until"]["threshold"]
_SETTLED = _BATH + 1e-9


def _rate(temperature, conductance):
    capacity = _MASS * _COPPER.specific_heat(np.array([temperature]))[0]
    integral = _STEEL.conductivity_integral(np.array([temperature, _BATH]))
    return capacity / (conductance * (integral[0] - integral[1]))


@cache
def _reference_time(temperature):
    single, both = _CONDUCTANCE["Rod A"], _CONDUCTANCE["Rod A"] + _CONDUCTANCE["Rod B"]
    if temperature >= _SWITCH:
        return quad(_rate, temperature, _INITIAL, args=(single,), epsabs=0, epsrel=1e-10, limit=200)[0]
    return (
        _reference_time(_SWITCH) + quad(_rate, temperature, _SWITCH, args=(both,), epsabs=0, epsrel=1e-10, limit=200)[0]
    )


def _reference_temperature(time):
    if time == 0.0:
        return _INITIAL
    if time >= _reference_time(_SETTLED):
        return _BATH
    return brentq(lambda temperature: _reference_time(temperature) - time, _SETTLED, _INITIAL, xtol=1e-12)


@pytest.fixture(scope="module")
def saved(tmp_path_factory):
    directory = tmp_path_factory.mktemp("toy")
    api, cli = directory / "api.nc", directory / "cli.nc"

    model = ptc.configure_model(_CONFIG)
    result = ptc.simulate_cooldown(model)
    ptc.save_output(result, api)

    completed = subprocess.run(
        [_RUN_PTC, "--config", str(_CONFIG), "--output", str(cli)], capture_output=True, text=True
    )
    assert completed.returncode == 0, completed.stderr

    with xr.open_dataset(api) as from_api, xr.open_dataset(cli) as from_cli:
        yield from_api.load(), from_cli.load()


def test_public_calls_and_cli_save_identical_results(saved):
    from_api, from_cli = saved

    xr.testing.assert_identical(from_api, from_cli)


def test_run_ends_at_the_time_limit(saved):
    data, _ = saved

    assert data.attrs["status"] == "max_time"
    assert data["time"].values[-1] == _MAX_TIME


def test_samples_fall_every_save_interval_from_each_phase_start(saved):
    data, _ = saved
    switch = data["time"].values[np.flatnonzero(data["phase"].values == "both rods")[0]]

    expected = np.concatenate(
        [np.arange(0.0, switch, _SAVE_EVERY), [switch], np.arange(switch, _MAX_TIME, _SAVE_EVERY), [_MAX_TIME]]
    )

    np.testing.assert_allclose(data["time"].values, expected)


def test_phase_switch_occurs_at_the_reference_time_with_both_sides_recorded(saved):
    data, _ = saved
    phase = data["phase"].values
    last_single = np.flatnonzero(phase == "single rod")[-1]

    switch = data.isel(time=[last_single, last_single + 1])

    assert list(switch["phase"].values) == ["single rod", "both rods"]
    np.testing.assert_allclose(switch["time"].values, _reference_time(_SWITCH), rtol=1e-6)
    np.testing.assert_allclose(switch["temperature"].sel(node="Plate").values, _SWITCH, rtol=1e-6)


def test_plate_follows_the_reference_cooldown_and_settles_at_the_bath(saved):
    data, _ = saved
    plate = data["temperature"].sel(node="Plate").values

    expected = [_reference_temperature(time) for time in data["time"].values]

    assert expected[-1] == _BATH
    np.testing.assert_allclose(plate, expected, rtol=1e-6)


def test_bath_stays_at_its_fixed_temperature(saved):
    data, _ = saved

    np.testing.assert_array_equal(data["temperature"].sel(node="Bath").values, _BATH)


def test_each_rod_conducts_heat_out_of_the_plate_while_active(saved):
    data, _ = saved
    single = data["phase"].values == "single rod"
    plate = data["temperature"].sel(node="Plate").values
    integral = _STEEL.conductivity_integral(np.stack([plate, np.full_like(plate, _BATH)]))
    conducted = integral[1] - integral[0]

    flow = data["heat_flow"].sel(node="Plate")

    np.testing.assert_allclose(flow.sel(link="Rod A"), _CONDUCTANCE["Rod A"] * conducted, rtol=1e-10, atol=1e-12)
    np.testing.assert_allclose(
        flow.sel(link="Rod B"), np.where(single, np.nan, _CONDUCTANCE["Rod B"] * conducted), rtol=1e-10, atol=1e-12
    )


def test_cooldown_matches_pinned_values(saved):
    data, _ = saved
    plate = data["temperature"].sel(node="Plate")
    both_rods = plate.where(data["phase"] == "both rods", drop=True)

    np.testing.assert_allclose(both_rods["time"].values[0], 49458.50005439385, rtol=1e-6)
    np.testing.assert_allclose(
        plate.where(data["time"].isin([21600.0, 43200.0]), drop=True), [179.7211747, 113.69263274], rtol=1e-6
    )
    np.testing.assert_allclose(both_rods.isel(time=6), 49.95611153, rtol=1e-6)
