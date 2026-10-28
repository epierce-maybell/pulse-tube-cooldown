"""The canonical toy run through the public calls and the `run-ptc` command, checked against a reference."""

import shutil
import subprocess
import sys
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

# The Plate (1 kg Cu) cools toward the Bath, fixed at 4 K, through one 304 SS rod of A / L = 1e-3 m, and through two
# once it reaches 100 K. With one free node, the time to cool from T0 to T is
# t(T) = integral from T to T0 of m * c(T') / (n * (A / L) * (K(T') - K(T_bath))) dT', where K is the conductivity
# integral and n the number of rods.
_COPPER = load_material("OFHC Cu RRR100")
_STEEL = load_material("304 SS")
_INITIAL, _BATH, _SWITCH = 300.0, 4.0, 100.0


def _rate(temperature, rods):
    capacity = _COPPER.specific_heat(np.array([temperature]))[0]
    integral = _STEEL.conductivity_integral(np.array([temperature, _BATH]))
    return capacity / (rods * 1.0e-3 * (integral[0] - integral[1]))


def _reference_time(temperature):
    if temperature >= _SWITCH:
        return quad(_rate, temperature, _INITIAL, args=(1,), epsabs=0, epsrel=1e-10, limit=200)[0]
    return _reference_time(_SWITCH) + quad(_rate, temperature, _SWITCH, args=(2,), epsabs=0, epsrel=1e-10, limit=200)[0]


def _reference_temperature(time):
    if time == 0.0:
        return _INITIAL
    return brentq(lambda temperature: _reference_time(temperature) - time, _BATH + 1e-9, _INITIAL, xtol=1e-12)


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
    assert data["time"].values[-1] == 75600.0


def test_phase_switch_occurs_at_the_reference_time_with_both_sides_recorded(saved):
    data, _ = saved
    phase = data["phase"].values
    last_single = np.flatnonzero(phase == "single rod")[-1]

    switch = data.isel(time=[last_single, last_single + 1])

    assert list(switch["phase"].values) == ["single rod", "both rods"]
    np.testing.assert_allclose(switch["time"].values, _reference_time(_SWITCH), rtol=1e-6)
    np.testing.assert_allclose(switch["temperature"].sel(node="Plate").values, _SWITCH, rtol=1e-6)


def test_plate_follows_the_reference_cooldown(saved):
    data, _ = saved

    expected = [_reference_temperature(time) for time in data["time"].values]

    np.testing.assert_allclose(data["temperature"].sel(node="Plate").values, expected, rtol=1e-6)


def test_bath_stays_at_its_fixed_temperature(saved):
    data, _ = saved

    np.testing.assert_array_equal(data["temperature"].sel(node="Bath").values, _BATH)


def test_second_rod_carries_heat_only_after_the_switch(saved):
    data, _ = saved
    single = data["phase"] == "single rod"
    rod_a = data["heat_flow"].sel(link="Rod A", node="Plate")
    rod_b = data["heat_flow"].sel(link="Rod B", node="Plate")

    assert rod_b.where(single, drop=True).isnull().all()
    np.testing.assert_allclose(rod_b.where(~single, drop=True), rod_a.where(~single, drop=True))
    assert (rod_a < 0).all()
