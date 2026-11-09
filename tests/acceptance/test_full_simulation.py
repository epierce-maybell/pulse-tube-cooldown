"""Test physical invariants of a full cooldown simulation."""

import pytest

pytestmark = pytest.mark.skip(reason="waiting on config file for an actual simulation")


@pytest.fixture
def config_path():
    ...


@pytest.fixture
def model(config_path):
    ...


@pytest.fixture
def result(model):
    ...


@pytest.fixture
def passive_model(model):
    ...


def test_cooldown_conserves_energy(result):
    ...


def test_uniform_temperature_stays_put(passive_model):
    ...


def test_total_entropy_never_decreases(passive_model):
    ...


def test_temperatures_stay_within_initial_and_boundary_extremes(passive_model):
    ...


def test_increasing_thermal_mass_delays_cooldown(model):
    ...


def test_decreasing_cooling_path_conductance_delays_cooldown(model):
    ...


def test_scaling_all_thermal_storage_rescales_time(model):
    ...
