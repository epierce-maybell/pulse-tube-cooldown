"""Requirements for the gas-pushing link from the gas-pushing ADR."""

from dataclasses import replace

import numpy as np
import pytest
from scipy.constants import R

from pulse_tube_cooldown import simulate_cooldown
from pulse_tube_cooldown.config import configure_model
from pulse_tube_cooldown.library.links.gas_pushing import GasPushing, GasVolume


@pytest.fixture
def pushing():
    # At 300, 100, 50 K: A = B = 4e-5 / R mol/Pa, tau_p = 4 s,
    # tau_d = 2 s. A discharge lasting 2 ln(2) s halves the pressure difference,
    # so the still-side excursion is one quarter of the peak difference.
    return GasPushing(
        name="Gas pushing",
        nodes=("Warm manifold", "STP", "MXP"),
        volumes=(GasVolume(nodes=("STP",), volume=0.001), GasVolume(nodes=("MXP",), volume=0.001)),
        warm_volume=0.003,
        buffer_volume=0.012,
        pump_speed=0.003,
        flow_conductance=1e-5 / R,
        fill_pressure=80000.0,
        wave_strength=40000.0,
        refill_pressure=20000.0,
        buffer_pressure_cap=125000.0,
        depressurization_time=2.0 * np.log(2.0),
    )


def test_wave_limited_cycle_matches_hand_calculation(pushing):
    # Delta P = 10000 Pa; P_max = 75000 Pa; P_min = 65000 Pa.
    # Downstream mole capacities are 3e-5/R and 2e-5/R mol/Pa.
    # Heat crossing the two element boundaries is 0.015 f Delta P and
    # 0.0025 f Delta P W, respectively.
    frequency = 1.0 / (4.0 * np.log(15.0 / 13.0) + 2.0 * np.log(2.0))
    flow = pushing.calc_heat_flow(np.array([300.0, 100.0, 50.0]))

    np.testing.assert_allclose(flow, np.array([-0.015, 0.0125, 0.0025]) * frequency * 10000.0, rtol=1e-12)


@pytest.mark.parametrize(
    ("refill_pressure", "excursion", "maximum", "minimum"),
    [
        # P_min + 4 Delta P = 95000 Pa at the buffer cap, while
        # P_max = 70000 + Delta P / 2 during the averaged refill cycle.
        (20000.0, 50000.0 / 7.0, 515000.0 / 7.0, 465000.0 / 7.0),
        # Once Delta P exceeds the refill threshold, P_max = P_fill.
        (2000.0, 5000.0, 80000.0, 75000.0),
    ],
)
def test_buffer_cap_stops_pumping_before_wave_limit(pushing, refill_pressure, excursion, maximum, minimum):
    pushing = replace(pushing, buffer_pressure_cap=95000.0, refill_pressure=refill_pressure)
    frequency = 1.0 / (4.0 * np.log(maximum / minimum) + 2.0 * np.log(2.0))

    flow = pushing.calc_heat_flow(np.array([300.0, 100.0, 50.0]))

    np.testing.assert_allclose(flow, np.array([-0.015, 0.0125, 0.0025]) * frequency * excursion, rtol=1e-12)


def test_volume_uses_log_mean_temperature_and_splits_heat_equally(pushing):
    # The log mean of 200 ln(2) and 100 ln(2) K is exactly 100 K.
    pushing = replace(
        pushing,
        nodes=("Warm manifold", "STP", "ICP", "MXP"),
        volumes=(GasVolume(nodes=("STP", "ICP"), volume=0.001), GasVolume(nodes=("MXP",), volume=0.001)),
    )
    temperature = np.array([300.0, 200.0 * np.log(2.0), 100.0 * np.log(2.0), 50.0])
    frequency = 1.0 / (4.0 * np.log(15.0 / 13.0) + 2.0 * np.log(2.0))

    flow = pushing.calc_heat_flow(temperature)

    np.testing.assert_allclose(flow, np.array([-0.015, 0.00625, 0.00625, 0.0025]) * frequency * 10000.0, rtol=1e-12)


@pytest.mark.parametrize("temperature", [[300.0, 100.0, 50.0], [100.0, 300.0, 50.0], [50.0, 100.0, 300.0]])
def test_signed_element_flows_conserve_energy(pushing, temperature):
    flow = pushing.calc_heat_flow(np.array(temperature))

    assert flow.sum() == pytest.approx(0.0, abs=1e-12)
    assert np.sign(flow[0]) == np.sign(temperature[1] - temperature[0])
    assert np.sign(flow[-1]) == np.sign(temperature[-2] - temperature[-1])


def test_uniform_temperatures_carry_no_heat_including_two_node_volumes(pushing):
    pushing = replace(pushing, volumes=(GasVolume(nodes=("STP", "MXP"), volume=0.002),))

    flow = pushing.calc_heat_flow(np.full(3, 100.0))

    np.testing.assert_array_equal(flow, np.zeros(3))


@pytest.fixture
def configured_model(tmp_path):
    path = tmp_path / "pushing.toml"
    path.write_text(
        """
[[fridge.nodes]]
name = "Warm manifold"
composition = [["304 SS", 1.0]]
[[fridge.nodes]]
name = "STP"
composition = [["OFHC Cu RRR100", 1.0]]
[[fridge.nodes]]
name = "MXP"
composition = [["OFHC Cu RRR100", 1.0]]

[[fridge.links]]
name = "Gas pushing"
type = "gas_pushing"
nodes = ["Warm manifold", "STP", "MXP"]
warm_volume = 0.003
buffer_volume = 0.012
pump_speed = 0.003
flow_conductance = 1.2027235504272605e-6
fill_pressure = 80000.0
wave_strength = 40000.0
refill_pressure = 20000.0
buffer_pressure_cap = 125000.0
depressurization_time = 1.3862943611198906
volumes = [{nodes = ["STP"], volume = 0.001}, {nodes = ["MXP"], volume = 0.001}]

[[procedure.phases]]
name = "gas pushing"

[simulation]
solver = "RK45"
initial_temperature = [300.0, 100.0, 50.0]
fixed_boundary_nodes = ["Warm manifold"]
max_time = 1.0
save_every = 0.5
""",
        encoding="utf-8",
    )

    return configure_model(path)


def test_toml_configuration_loads_the_registered_link(configured_model):
    model = configured_model
    (link,) = model.graphs[0].links

    assert isinstance(link, GasPushing)
    assert link.nodes == ("Warm manifold", "STP", "MXP")
    assert model.fixed_boundary_nodes == (0,)
    frequency = 1.0 / (4.0 * np.log(15.0 / 13.0) + 2.0 * np.log(2.0))
    np.testing.assert_allclose(
        link.calc_heat_flow(model.initial_temperature),
        np.array([-0.015, 0.0125, 0.0025]) * frequency * 10000.0,
        rtol=1e-12,
    )


def test_shared_node_accumulates_heat_from_both_volume_elements(pushing):
    pushing = replace(
        pushing,
        nodes=("Warm manifold", "STP", "ICP", "MXP"),
        volumes=(GasVolume(nodes=("STP", "ICP"), volume=0.001), GasVolume(nodes=("ICP", "MXP"), volume=0.001)),
    )
    # The two log-mean gas temperatures remain 100 and 50 K.
    temperature = np.array([300.0, 200.0 * np.log(2.0), 100.0 * np.log(2.0), 50.0 * np.log(2.0)])
    frequency = 1.0 / (4.0 * np.log(15.0 / 13.0) + 2.0 * np.log(2.0))

    flow = pushing.calc_heat_flow(temperature)

    np.testing.assert_allclose(flow, np.array([-0.015, 0.00625, 0.0075, 0.00125]) * frequency * 10000.0, rtol=1e-12)


def test_zero_wave_strength_carries_no_heat(pushing):
    flow = replace(pushing, wave_strength=0.0).calc_heat_flow(np.array([300.0, 100.0, 50.0]))

    np.testing.assert_array_equal(flow, np.zeros(3))


def test_configured_run_preserves_warm_boundary_and_conserves_link_heat(configured_model):
    result = simulate_cooldown(configured_model)

    assert result.status == "max_time"
    temperature = result.data.temperature
    np.testing.assert_array_equal(temperature.sel(node="Warm manifold"), 300.0)
    assert float(temperature.sel(node="STP").isel(time=-1)) > 100.0
    assert float(temperature.sel(node="MXP").isel(time=-1)) > 50.0
    np.testing.assert_allclose(result.data.heat_flow.sum(dim="node"), 0.0, atol=1e-12)
