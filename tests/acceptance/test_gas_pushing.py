"""Test conservation, inactive limits, and equivalent representations for gas pushing."""

from dataclasses import replace

import numpy as np
import pytest

from pulse_tube_cooldown.library.links.gas_pushing import GasPushing, GasVolume


@pytest.fixture
def link():
    """Return an ordered gas path with three plate volumes and a warm manifold."""
    return GasPushing(
        name="Gas pushing",
        nodes=("Manifold", "STP", "ICP", "MXP"),
        volumes=(
            GasVolume(nodes=("STP",), volume=0.001),
            GasVolume(nodes=("ICP",), volume=0.0006),
            GasVolume(nodes=("MXP",), volume=0.0004),
        ),
        warm_volume=0.003,
        buffer_volume=0.012,
        pump_speed=0.003,
        flow_conductance=1.0e-6,
        fill_pressure=80000.0,
        wave_strength=40000.0,
        refill_pressure=20000.0,
        buffer_pressure_cap=125000.0,
        depressurization_time=1.0,
    )


@pytest.fixture(params=["separate", "paired", "overlapping"])
def volume_link(request, link):
    """Return separate plate volumes, a two-node volume, or volumes sharing a plate."""
    volumes = {
        "separate": link.volumes,
        "paired": (
            GasVolume(nodes=("STP", "ICP"), volume=0.0016),
            GasVolume(nodes=("MXP",), volume=0.0004),
        ),
        "overlapping": (
            GasVolume(nodes=("STP", "ICP"), volume=0.001),
            GasVolume(nodes=("ICP", "MXP"), volume=0.001),
        ),
    }
    return replace(link, volumes=volumes[request.param])


@pytest.fixture(params=[(125000.0, 20000.0), (95000.0, 20000.0), (95000.0, 2000.0)])
def cycle_link(request, link):
    """Return wave-limited cycles or cap-limited cycles with two refill thresholds, all in Pa."""
    cap, refill = request.param
    return replace(link, buffer_pressure_cap=cap, refill_pressure=refill)


@pytest.mark.parametrize(
    "temperature",
    [(300.0, 150.0, 80.0, 50.0), (50.0, 150.0, 80.0, 300.0), (100.0, 300.0, 50.0, 80.0)],
)
def test_gas_transport_conserves_energy_including_manifold(volume_link, temperature):
    """Verify that modeled gas transport balances heat across the manifold and all plate nodes."""
    terms = volume_link.calc_energy_balance_terms(np.array(temperature))

    assert terms.heat_flow.sum() == pytest.approx(0.0, abs=1e-12)


@pytest.mark.parametrize("temperature", [20.0, 100.0, 300.0])
def test_uniform_temperatures_carry_no_heat(volume_link, temperature):
    """Verify zero heat transport when every node shares one temperature."""
    terms = volume_link.calc_energy_balance_terms(np.full(4, temperature))

    np.testing.assert_array_equal(terms.heat_flow, [0.0, 0.0, 0.0, 0.0])


@pytest.mark.parametrize(
    "temperature",
    [(300.0, 150.0, 80.0, 50.0), (50.0, 150.0, 80.0, 300.0), (100.0, 300.0, 50.0, 80.0)],
)
def test_hottest_node_loses_heat_and_coldest_gains_heat_for_single_node_volumes(cycle_link, temperature):
    """Verify passive heat exchange at the temperature extremes when each volume touches one node."""
    temperature = np.array(temperature)

    terms = cycle_link.calc_energy_balance_terms(temperature)

    assert terms.heat_flow[np.argmax(temperature)] < 0.0
    assert terms.heat_flow[np.argmin(temperature)] > 0.0


@pytest.mark.parametrize("parameter", ["pump_speed", "flow_conductance", "wave_strength"])
def test_no_circulation_carries_no_heat(cycle_link, parameter):
    """Verify zero heat transport with a stopped pump, closed path, or zero pressure wave."""
    inactive_link = replace(cycle_link, **{parameter: 0.0})

    terms = inactive_link.calc_energy_balance_terms(np.array([300.0, 150.0, 80.0, 50.0]))

    np.testing.assert_array_equal(terms.heat_flow, [0.0, 0.0, 0.0, 0.0])


@pytest.mark.parametrize("parameter", ["pump_speed", "flow_conductance"])
def test_heat_transport_tends_continuously_to_zero_with_circulation(cycle_link, parameter):
    """Verify diminishing nonzero transport as pumping speed or path conductance approaches zero."""
    initial = getattr(cycle_link, parameter)
    temperature = np.array([300.0, 150.0, 80.0, 50.0])

    transport = np.array(
        [
            np.abs(
                replace(cycle_link, **{parameter: initial * factor}).calc_energy_balance_terms(temperature).heat_flow
            ).sum()
            for factor in [0.01, 0.0001, 0.000001]
        ]
    )

    assert np.all(transport > 0.0)
    assert np.all(transport[1:] < 0.05 * transport[:-1])


@pytest.mark.parametrize("volume_index", [0, 1, 2])
def test_partitioning_an_isothermal_volume_preserves_heat_flows(cycle_link, volume_index):
    """Verify that splitting a gas volume into consecutive isothermal pieces preserves node heat flows."""
    volume = cycle_link.volumes[volume_index]
    partitioned_link = replace(
        cycle_link,
        volumes=(
            *cycle_link.volumes[:volume_index],
            replace(volume, volume=volume.volume / 4),
            replace(volume, volume=3 * volume.volume / 4),
            *cycle_link.volumes[volume_index + 1 :],
        ),
    )
    temperature = np.array([300.0, 150.0, 80.0, 50.0])

    original = cycle_link.calc_energy_balance_terms(temperature)
    partitioned = partitioned_link.calc_energy_balance_terms(temperature)

    np.testing.assert_allclose(partitioned.heat_flow, original.heat_flow, rtol=1e-12, atol=0.0)


@pytest.mark.parametrize("time_scale", [0.5, 2.0, 10.0])
def test_stretching_cycle_time_reduces_average_heat_flow_inversely(cycle_link, time_scale):
    """Verify inverse average-power scaling when every part of the cycle takes proportionally longer."""
    stretched_link = replace(
        cycle_link,
        pump_speed=cycle_link.pump_speed / time_scale,
        flow_conductance=cycle_link.flow_conductance / time_scale,
        depressurization_time=cycle_link.depressurization_time * time_scale,
    )
    temperature = np.array([300.0, 150.0, 80.0, 50.0])

    original = cycle_link.calc_energy_balance_terms(temperature)
    stretched = stretched_link.calc_energy_balance_terms(temperature)

    np.testing.assert_allclose(stretched.heat_flow, original.heat_flow / time_scale, rtol=1e-12, atol=0.0)


def test_loosening_buffer_cap_preserves_wave_limited_heat_flows(link):
    """Verify that the wave threshold controls transport when the buffer cap is sufficiently high."""
    wave_limited_link = replace(link, wave_strength=20000.0, buffer_pressure_cap=1.0e6)
    relaxed_cap_link = replace(wave_limited_link, buffer_pressure_cap=2.0e6)
    changed_wave_link = replace(wave_limited_link, wave_strength=40000.0)
    temperature = np.array([300.0, 150.0, 80.0, 50.0])

    original = wave_limited_link.calc_energy_balance_terms(temperature)
    relaxed_cap = relaxed_cap_link.calc_energy_balance_terms(temperature)
    changed_wave = changed_wave_link.calc_energy_balance_terms(temperature)

    np.testing.assert_allclose(relaxed_cap.heat_flow, original.heat_flow, rtol=1e-12, atol=0.0)
    assert not np.allclose(changed_wave.heat_flow, original.heat_flow, rtol=1e-6, atol=0.0)


@pytest.mark.parametrize("refill_pressure", [20000.0, 2000.0])
def test_loosening_wave_threshold_preserves_cap_limited_heat_flows(link, refill_pressure):
    """Verify that the buffer cap controls transport when the wave threshold is sufficiently high."""
    cap_limited_link = replace(
        link, buffer_pressure_cap=95000.0, wave_strength=200000.0, refill_pressure=refill_pressure
    )
    relaxed_wave_link = replace(cap_limited_link, wave_strength=400000.0)
    changed_cap_link = replace(cap_limited_link, buffer_pressure_cap=105000.0)
    temperature = np.array([300.0, 150.0, 80.0, 50.0])

    original = cap_limited_link.calc_energy_balance_terms(temperature)
    relaxed_wave = relaxed_wave_link.calc_energy_balance_terms(temperature)
    changed_cap = changed_cap_link.calc_energy_balance_terms(temperature)

    np.testing.assert_allclose(relaxed_wave.heat_flow, original.heat_flow, rtol=1e-12, atol=0.0)
    assert not np.allclose(changed_cap.heat_flow, original.heat_flow, rtol=1e-6, atol=0.0)
