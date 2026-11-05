"""A bare PT1.5 head on one compressor cools from room temperature to a cold steady state."""

import numpy as np
from scipy.integrate import solve_ivp

from pulse_tube_cooldown.library.gases import helium
from pulse_tube_cooldown.library.heads import load_head_model
from pulse_tube_cooldown.library.links.cooler import (
    Compressor,
    Head,
    calc_charge,
    calc_gas_mass,
    calc_heat_flow,
    calc_mass_swing,
    calc_mean_pressure,
    calc_pressure_swing,
)
from pulse_tube_cooldown.library.materials.load import load_material

_ROOM_TEMPERATURE = 295.0
_WARM_VOLUME = 19.6e-3  # m^3, data/compressors.md

# Placeholder packing coefficients and valve resistances pending task 2.3's fit.
_HEAD = Head(model=load_head_model("PT1.5"), packing=(0.15, 0.05), resistance=(1e10, 3e10), frequency=1.4)


def _stage_heat_capacity(stage_1, stage_2):
    """Heat capacities in J/K of the test stand's stage 1 and stage 2, from the masses in data/test_stand.md."""
    copper, aluminum, steel = (load_material(name) for name in ("OFHC Cu RRR100", "Al 6061-T6", "304 SS"))
    first, second = np.array([stage_1]), np.array([stage_2])
    return (
        (
            1.5838 * copper.specific_heat(first)
            + 0.4536 * aluminum.specific_heat(first)
            + 0.0709 * steel.specific_heat(first)
        ).item(),
        (0.4218 * steel.specific_heat(second) + 0.0397 * copper.specific_heat(second)).item(),
    )


def test_bare_head_reaches_a_cold_steady_state():
    model = _HEAD.model
    layers = [layer for section in model.regenerators for layer in section.layers for _ in range(layer.node_count)]
    stage_1 = 1 + sum(layer.node_count for layer in model.regenerators[0].layers)
    stage_2 = len(layers) + 2
    nodes = [index for index in range(1, stage_2) if index != stage_1]
    compressors = [Compressor(displacement_rate=5.72e-3, relief_swing=1.93e6)]
    charge = calc_charge([model], warm_volume=_WARM_VOLUME, fill_pressure=1.70e6, fill_temperature=_ROOM_TEMPERATURE)

    def heat_flow(temperature):
        mean_pressure = calc_mean_pressure(
            lambda pressure: calc_gas_mass(model, temperature, pressure),
            charge,
            warm_volume=_WARM_VOLUME,
            warm_temperature=_ROOM_TEMPERATURE,
        )
        swing = calc_pressure_swing(
            lambda s: _HEAD.frequency * calc_mass_swing(_HEAD, temperature, mean_pressure, s).sum(),
            mean_pressure,
            _ROOM_TEMPERATURE,
            compressors=compressors,
        )
        return calc_heat_flow(_HEAD, temperature, mean_pressure, swing), mean_pressure

    def rate(_, cold):
        temperature = np.concatenate(([_ROOM_TEMPERATURE], cold))
        heat, mean_pressure = heat_flow(temperature)
        capacity = np.empty_like(temperature)
        for index, layer in zip(nodes, layers, strict=True):
            t = temperature[index]
            matrix = layer.mass / layer.node_count * layer.material.specific_heat(np.array([t])).item()
            gas = (
                layer.pore_volume
                / layer.node_count
                * float(helium.density(t, mean_pressure) * helium.isobaric_specific_heat(t, mean_pressure))
            )
            capacity[index] = matrix + gas
        capacity[stage_1], capacity[stage_2] = _stage_heat_capacity(temperature[stage_1], temperature[stage_2])
        return heat[1:] / capacity[1:]

    end = 2e4  # s
    solution = solve_ivp(
        rate, (0.0, end), np.full(stage_2, _ROOM_TEMPERATURE), method="BDF", rtol=1e-4, atol=1e-3, t_eval=[end / 2, end]
    )
    final = np.concatenate(([_ROOM_TEMPERATURE], solution.y[:, -1]))

    assert solution.success
    np.testing.assert_allclose(solution.y[:, 0], solution.y[:, 1], rtol=1e-4)
    assert np.max(np.abs(heat_flow(final)[0][1:])) < 1e-3
    assert 20.0 < final[stage_1] < 80.0
    assert final[stage_2] < 12.0
