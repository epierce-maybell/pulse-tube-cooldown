"""A configured switch carries heat early and shuts off without changing phase."""

from pathlib import Path

import numpy as np

import pulse_tube_cooldown as ptc


def test_configured_switch_shuts_off_as_plate_cools():
    config = Path(__file__).parent / "fixtures/gas_gap_switch.toml"
    model = ptc.configure_model(config)
    result = ptc.simulate_cooldown(model)
    data = result.data
    plate = data.temperature.sel(node="Plate").values
    switch_flow = data.heat_flow.sel(link="Switch", node="Plate").values
    support_flow = data.heat_flow.sel(link="Support", node="Plate").values
    conductance = -switch_flow / (plate - 6.0)

    assert result.status == "max_time"
    assert set(data.phase.values) == {"cooldown"}
    assert switch_flow[0] < -1.0
    assert abs(switch_flow[0]) > abs(support_flow[0])
    assert abs(plate[-1] - 6.0) < 1e-6
    assert conductance[-1] < 1e-6 * conductance[0]
    assert abs(switch_flow[-1]) < 1e-3 * abs(support_flow[-1])
    np.testing.assert_allclose(data.heat_flow.sel(link="Switch").sum(dim="node"), 0.0, atol=1e-12)
