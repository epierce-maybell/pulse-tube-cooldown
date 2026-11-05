"""A configured cooler builds its regenerator nodes and cools its stages."""

from pathlib import Path

import numpy as np

import pulse_tube_cooldown as ptc


def test_configured_cooler_cools_its_stages():
    model = ptc.configure_model(Path(__file__).parent / "fixtures/cooler.toml")
    result = ptc.simulate_cooldown(model)
    temperature = result.data.temperature

    assert result.status == "max_time"
    assert "PT stage 2 Gd2O2S 8" in temperature.node.values
    assert temperature.sel(node="Stage 1").values[-1] < 290.0
    assert temperature.sel(node="Stage 2").values[-1] < temperature.sel(node="Stage 1").values[-1]
    assert np.all(temperature.sel(node="Warm flange").values == 295.0)
    assert np.all(np.diff(temperature.sel(node="PT stage 2 Gd2O2S 8").values) <= 0.0)
