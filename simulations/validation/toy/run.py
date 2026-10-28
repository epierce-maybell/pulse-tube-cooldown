"""Runs the toy cooldown and saves its Result beside the configuration."""

from pathlib import Path

import pulse_tube_cooldown as ptc

here = Path(__file__).parent
model = ptc.configure_model(here / "toy.toml")
result = ptc.simulate_cooldown(model)
ptc.save_output(result, here / "toy.nc")
