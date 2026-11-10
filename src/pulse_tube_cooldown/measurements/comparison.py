"""Log-temperature residuals between simulations and measured cooldowns."""

import numpy as np
import xarray as xr

__all__ = ["log_temperature_residual"]


def log_temperature_residual(simulated: xr.DataArray, measured: xr.Dataset) -> tuple[float, dict[str, float]]:
    """Returns pooled and per-node RMS values of ln(T_simulated / T_measured).

    Args:
        simulated: Temperature(time, node) in K with time in s since compressor start. Temperatures are interpolated
            linearly in log temperature at each valid measured sample.
        measured: A temperature(time, node) dataset in K with a scored(node) mask. NaN readings are omitted.

    Raises:
        ValueError: The simulated run ends before the measured window.
    """
    if simulated.time[-1] < measured.time[-1]:
        raise ValueError(
            f"Run ends at {float(simulated.time[-1])} s, before the window's {float(measured.time[-1])} s."
        )

    log_ratios = {}
    for node in measured.node.values[measured.scored.values]:
        temperature = measured.temperature.sel(node=node).values
        valid = ~np.isnan(temperature)
        log_simulated = np.interp(measured.time.values[valid], simulated.time.values, np.log(simulated.sel(node=node)))
        log_ratios[str(node)] = log_simulated - np.log(temperature[valid])

    def rms(values: np.ndarray) -> float:
        return float(np.sqrt(np.mean(values**2)))

    return rms(np.concatenate(list(log_ratios.values()))), {node: rms(ratio) for node, ratio in log_ratios.items()}
