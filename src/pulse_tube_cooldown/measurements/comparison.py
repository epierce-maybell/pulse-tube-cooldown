"""Log-temperature residuals between simulations and measured cooldowns."""

import numpy as np
import xarray as xr

__all__ = ["log_temperature_ratio", "log_temperature_residual"]


def log_temperature_ratio(simulated: xr.DataArray, measured: xr.Dataset) -> xr.DataArray:
    """Returns ln(T_simulated / T_measured) at each measured sample.

    Args:
        simulated: Temperature(time, node) in K with time in s since compressor start, covering the measured nodes.
            Temperatures are interpolated linearly in log temperature at each measured sample.
        measured: A temperature(time, node) dataset in K.

    Returns:
        The ratio(time, node) on the measured samples. It is NaN where the reading is NaN or the sample lies outside
        the simulated time.
    """
    log_simulated = np.column_stack(
        [
            np.interp(
                measured.time.values,
                simulated.time.values,
                np.log(simulated.sel(node=node).values),
                left=np.nan,
                right=np.nan,
            )
            for node in measured.node.values
        ]
    )
    log_measured = np.log(measured.temperature.transpose("time", "node").values)
    return xr.DataArray(
        log_simulated - log_measured, coords={"time": measured.time, "node": measured.node}, dims=("time", "node")
    )


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

    scored = measured.isel(node=measured.scored.values)
    ratio = log_temperature_ratio(simulated, scored)
    log_ratios = {
        str(node): ratio.sel(node=node).values[~np.isnan(scored.temperature.sel(node=node).values)]
        for node in scored.node.values
    }

    def rms(values: np.ndarray) -> float:
        return float(np.sqrt(np.mean(values**2)))

    return rms(np.concatenate(list(log_ratios.values()))), {node: rms(ratio) for node, ratio in log_ratios.items()}
