"""Saving simulation results."""

from os import PathLike

from pulse_tube_cooldown.core.types import Result

__all__ = ["save_output"]


def save_output(result: Result, path_to_output: str | PathLike[str]) -> None:
    """Writes a Result's data to a netCDF file.

    Args:
        result: The Result to save.
        path_to_output: Path of the netCDF file to write, replacing any existing file.
    """
    result.data.to_netcdf(path_to_output, engine="netcdf4")
