"""The `run-ptc` command, which simulates a configured cooldown and saves its Result."""

import argparse
import sys
import tomllib
from collections.abc import Sequence

import pulse_tube_cooldown as ptc
from pulse_tube_cooldown.core.errors import PTCError

__all__ = ["main"]


def main(argv: Sequence[str] | None = None) -> int:
    """Runs `run-ptc --config path_to_config --output path_to_output`.

    Args:
        argv: Command-line arguments, defaulting to `sys.argv[1:]`.

    Returns:
        The exit status: 0 when the run is saved, or 1 when the inputs fail or the solver fails. A solver failure still
        saves the run.
    """
    parser = argparse.ArgumentParser(prog="run-ptc", description="Simulates a cooldown and saves its Result as netCDF.")
    parser.add_argument("--config", required=True, help="path to the TOML configuration file")
    parser.add_argument("--output", required=True, help="path of the netCDF file to write")
    args = parser.parse_args(argv)

    try:
        model = ptc.configure_model(args.config)
        result = ptc.simulate_cooldown(model)
        ptc.save_output(result, args.output)
    except (PTCError, OSError, tomllib.TOMLDecodeError) as error:
        print(f"run-ptc: error: {error}", file=sys.stderr)
        return 1

    if result.status == "solver_failure":
        print(f"run-ptc: error: the solver failed: {result.data.attrs['message']}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
