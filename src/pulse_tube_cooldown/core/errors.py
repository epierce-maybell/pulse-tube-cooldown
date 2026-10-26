"""Exceptions raised by the simulator.

Every exception derives from `PTCError`, and specific exceptions also derive from the matching builtin.
"""

__all__ = ["ConfigurationError", "PTCError", "SimulationError"]


class PTCError(Exception):
    """Base class for all simulator exceptions."""


class ConfigurationError(PTCError, ValueError):
    """Inputs do not describe a valid model, such as a malformed graph or an unknown name."""


class SimulationError(PTCError, RuntimeError):
    """A valid model fails to simulate, such as when the solver fails."""
