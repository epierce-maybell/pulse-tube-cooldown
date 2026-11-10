"""Five-coefficient, two-temperature cooling law for the unloaded stand.

P1 = A*(1 - B1/T1) - K*(1 - B2/T2)/T1
P2 = g2*(T2 - B2)

P1 and P2 are net cooling powers in W. A is in W, B1 and B2 in K,
K in W K, and g2 in W/K. The stage-1 cold slope and warm limit share A.
The stage-2 cold slope also sets its warm response.
"""

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray


@dataclass(frozen=True, slots=True)
class FiveCoefficientLaw:
    """Stage-1 scale, base, coupling, and stage-2 slope and base."""

    power_1: float
    base_1: float
    coupling: float
    slope_2: float
    base_2: float

    @classmethod
    def from_array(cls, values: NDArray) -> "FiveCoefficientLaw":
        """Reads A, B1, K, g2, B2 in the configuration's order."""
        return cls(*values)

    def power(self, temperature: NDArray) -> NDArray:
        """Returns net cooling powers in W, with stages on the last axis."""
        stage_1, stage_2 = temperature[..., 0], temperature[..., 1]
        first = self.power_1 * (1 - self.base_1 / stage_1)
        first -= self.coupling * (1 - self.base_2 / stage_2) / stage_1
        second = self.slope_2 * (stage_2 - self.base_2)
        return np.stack((first, second), axis=-1)


def calibrate_map(temperatures: NDArray, powers: NDArray) -> NDArray:
    """Fits five coefficients to net head cooling powers by linear least squares.

    P1 = a + b/T1 + c/(T1*T2); P2 = d*T2 + e. Converting these five
    linear coefficients gives the physical-coordinate form above.
    """
    first, second = temperatures.T
    a, b, c = np.linalg.lstsq(
        np.column_stack((np.ones_like(first), 1 / first, 1 / (first * second))), powers[:, 0], rcond=None
    )[0]
    d, e = np.linalg.lstsq(np.column_stack((second, np.ones_like(second))), powers[:, 1], rcond=None)[0]
    base_2 = -e / d
    coupling = c / base_2
    return np.array([a, -(b + coupling) / a, coupling, d, base_2])
