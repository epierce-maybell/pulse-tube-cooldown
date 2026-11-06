"""Checks of stage coupling after the internal matrix has re-equilibrated."""

from dataclasses import dataclass

import numpy as np
import pytest

from pulse_tube_cooldown.core.errors import SimulationError
from pulse_tube_cooldown.core.types import EnergyBalanceTerms, Graph, Material, Node
from simulations.validation.test_stand_steady_state.coupling import analyze_coupling


@dataclass(frozen=True)
class Conductance:
    """Constant thermal conductance between two nodes."""

    name: str
    nodes: tuple[str, str]
    value: float

    def calc_energy_balance_terms(self, temperature):
        """Returns the heat flow into the two nodes."""
        flow = self.value * (temperature[1] - temperature[0])
        return EnergyBalanceTerms(heat_flow=np.array([flow, -flow]))


def graph():
    """Builds two stages sharing a matrix node linked to a bath."""
    material = Material(
        name="Constant",
        density=1.0,
        specific_heat=lambda t: np.ones_like(t),
        specific_heat_integral=lambda t: np.asarray(t, dtype=float),
        thermal_conductivity=None,
        conductivity_integral=None,
    )
    return Graph(
        nodes=tuple(Node(name, ((material, 1.0),)) for name in ("Bath", "First", "Second", "Matrix")),
        links=(
            Conductance("First matrix", ("First", "Matrix"), 2.0),
            Conductance("Second matrix", ("Second", "Matrix"), 3.0),
            Conductance("Matrix bath", ("Matrix", "Bath"), 5.0),
        ),
    )


def test_stage_derivatives_include_matrix_response():
    """Internal equilibration produces the analytic cross-coupling."""
    # Matrix temperature = (2 T_first + 3 T_second + 5 T_bath) / 10.
    result = analyze_coupling(
        graph(), np.array([4.0, 14.0, 24.0, 12.0]), stages=("First", "Second"), fixed_boundary_nodes=(0,)
    )
    np.testing.assert_allclose(result.stage_derivative, [[-1.6, 0.6], [0.6, -2.1]], atol=1e-9)
    np.testing.assert_allclose(result.stage_heat_flow, [-4.0, -36.0], atol=1e-9)
    assert result.matrix_growth_rate == pytest.approx(-10.0)
    assert result.free_node_growth_rate < 0.0


def test_unbalanced_matrix_is_rejected_before_reporting_coupling():
    """A non-equilibrium matrix cannot supply a load-map derivative."""
    with pytest.raises(SimulationError, match="Matrix"):
        analyze_coupling(
            graph(), np.array([4.0, 14.0, 24.0, 12.1]), stages=("First", "Second"), fixed_boundary_nodes=(0,)
        )
