"""Tests for the engine's steady-state solve."""

import dataclasses

import numpy as np
import pytest
from numpy.typing import NDArray

from pulse_tube_cooldown.core.errors import SimulationError
from pulse_tube_cooldown.core.types import EnergyBalanceTerms, Graph, Material, Node
from pulse_tube_cooldown.engine.steady_state import solve_steady_state

_MATERIAL = Material(
    name="Constant",
    density=1.0,
    specific_heat=lambda t: 1.0 + 0.0 * t,
    specific_heat_integral=lambda t: 1.0 * t,
    thermal_conductivity=None,
    conductivity_integral=None,
)


def _nodes(*names: str) -> tuple[Node, ...]:
    return tuple(Node(name, ((_MATERIAL, 1.0),)) for name in names)


@dataclasses.dataclass(frozen=True, slots=True)
class _Conductance:
    """Constant conductance G from the first node to the second."""

    name: str
    nodes: tuple[str, str]
    conductance: float

    def calc_energy_balance_terms(self, temperature: NDArray[np.float64]) -> EnergyBalanceTerms:
        flow = self.conductance * (temperature[1] - temperature[0])
        return EnergyBalanceTerms(heat_flow=np.array([flow, -flow]))


@dataclasses.dataclass(frozen=True, slots=True)
class _ProportionalConductance:
    """Conductance g T from the first node to the second, carrying g (T_2^2 - T_1^2) / 2 into the first."""

    name: str
    nodes: tuple[str, str]
    coefficient: float

    def calc_energy_balance_terms(self, temperature: NDArray[np.float64]) -> EnergyBalanceTerms:
        flow = self.coefficient * (temperature[1] ** 2 - temperature[0] ** 2) / 2.0
        return EnergyBalanceTerms(heat_flow=np.array([flow, -flow]))


@dataclasses.dataclass(frozen=True, slots=True)
class _CountedConductance:
    """Constant conductance G from the first node to the second, recording each evaluation in `calls`."""

    name: str
    nodes: tuple[str, str]
    conductance: float
    calls: list[None] = dataclasses.field(default_factory=list)

    def calc_energy_balance_terms(self, temperature: NDArray[np.float64]) -> EnergyBalanceTerms:
        self.calls.append(None)
        flow = self.conductance * (temperature[1] - temperature[0])
        return EnergyBalanceTerms(heat_flow=np.array([flow, -flow]))


def test_heater_power_crosses_a_chain_of_conductances():
    graph = Graph(
        nodes=_nodes("Bath", "Plate", "Sample"),
        links=(_Conductance("Rod", ("Bath", "Plate"), 2.0), _Conductance("Wire", ("Plate", "Sample"), 0.5)),
    )

    temperature = solve_steady_state(
        graph, np.array([4.0, 10.0, 10.0]), fixed_boundary_nodes=(0,), heat_input={"Sample": 1.0}, tolerance=1e-12
    )

    # 1 W crosses both links: the plate is 1 / 2 K above the bath, the sample 1 / 0.5 K above the plate.
    np.testing.assert_allclose(temperature, [4.0, 4.5, 6.5], rtol=1e-8)


def test_nonlinear_link_from_a_distant_guess():
    graph = Graph(nodes=_nodes("Bath", "Stage"), links=(_ProportionalConductance("Strap", ("Stage", "Bath"), 0.1),))

    temperature = solve_steady_state(
        graph, np.array([4.0, 300.0]), fixed_boundary_nodes=(0,), heat_input={"Stage": 1.0}, tolerance=1e-12
    )

    # 1 W = 0.1 (T^2 - 4^2) / 2, so T = 6 K.
    np.testing.assert_allclose(temperature, [4.0, 6.0], rtol=1e-8)


def test_without_heat_input_free_nodes_settle_at_the_boundary_temperature():
    graph = Graph(
        nodes=_nodes("Room", "Shield", "Stage"),
        links=(_Conductance("Post", ("Room", "Shield"), 1.0), _Conductance("Rod", ("Shield", "Stage"), 3.0)),
    )

    temperature = solve_steady_state(graph, np.array([295.0, 40.0, 4.0]), fixed_boundary_nodes=(0,))

    np.testing.assert_allclose(temperature, [295.0, 295.0, 295.0], rtol=1e-8)


def test_no_steady_state_raises():
    graph = Graph(nodes=_nodes("Bath", "Isolated"), links=())

    with pytest.raises(SimulationError):
        solve_steady_state(graph, np.array([4.0, 10.0]), fixed_boundary_nodes=(0,), heat_input={"Isolated": 1.0})


def test_state_below_the_floor_raises_promptly():
    link = _CountedConductance("Rod", ("Bath", "Stage"), 1.0)
    graph = Graph(nodes=_nodes("Bath", "Stage"), links=(link,))

    # 1 W/K (3 - T) = 2 W gives T = 1 K, below the 2.1 K floor.
    with pytest.raises(SimulationError, match="Stage"):
        solve_steady_state(graph, np.array([3.0, 10.0]), fixed_boundary_nodes=(0,), heat_input={"Stage": -2.0})

    # A tenth of the heat-flow evaluations a single free node is allowed by the default 500 steps.
    assert len(link.calls) <= 100


def test_lower_floor_admits_a_colder_state():
    graph = Graph(nodes=_nodes("Bath", "Stage"), links=(_Conductance("Rod", ("Bath", "Stage"), 1.0),))

    temperature = solve_steady_state(
        graph,
        np.array([3.0, 10.0]),
        fixed_boundary_nodes=(0,),
        heat_input={"Stage": -2.0},
        temperature_floor=0.5,
        tolerance=1e-12,
    )

    np.testing.assert_allclose(temperature, [3.0, 1.0], rtol=1e-8)
