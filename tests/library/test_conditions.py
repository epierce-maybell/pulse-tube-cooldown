"""Tests for conditions that end phases and procedures."""

import numpy as np
import pytest

from pulse_tube_cooldown.core.errors import ConfigurationError
from pulse_tube_cooldown.core.types import Graph, Material, Node
from pulse_tube_cooldown.library.conditions import CONDITIONS, condensing_begins, gas_pushing_begins

_TOY = Material(
    name="Toy",
    density=1.0,
    specific_heat=lambda t: t,
    thermal_conductivity=lambda t: t,
    conductivity_integral=lambda t: t**2 / 2.0,
)

_GRAPH = Graph(nodes=tuple(Node(name, ((_TOY, 1.0),)) for name in ("PRP", "CFP", "MXP")), links=())


def test_registers_gas_pushing_begins():
    assert CONDITIONS["gas_pushing_begins"] is gas_pushing_begins


@pytest.mark.parametrize(
    ("cfp_temperature", "sign"),
    [(120.0, 1.0), (70.0, 0.0), (40.0, -1.0)],
)
def test_gas_pushing_begins_changes_sign_as_the_node_cools_through_its_threshold(cfp_temperature, sign):
    condition = gas_pushing_begins(_GRAPH, node="CFP", threshold=70.0)

    assert np.sign(condition(0.0, np.array([300.0, cfp_temperature, 300.0]))) == sign


def test_gas_pushing_begins_ignores_other_nodes():
    condition = gas_pushing_begins(_GRAPH, node="CFP", threshold=70.0)

    assert condition(0.0, np.array([1.0, 70.0, 500.0])) == 0.0


def test_gas_pushing_begins_rejects_an_unknown_node_when_built():
    with pytest.raises(ConfigurationError, match="STP"):
        gas_pushing_begins(_GRAPH, node="STP", threshold=70.0)


def test_registers_condensing_begins():
    assert CONDITIONS["condensing_begins"] is condensing_begins


@pytest.mark.parametrize(
    ("mxp_temperature", "sign"),
    [(40.0, 1.0), (10.0, 0.0), (4.0, -1.0)],
)
def test_condensing_begins_changes_sign_as_the_node_cools_through_its_threshold(mxp_temperature, sign):
    condition = condensing_begins(_GRAPH, node="MXP", threshold=10.0)

    assert np.sign(condition(0.0, np.array([300.0, 300.0, mxp_temperature]))) == sign


def test_condensing_begins_ignores_other_nodes():
    condition = condensing_begins(_GRAPH, node="MXP", threshold=10.0)

    assert condition(0.0, np.array([1.0, 500.0, 10.0])) == 0.0


def test_condensing_begins_rejects_an_unknown_node_when_built():
    with pytest.raises(ConfigurationError, match="STP"):
        condensing_begins(_GRAPH, node="STP", threshold=10.0)
