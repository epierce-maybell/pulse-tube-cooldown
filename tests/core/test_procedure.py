"""Tests for procedure validation and phase graphs."""

import dataclasses

import numpy as np
import pytest
from numpy.typing import NDArray

from pulse_tube_cooldown.core.errors import ConfigurationError, PTCError
from pulse_tube_cooldown.core.types import EnergyBalanceTerms, Material, Model, Node, Phase, Procedure

_TOY = Material(
    name="Toy",
    density=1.0,
    specific_heat=lambda t: t,
    specific_heat_integral=lambda t: t**2 / 2.0,
    thermal_conductivity=lambda t: t,
    conductivity_integral=lambda t: t**2 / 2.0,
)


@dataclasses.dataclass(frozen=True, slots=True)
class _StubLink:
    """Link whose flow is irrelevant to procedure structure."""

    name: str
    nodes: tuple[str, ...]

    def calc_energy_balance_terms(self, temperature: NDArray[np.float64]) -> EnergyBalanceTerms:
        return EnergyBalanceTerms(heat_flow=np.zeros_like(temperature))


_NODES = tuple(Node(name, ((_TOY, 1.0),)) for name in ("PRP", "RGP", "Can"))
_EXCHANGE_GAS = _StubLink("exchange gas", ("PRP", "RGP", "Can"))
_SWITCH = _StubLink("switch", ("PRP", "RGP"))


def _cutoff(time: float, temperature: NDArray[np.float64]) -> float:
    return float(temperature[0] - 70.0)


def _two_phase_procedure() -> Procedure:
    return Procedure(
        phases=(Phase("precool", (_EXCHANGE_GAS,), _cutoff), Phase("cooldown", (_SWITCH,), None)),
        stop_condition=_cutoff,
    )


def _model(procedure: Procedure | None = None, nodes: tuple[Node, ...] = _NODES) -> Model:
    return Model({}, nodes, procedure or _two_phase_procedure(), np.full(len(nodes), 300.0), (2,), {})


@pytest.mark.parametrize(
    ("phases", "mentions"),
    [
        ((), []),
        ((Phase("cool", (), _cutoff), Phase("cool", (), None)), ["cool"]),
        ((Phase("precool", (), None), Phase("cooldown", (), None)), ["precool"]),
        ((Phase("precool", (), _cutoff), Phase("cooldown", (), _cutoff)), ["cooldown"]),
        (
            (Phase("precool", (_SWITCH,), _cutoff), Phase("cooldown", (_StubLink("switch", ("RGP", "Can")),), None)),
            ["switch"],
        ),
        (
            (Phase("precool", (_SWITCH,), _cutoff), Phase("cooldown", (_StubLink("switch", ("RGP", "PRP")),), None)),
            ["switch"],
        ),
    ],
    ids=[
        "no phases",
        "repeated phase name",
        "non-final phase without until",
        "final phase with until",
        "link name with different nodes",
        "link name with reordered nodes",
    ],
)
def test_invalid_procedure_raises_value_error_naming_offenders(phases, mentions):
    with pytest.raises(ConfigurationError) as error:
        Procedure(phases=phases, stop_condition=None)

    assert isinstance(error.value, ValueError)
    assert isinstance(error.value, PTCError)
    for mention in mentions:
        assert repr(mention) in str(error.value)


def test_link_with_reordered_nodes_reports_both_orders():
    reversed_switch = _StubLink("switch", ("RGP", "PRP"))

    with pytest.raises(ConfigurationError) as error:
        Procedure(
            phases=(Phase("precool", (_SWITCH,), _cutoff), Phase("cooldown", (reversed_switch,), None)),
            stop_condition=None,
        )

    assert "order" in str(error.value)
    assert repr(("PRP", "RGP")) in str(error.value)
    assert repr(("RGP", "PRP")) in str(error.value)


def test_model_builds_one_graph_per_phase_from_shared_nodes():
    model = _model()

    assert len(model.graphs) == 2
    for graph in model.graphs:
        assert [node.name for node in graph.nodes] == ["PRP", "RGP", "Can"]
    assert [link.name for link in model.graphs[0].links] == ["exchange gas"]
    assert model.graphs[0].connectivity == ((0, 1, 2),)
    assert [link.name for link in model.graphs[1].links] == ["switch"]
    assert model.graphs[1].connectivity == ((0, 1),)


@pytest.mark.parametrize(
    ("links", "mentions"),
    [
        ((_StubLink("strut", ("PRP", "Still")),), ["cooldown", "strut", "Still"]),
        ((_SWITCH, _SWITCH), ["cooldown", "switch"]),
    ],
    ids=["unknown node", "duplicate link"],
)
def test_invalid_phase_graph_raises_value_error_naming_phase_and_offenders(links, mentions):
    procedure = Procedure(
        phases=(Phase("precool", (_EXCHANGE_GAS,), _cutoff), Phase("cooldown", links, None)),
        stop_condition=None,
    )

    with pytest.raises(ConfigurationError) as error:
        _model(procedure)

    assert isinstance(error.value, ValueError)
    assert isinstance(error.value, PTCError)
    for mention in mentions:
        assert repr(mention) in str(error.value)
