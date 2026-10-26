"""Tests for phases, procedures, models, and results."""

import dataclasses

import numpy as np
import pytest
import xarray as xr
from numpy.typing import NDArray

from pulse_tube_cooldown.core.errors import ConfigurationError, PTCError
from pulse_tube_cooldown.core.types import Model, Node, Phase, Procedure, Result


@dataclasses.dataclass(frozen=True, slots=True)
class _StubLink:
    """Link whose flow is irrelevant to procedure structure."""

    name: str
    nodes: tuple[str, ...]
    conductance: float = 1.0

    def calc_heat_flow(self, temperature: NDArray[np.float64]) -> NDArray[np.float64]:
        return np.zeros_like(temperature)


_NODES = tuple(Node(name, ()) for name in ("PRP", "RGP", "Can"))
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
    return Model(nodes, procedure or _two_phase_procedure(), np.full(len(nodes), 300.0), (2,), {})


def test_multi_phase_procedure_keeps_phases_in_order():
    phases = (
        Phase("precool", (_EXCHANGE_GAS,), _cutoff),
        Phase("transfer", (_EXCHANGE_GAS, _SWITCH), _cutoff),
        Phase("cooldown", (_SWITCH,), None),
    )

    procedure = Procedure(phases=phases, stop_condition=_cutoff)

    assert [phase.name for phase in procedure.phases] == ["precool", "transfer", "cooldown"]
    assert procedure.stop_condition is _cutoff
    assert procedure.phases[0].until(0.0, np.array([80.0, 90.0, 300.0])) == pytest.approx(10.0)


def test_link_may_change_parameters_between_phases():
    weak = _StubLink("switch", ("PRP", "RGP"), conductance=0.01)
    strong = _StubLink("switch", ("PRP", "RGP"), conductance=10.0)

    procedure = Procedure(
        phases=(Phase("open", (weak,), _cutoff), Phase("closed", (strong,), None)),
        stop_condition=None,
    )

    assert [phase.links[0].conductance for phase in procedure.phases] == [0.01, 10.0]


def test_single_phase_procedure_without_stop_condition_is_accepted():
    procedure = Procedure(phases=(Phase("only", (), None),), stop_condition=None)

    assert procedure.phases[0].until is None
    assert procedure.stop_condition is None


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


def test_model_holds_run_inputs():
    procedure = _two_phase_procedure()
    initial_temperature = np.array([300.0, 300.0, 295.0])

    model = Model(_NODES, procedure, initial_temperature, (2,), {"method": "Radau", "max_step": 100.0})

    assert model.nodes == _NODES
    assert model.procedure is procedure
    np.testing.assert_array_equal(model.initial_temperature, [300.0, 300.0, 295.0])
    assert model.fixed_boundary_nodes == (2,)
    assert model.solver_options["method"] == "Radau"


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


def test_result_holds_status_model_and_dataset():
    model = _model()
    data = xr.Dataset(
        {"temperature": (("time", "node"), np.full((2, 3), 300.0))},
        coords={"time": [0.0, 1.0], "node": ["PRP", "RGP", "Can"]},
    )

    result = Result("success", model, data)

    assert result.status == "success"
    assert result.model is model
    assert result.data is data


@pytest.mark.parametrize(
    ("component", "field"),
    [
        (Phase("only", (), None), "links"),
        (_two_phase_procedure(), "phases"),
        (_model(), "graphs"),
    ],
)
def test_run_types_are_immutable(component, field):
    with pytest.raises(dataclasses.FrozenInstanceError):
        setattr(component, field, None)
