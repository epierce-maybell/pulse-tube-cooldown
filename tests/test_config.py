"""Tests for loading configuration files."""

import tomllib

import numpy as np
import pytest

from pulse_tube_cooldown.config import configure_model, load_fridge, load_procedure
from pulse_tube_cooldown.core.errors import ConfigurationError
from pulse_tube_cooldown.core.types import Graph, Model
from pulse_tube_cooldown.library.links.conduction import Conduction

_FRIDGE = tomllib.loads(
    """
    [[nodes]]
    name = "Plate"
    composition = [["OFHC Cu RRR100", 1.5]]

    [[nodes]]
    name = "Mixed"
    composition = [["OFHC Cu RRR100", 2.0], ["304 SS", 0.5]]

    [[nodes]]
    name = "Bath"
    composition = [["304 SS", 3.0]]

    [[links]]
    name = "Rod A"
    type = "conduction"
    nodes = ["Plate", "Bath"]
    material = "304 SS"
    area = 1.0e-4
    length = 0.1

    [[links]]
    name = "Rod B"
    type = "conduction"
    nodes = ["Bath", "Mixed"]
    material = "OFHC Cu RRR100"
    area = 2.0e-5
    length = 0.3
    """
)


def _fridge(**changes):
    return _FRIDGE | changes


def test_loads_nodes_in_order_with_compositions():
    graph = load_fridge(_FRIDGE)

    assert isinstance(graph, Graph)
    assert [node.name for node in graph.nodes] == ["Plate", "Mixed", "Bath"]
    assert [[(material.name, mass) for material, mass in node.composition] for node in graph.nodes] == [
        [("OFHC Cu RRR100", 1.5)],
        [("OFHC Cu RRR100", 2.0), ("304 SS", 0.5)],
        [("304 SS", 3.0)],
    ]


def test_loads_links_with_parameters_and_connectivity():
    graph = load_fridge(_FRIDGE)

    assert all(isinstance(link, Conduction) for link in graph.links)
    assert [(link.name, link.nodes, link.material.name, link.area, link.length) for link in graph.links] == [
        ("Rod A", ("Plate", "Bath"), "304 SS", 1.0e-4, 0.1),
        ("Rod B", ("Bath", "Mixed"), "OFHC Cu RRR100", 2.0e-5, 0.3),
    ]
    assert graph.connectivity == ((0, 2), (2, 1))


def test_fridge_without_links_has_only_nodes():
    graph = load_fridge({"nodes": _FRIDGE["nodes"]})

    assert len(graph.nodes) == 3
    assert graph.links == ()


def test_link_to_unknown_node_names_link_and_node():
    links = [_FRIDGE["links"][0] | {"nodes": ["Plate", "Tub"]}]

    with pytest.raises(ConfigurationError, match=r"(?s)Rod A.*Tub"):
        load_fridge(_fridge(links=links))


def test_unknown_link_type_names_the_link():
    links = [_FRIDGE["links"][0] | {"type": "strut"}]

    with pytest.raises(ConfigurationError, match=r"(?s)Rod A.*strut"):
        load_fridge(_fridge(links=links))


def test_unknown_node_material_names_the_material():
    nodes = [{"name": "Plate", "composition": [["Unobtainium", 1.0]]}]

    with pytest.raises(ConfigurationError, match="Unobtainium"):
        load_fridge({"nodes": nodes})


def test_duplicate_node_names_are_rejected():
    nodes = _FRIDGE["nodes"] + [_FRIDGE["nodes"][0]]

    with pytest.raises(ConfigurationError, match="Plate"):
        load_fridge(_fridge(nodes=nodes))


@pytest.mark.parametrize("node", [{"composition": [["304 SS", 1.0]]}, {"name": "Plate"}])
def test_node_missing_a_key_is_rejected(node):
    with pytest.raises(ConfigurationError, match="name|composition"):
        load_fridge({"nodes": [node]})


_GRAPH = load_fridge(
    _FRIDGE
    | {
        "links": _FRIDGE["links"]
        + [
            {
                "name": "Rod C",
                "type": "conduction",
                "nodes": ["Plate", "Mixed"],
                "material": "304 SS",
                "area": 1e-4,
                "length": 0.1,
            }
        ]
    }
)

_PROCEDURE = tomllib.loads(
    """
    until = { condition = "gas_pushing_begins", node = "Mixed", threshold = 10.0 }

    [[phases]]
    name = "precool"
    disable = ["Rod B"]
    until = { condition = "gas_pushing_begins", node = "Plate", threshold = 100.0 }

    [[phases]]
    name = "swap"
    enable = ["Rod B"]
    disable = ["Rod C"]
    """
)


def _phases(*phases):
    return {"phases": list(phases)}


def _temperatures(**by_node):
    return np.array([by_node.get(node.name, 300.0) for node in _GRAPH.nodes])


def test_phases_activate_links_cumulatively_in_fridge_order():
    procedure = load_procedure(_PROCEDURE, _GRAPH)
    rod_a, rod_b, rod_c = _GRAPH.links

    assert [phase.name for phase in procedure.phases] == ["precool", "swap"]
    assert procedure.phases[0].links == (rod_a, rod_c)
    assert procedure.phases[1].links == (rod_a, rod_b)


def test_phase_without_changes_keeps_the_previous_links():
    until = _PROCEDURE["phases"][0]["until"]
    procedure = load_procedure(
        _phases(
            {"name": "precool", "disable": ["Rod B"], "until": until},
            {"name": "hold", "until": until},
            {"name": "final", "enable": ["Rod B"]},
        ),
        _GRAPH,
    )

    assert [link.name for link in procedure.phases[1].links] == ["Rod A", "Rod C"]
    assert [link.name for link in procedure.phases[2].links] == ["Rod A", "Rod B", "Rod C"]


@pytest.mark.parametrize(("plate", "sign"), [(120.0, 1.0), (100.0, 0.0), (80.0, -1.0)])
def test_phase_until_is_bound_to_its_node_and_threshold(plate, sign):
    procedure = load_procedure(_PROCEDURE, _GRAPH)

    assert np.sign(procedure.phases[0].until(0.0, _temperatures(Plate=plate, Mixed=1.0))) == sign
    assert procedure.phases[1].until is None


@pytest.mark.parametrize(("mixed", "sign"), [(20.0, 1.0), (10.0, 0.0), (5.0, -1.0)])
def test_procedure_until_is_bound_as_the_stop_condition(mixed, sign):
    procedure = load_procedure(_PROCEDURE, _GRAPH)

    assert np.sign(procedure.stop_condition(0.0, _temperatures(Mixed=mixed, Plate=1.0))) == sign


def test_procedure_without_until_has_no_stop_condition():
    procedure = load_procedure(_phases({"name": "only"}), _GRAPH)

    assert procedure.stop_condition is None
    assert procedure.phases[0].links == _GRAPH.links


@pytest.mark.parametrize(
    "phases",
    [
        pytest.param([{"name": "swap", "enable": ["Rod B"]}], id="enable-in-first-phase"),
        pytest.param(
            [
                {"name": "precool", "disable": ["Rod B"], "until": _PROCEDURE["phases"][0]["until"]},
                {"name": "swap", "disable": ["Rod B"]},
            ],
            id="disable-inactive",
        ),
        pytest.param(
            [
                {"name": "precool", "disable": ["Rod B"], "until": _PROCEDURE["phases"][0]["until"]},
                {"name": "swap", "enable": ["Rod B"], "disable": ["Rod B"]},
            ],
            id="enable-and-disable",
        ),
    ],
)
def test_redundant_activation_names_the_phase_and_link(phases):
    with pytest.raises(ConfigurationError, match=r"(?s)swap.*Rod B"):
        load_procedure(_phases(*phases), _GRAPH)


def test_unknown_link_names_the_phase_and_link():
    with pytest.raises(ConfigurationError, match=r"(?s)swap.*Rod D"):
        load_procedure(_phases({"name": "swap", "disable": ["Rod D"]}), _GRAPH)


@pytest.mark.parametrize(
    ("until", "named"),
    [
        ({"condition": "plate_melts", "node": "Plate", "threshold": 1.0}, "plate_melts"),
        ({"condition": "gas_pushing_begins", "node": "Tub", "threshold": 1.0}, "Tub"),
        ({"node": "Plate", "threshold": 1.0}, "condition"),
    ],
)
def test_bad_condition_names_the_problem(until, named):
    phases = _phases({"name": "precool", "until": until}, {"name": "final"})

    with pytest.raises(ConfigurationError, match=named):
        load_procedure(phases, _GRAPH)


_CONFIG = """
[meta]
experiment = "Toy conduction with a phase switch"
version = "abc1234"

[[fridge.nodes]]
name = "Plate"
composition = [["OFHC Cu RRR100", 1.0]]

[[fridge.nodes]]
name = "Bath"
composition = [["OFHC Cu RRR100", 1.0]]

[[fridge.links]]
name = "Rod A"
type = "conduction"
nodes = ["Plate", "Bath"]
material = "304 SS"
area = 1.0e-4
length = 0.1

[[fridge.links]]
name = "Rod B"
type = "conduction"
nodes = ["Plate", "Bath"]
material = "304 SS"
area = 1.0e-4
length = 0.1

[procedure]
until = { condition = "gas_pushing_begins", node = "Plate", threshold = 10.0 }

[[procedure.phases]]
name = "single rod"
disable = ["Rod B"]
until = { condition = "gas_pushing_begins", node = "Plate", threshold = 100.0 }

[[procedure.phases]]
name = "both rods"
enable = ["Rod B"]

[simulation]
solver = "Radau"
max_time = 86400.0
save_every = 100.0
rtol = 1.0e-6
initial_temperature = [300.0, 4.0]
fixed_boundary_nodes = ["Bath"]
"""


def _configure(tmp_path, text=_CONFIG):
    path = tmp_path / "toy.toml"
    path.write_text(text)
    return configure_model(path)


def test_configures_the_complete_model(tmp_path):
    model = _configure(tmp_path)

    assert isinstance(model, Model)
    assert model.meta == {"experiment": "Toy conduction with a phase switch", "version": "abc1234"}
    assert [node.name for node in model.nodes] == ["Plate", "Bath"]
    assert [phase.name for phase in model.procedure.phases] == ["single rod", "both rods"]
    assert [[link.name for link in graph.links] for graph in model.graphs] == [["Rod A"], ["Rod A", "Rod B"]]
    np.testing.assert_array_equal(model.initial_temperature, [300.0, 4.0])
    assert model.fixed_boundary_nodes == (1,)
    assert model.solver_options == {"solver": "Radau", "max_time": 86400.0, "save_every": 100.0, "rtol": 1.0e-6}


def test_procedure_conditions_are_bound_to_the_configured_nodes(tmp_path):
    model = _configure(tmp_path)

    assert model.procedure.phases[0].until(0.0, np.array([100.0, 4.0])) == 0.0
    assert model.procedure.stop_condition(0.0, np.array([10.0, 4.0])) == 0.0


def test_scalar_initial_temperature_applies_to_every_node(tmp_path):
    model = _configure(tmp_path, _CONFIG.replace("initial_temperature = [300.0, 4.0]", "initial_temperature = 295.0"))

    np.testing.assert_array_equal(model.initial_temperature, [295.0, 295.0])


def test_meta_and_fixed_boundary_nodes_are_optional(tmp_path):
    text = _CONFIG.replace('fixed_boundary_nodes = ["Bath"]', "")
    text = text.replace('[meta]\nexperiment = "Toy conduction with a phase switch"\nversion = "abc1234"', "")

    model = _configure(tmp_path, text)

    assert model.meta == {}
    assert model.fixed_boundary_nodes == ()


def test_initial_temperature_of_the_wrong_length_names_the_nodes(tmp_path):
    text = _CONFIG.replace("initial_temperature = [300.0, 4.0]", "initial_temperature = [300.0, 4.0, 4.0]")

    with pytest.raises(ConfigurationError, match=r"(?s)initial_temperature.*3.*Plate.*Bath"):
        _configure(tmp_path, text)


def test_unknown_fixed_boundary_node_is_named(tmp_path):
    text = _CONFIG.replace('fixed_boundary_nodes = ["Bath"]', 'fixed_boundary_nodes = ["Tub"]')

    with pytest.raises(ConfigurationError, match="Tub"):
        _configure(tmp_path, text)
