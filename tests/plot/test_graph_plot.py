"""Tests for drawing a Graph."""

import dataclasses
import itertools

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest
from matplotlib.text import Text
from numpy.typing import NDArray

from pulse_tube_cooldown.config import load_link
from pulse_tube_cooldown.core.types import EnergyBalanceTerms, Graph, Material, Node
from pulse_tube_cooldown.plot.graph import plot_graph

_TOY = Material(
    name="Toy",
    density=1.0,
    specific_heat=lambda t: t,
    specific_heat_integral=lambda t: t**2 / 2.0,
    thermal_conductivity=lambda t: t,
    conductivity_integral=lambda t: t**2 / 2.0,
)


@dataclasses.dataclass(frozen=True, slots=True)
class Strap:
    name: str
    nodes: tuple[str, ...]

    def calc_energy_balance_terms(self, temperature: NDArray[np.float64]) -> EnergyBalanceTerms:
        return EnergyBalanceTerms(heat_flow=np.zeros_like(temperature))


@dataclasses.dataclass(frozen=True, slots=True)
class Bath(Strap):
    pass


_LAYOUT = {
    "Top": (0, 0),
    "PRP": (0, 1),
    "RGP": (0, 2),
    "CFP": (0, 3),
    "STP": (0, 4),
    "PRP shield": (1, 1),
    "STP shield": (1, 4),
}

_LINKS = (
    Strap("Top supports", ("Top", "PRP")),
    *(Strap(f"PRP-RGP strap with a long name {k}", ("PRP", "RGP")) for k in range(4)),
    Bath("PRP-RGP bath", ("RGP", "PRP")),
    *(Strap(f"RGP-CFP strap {k}", ("RGP", "CFP")) for k in range(3)),
    Strap("CFP-STP strap", ("CFP", "STP")),
    Strap("PRP-CFP bypass", ("PRP", "CFP")),
    Strap("Top-STP suspension", ("Top", "STP")),
    *(Strap(f"PRP shield strap {k}", ("PRP", "PRP shield")) for k in range(2)),
    Strap("STP shield strap", ("STP", "STP shield")),
    Bath("Shield radiation", ("STP shield", "PRP shield")),
    Bath("Top-shield blanket", ("Top", "PRP shield")),
    Bath("Gas", ("RGP", "CFP", "STP")),
)


def _graph(*links) -> Graph:
    names = [name for name in _LAYOUT if any(name in link.nodes for link in links)]
    nodes = [Node(name, ((_TOY, 1.0),)) for name in names]
    nodes += [node for link in links for node in getattr(link, "internal_nodes", ())]
    return Graph(nodes=tuple(nodes), links=links)


@pytest.fixture(scope="module")
def head():
    return load_link(
        {
            "name": "Head",
            "type": "pulse_tube_head",
            "from_archetype": True,
            "stages": ["PRP", "STP"],
            "clamp": "RGP",
            "contact_resistance": 1.0,
            "warm_flange_temperature": 300.0,
        }
    )


@pytest.fixture
def crowded(head):
    figure = plot_graph(_graph(*_LINKS, head), _LAYOUT)
    figure.canvas.draw()
    yield figure
    plt.close(figure)


def _texts(figure) -> list[Text]:
    return [text for text in figure.findobj(Text) if text.get_visible() and text.get_text()]


def _names(figure) -> set[str]:
    return {text.get_text() for text in _texts(figure)}


def test_names_every_node_and_link(crowded, head):
    names = _names(crowded)

    assert set(_LAYOUT) <= names
    assert {link.name for link in (*_LINKS, head)} <= names


def test_internal_nodes_are_not_drawn_as_named_nodes(crowded, head):
    assert not {node.name for node in head.internal_nodes} & _names(crowded)


def test_texts_do_not_overlap(crowded):
    renderer = crowded.canvas.get_renderer()
    boxes = [text.get_window_extent(renderer) for text in _texts(crowded)]

    for first, second in itertools.combinations(boxes, 2):
        assert not first.overlaps(second)


def test_link_texts_stay_clear_of_node_boxes(crowded):
    renderer = crowded.canvas.get_renderer()
    ax = crowded.axes[0]
    boxes = [patch.get_window_extent(renderer) for patch in ax.patches if patch.get_gid() is None]
    link_texts = [text for text in _texts(crowded) if text.get_text() not in _LAYOUT]

    for text in link_texts:
        extent = text.get_window_extent(renderer)
        assert not any(extent.overlaps(box) for box in boxes), text.get_text()


def test_figure_holds_every_text(crowded):
    renderer = crowded.canvas.get_renderer()
    frame = crowded.bbox

    for text in _texts(crowded):
        box = text.get_window_extent(renderer)
        assert frame.x0 <= box.x0 and box.x1 <= frame.x1 and frame.y0 <= box.y0 and box.y1 <= frame.y1


def test_regenerator_nodes_lie_between_the_boxes_of_their_ends(crowded, head):
    renderer = crowded.canvas.get_renderer()
    centers = {text.get_text(): text.get_window_extent(renderer) for text in _texts(crowded)}
    boxes = {}
    for patch in crowded.axes[0].patches:
        extent = patch.get_window_extent(renderer)
        for name in ("PRP", "RGP", "STP"):
            center = centers[name]
            if patch.get_gid() is None and extent.contains((center.x0 + center.x1) / 2, (center.y0 + center.y1) / 2):
                boxes[name] = extent
    segments = [
        extent
        for patch in crowded.axes[0].patches
        if patch.get_gid() == head.name
        and not any((extent := patch.get_window_extent(renderer)).overlaps(box) for box in boxes.values())
    ]

    above_stage_1 = [s for s in segments if s.y0 >= boxes["PRP"].y1]
    above_clamp = [s for s in segments if boxes["RGP"].y1 <= s.y0 and s.y1 <= boxes["PRP"].y0]
    below_clamp = [s for s in segments if boxes["STP"].y1 <= s.y0 and s.y1 <= boxes["RGP"].y0]
    assert len(segments) == len(head.internal_nodes)
    assert len(above_stage_1) == len(head.regenerator_1.nodes)
    assert len(above_clamp) == head.memory_node
    assert len(below_clamp) == len(head.regenerator_2.nodes) - head.memory_node


def test_without_layout_nodes_form_a_column_in_graph_order():
    graph = _graph(Strap("Rod", ("PRP", "RGP")), Strap("Bar", ("RGP", "CFP")))
    figure = plot_graph(graph)
    figure.canvas.draw()
    try:
        renderer = figure.canvas.get_renderer()
        positions = {text.get_text(): text.get_window_extent(renderer) for text in _texts(figure)}
        centers = np.array([[(positions[n].x0 + positions[n].x1) / 2, positions[n].y0] for n in ("PRP", "RGP", "CFP")])
        np.testing.assert_allclose(centers[:, 0], centers[0, 0])
        assert np.all(np.diff(centers[:, 1]) < 0.0)
    finally:
        plt.close(figure)
