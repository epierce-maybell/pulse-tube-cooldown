"""Tests for drawing a Graph."""

import dataclasses

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest
from numpy.typing import NDArray

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


def _graph(*links) -> Graph:
    nodes = tuple(Node(name, ((_TOY, 1.0),)) for name in ("PRP", "RGP", "CFP", "STP"))
    return Graph(nodes=nodes, links=links)


@pytest.fixture
def ax():
    fig, ax = plt.subplots()
    yield ax
    plt.close(fig)


def _text_position(ax, text):
    (match,) = [artist for artist in ax.texts if artist.get_text() == text]
    return np.array(match.get_position())


def _edges(ax, link_name):
    return [line for line in ax.get_lines() if line.get_gid() == link_name]


def _ends(line):
    xs, ys = line.get_xdata(), line.get_ydata()
    return np.array([xs[0], ys[0]]), np.array([xs[-1], ys[-1]])


def test_nodes_form_a_column_in_graph_order_first_on_top(ax):
    plot_graph(_graph(), ax=ax)

    positions = np.array([_text_position(ax, name) for name in ("PRP", "RGP", "CFP", "STP")])
    np.testing.assert_allclose(positions[:, 0], positions[0, 0])
    assert np.all(np.diff(positions[:, 1]) < 0.0)


def test_labels_each_link_with_its_name_and_class_name(ax):
    plot_graph(_graph(Strap("Rod", ("PRP", "RGP")), Bath("Gas", ("RGP", "CFP", "STP"))), ax=ax)

    texts = [text.get_text() for text in ax.texts]
    assert any("Rod" in text and "Strap" in text for text in texts)
    assert any("Gas" in text and "Bath" in text for text in texts)


def test_two_node_link_is_one_edge_between_its_nodes(ax):
    plot_graph(_graph(Strap("Rod", ("PRP", "CFP"))), ax=ax)

    (edge,) = _edges(ax, "Rod")
    start, end = _ends(edge)
    np.testing.assert_allclose(start, _text_position(ax, "PRP"))
    np.testing.assert_allclose(end, _text_position(ax, "CFP"))


def test_multi_node_link_joins_one_point_to_each_of_its_nodes(ax):
    plot_graph(_graph(Bath("Gas", ("RGP", "CFP", "STP"))), ax=ax)

    edges = _edges(ax, "Gas")
    assert len(edges) == 3
    ends = [_ends(edge) for edge in edges]
    hubs = np.array([hub for hub, _ in ends])
    np.testing.assert_allclose(hubs, [hubs[0]] * 3)
    for (_, node_end), name in zip(ends, ("RGP", "CFP", "STP"), strict=True):
        np.testing.assert_allclose(node_end, _text_position(ax, name))


def test_single_node_link_is_an_edge_from_its_node_to_its_label(ax):
    plot_graph(_graph(Strap("Heater", ("CFP",))), ax=ax)

    (edge,) = _edges(ax, "Heater")
    start, end = _ends(edge)
    np.testing.assert_allclose(start, _text_position(ax, "CFP"))
    np.testing.assert_allclose(end, _text_position(ax, "Heater\nStrap"))


def test_parallel_links_are_drawn_apart(ax):
    plot_graph(_graph(Strap("Rod", ("PRP", "RGP")), Bath("Glow", ("PRP", "RGP"))), ax=ax)

    rod, glow = _text_position(ax, "Rod\nStrap"), _text_position(ax, "Glow\nBath")
    assert not np.allclose(rod, glow)
    (rod_edge,), (glow_edge,) = _edges(ax, "Rod"), _edges(ax, "Glow")
    assert not np.allclose(rod_edge.get_xydata(), glow_edge.get_xydata())


@pytest.mark.parametrize("link", [Strap("Rod", ("PRP", "RGP")), Strap("Heater", ("CFP",))])
def test_lane_gap_is_in_points(ax, link):
    distances = []
    for lane_gap in (0.0, 36.0):
        ax.clear()
        plot_graph(_graph(link), ax=ax, lane_gap=lane_gap)
        label = _text_position(ax, f"{link.name}\nStrap")
        node = _text_position(ax, link.nodes[0])
        distances.append(abs(label[0] - node[0]))

    half_inch = 0.5 * ax.figure.dpi / ax.get_window_extent().width
    assert distances[1] - distances[0] == pytest.approx(half_inch)


class MultiLayerInsulationBlanket(Strap):
    pass


_CROWDED = (
    Strap("PRP-RGP supports", ("PRP", "RGP")),
    MultiLayerInsulationBlanket("PRP blanket with a long name", ("PRP", "RGP")),
    Bath("PRP shield", ("RGP", "PRP")),
    Strap("RGP-CFP supports", ("RGP", "CFP")),
    Bath("RGP-STP switch", ("RGP", "STP")),
    Bath("Exchange gas", ("PRP", "RGP", "CFP", "STP")),
    Strap("Heater", ("STP",)),
    Bath("Second heater", ("STP",)),
)


@pytest.fixture
def crowded():
    ax = plot_graph(_graph(*_CROWDED))
    ax.figure.canvas.draw()
    yield ax
    plt.close(ax.figure)


def test_labels_do_not_overlap_however_long(crowded):
    boxes = [text.get_bbox_patch().get_window_extent() for text in crowded.texts]
    for i, first in enumerate(boxes):
        for second in boxes[i + 1 :]:
            assert not first.overlaps(second)


def test_no_edge_passes_behind_another_links_label(crowded):
    for link in _CROWDED:
        (label,) = [text for text in crowded.texts if text.get_text().startswith(f"{link.name}\n")]
        box = label.get_bbox_patch().get_window_extent()
        for line in crowded.get_lines():
            if line.get_gid() not in {None, link.name}:
                path = line.get_transform().transform_path(line.get_path())
                assert not path.intersects_bbox(box, filled=False), (line.get_gid(), link.name)


def test_returns_the_axes_it_drew_on(ax):
    assert plot_graph(_graph(), ax=ax) is ax


def test_draws_on_a_new_figure_without_axes():
    ax = plot_graph(_graph(Strap("Rod", ("PRP", "RGP"))))
    try:
        assert len(_edges(ax, "Rod")) == 1
    finally:
        plt.close(ax.figure)
