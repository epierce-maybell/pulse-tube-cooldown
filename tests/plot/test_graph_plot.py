"""Tests for drawing a Model."""

import dataclasses
import itertools

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest
from matplotlib.colors import to_hex
from matplotlib.patches import Patch
from matplotlib.text import Text
from numpy.typing import NDArray

from pulse_tube_cooldown.config import load_link
from pulse_tube_cooldown.core.types import EnergyBalanceTerms, Material, Model, Node, Phase, Procedure
from pulse_tube_cooldown.plot.graph import GraphColors, plot_graph

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
class Conduction(Strap):
    pass


@dataclasses.dataclass(frozen=True, slots=True)
class Radiation(Strap):
    pass


@dataclasses.dataclass(frozen=True, slots=True)
class GasGapSwitch(Strap):
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
    Radiation("PRP-RGP bath", ("RGP", "PRP")),
    *(Conduction(f"RGP-CFP strap {k}", ("RGP", "CFP")) for k in range(3)),
    GasGapSwitch("CFP-STP switch", ("CFP", "STP")),
    Strap("PRP-CFP bypass", ("PRP", "CFP")),
    Strap("Top-STP suspension", ("Top", "STP")),
    *(Strap(f"PRP shield strap {k}", ("PRP", "PRP shield")) for k in range(2)),
    Strap("STP shield strap", ("STP", "STP shield")),
    Radiation("Shield radiation", ("STP shield", "PRP shield")),
    Radiation("Top-shield blanket", ("Top", "PRP shield")),
    Radiation("Gas", ("RGP", "CFP", "STP")),
    Radiation("Shield gas", ("PRP", "PRP shield", "STP shield")),
)


def _model(*phase_links, fixed: tuple[str, ...] = (), layout=_LAYOUT) -> Model:
    """A model with one phase per tuple of links; nodes in `layout` order plus each head's internal nodes."""
    links = [link for links in phase_links for link in links]
    nodes = [Node(name, ((_TOY, 1.0),)) for name in layout if any(name in link.nodes for link in links)]
    nodes += [node for link in links for node in getattr(link, "internal_nodes", ())]
    count = len(phase_links)
    phases = tuple(
        Phase(f"phase {i}", tuple(links), None if i == count - 1 else (lambda t, temperature: 1.0))
        for i, links in enumerate(phase_links)
    )
    names = [node.name for node in nodes]
    return Model(
        meta={},
        nodes=tuple(nodes),
        procedure=Procedure(phases=phases, stop_condition=None),
        initial_temperature=np.full(len(nodes), 300.0),
        fixed_boundary_nodes=tuple(names.index(name) for name in fixed),
        solver_options={},
    )


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
    figure = plot_graph(_model((*_LINKS, head)), _LAYOUT)
    figure.canvas.draw()
    yield figure
    plt.close(figure)


def _texts(figure) -> list[Text]:
    return [text for text in figure.findobj(Text) if text.get_visible() and text.get_text()]


def _names(figure) -> set[str]:
    return {text.get_text() for text in _texts(figure)}


def _patches(figure, gid: str) -> list[Patch]:
    return [patch for patch in figure.axes[0].patches if patch.get_gid() == gid]


def test_names_every_node_and_no_regenerator_node(crowded, head):
    """Catches a node left unlabeled, or a link's internal node drawn as if it were a node."""
    names = _names(crowded)

    assert set(_LAYOUT) <= names
    assert not {node.name for node in head.internal_nodes} & names


def test_texts_do_not_overlap(crowded):
    """Catches labels colliding in a crowded drawing."""
    renderer = crowded.canvas.get_renderer()
    boxes = [text.get_window_extent(renderer) for text in _texts(crowded)]

    for first, second in itertools.combinations(boxes, 2):
        assert not first.overlaps(second)


def test_figure_holds_every_text(crowded):
    """Catches a figure sized too small, clipping labels or the legend."""
    renderer = crowded.canvas.get_renderer()
    frame = crowded.bbox

    for text in _texts(crowded):
        box = text.get_window_extent(renderer)
        assert frame.x0 <= box.x0 and box.x1 <= frame.x1 and frame.y0 <= box.y0 and box.y1 <= frame.y1


def test_head_draws_one_patch_per_regenerator_node(crowded, head):
    """Catches missing or extra regenerator cells in a pulse-tube head."""
    assert len(_patches(crowded, head.name)) == len(head.internal_nodes)


def test_pipe_has_a_joint_at_each_node(crowded):
    """Catches a same-column multi-node link drawn without one joint per node."""
    assert len(_patches(crowded, "Gas")) == 3


def test_region_encloses_the_boxes_of_its_nodes(crowded):
    """Catches a cross-column multi-node link whose outline misses one of its nodes."""
    renderer = crowded.canvas.get_renderer()
    outlines = _patches(crowded, "Shield gas")
    boxes = [_patches(crowded, name)[0].get_window_extent(renderer) for name in ("PRP", "PRP shield", "STP shield")]

    assert any(
        all(o.x0 <= b.x0 and b.x1 <= o.x1 and o.y0 <= b.y0 and b.y1 <= o.y1 for b in boxes)
        for o in (outline.get_window_extent(renderer) for outline in outlines)
    )


def test_links_of_every_phase_are_drawn():
    """Catches drawing only the first phase's links."""
    model = _model((Strap("Rod", ("PRP", "RGP")),), (Strap("Rod", ("PRP", "RGP")), Radiation("Late gas", ("PRP", "RGP", "CFP"))))
    figure = plot_graph(model, _LAYOUT)
    try:
        assert figure.findobj(lambda artist: artist.get_gid() == "Late gas")
    finally:
        plt.close(figure)


def test_sensor_nodes_have_heavier_outline_and_fixed_nodes_are_hatched():
    """Catches sensors and fixed-temperature nodes drawn like plain nodes."""
    model = _model((Strap("Rod", ("PRP", "RGP")), Strap("Bar", ("RGP", "CFP"))), fixed=("CFP",))
    figure = plot_graph(model, _LAYOUT, sensors=("PRP",))
    try:
        box = {name: _patches(figure, name)[0] for name in ("PRP", "RGP", "CFP")}
        assert box["PRP"].get_linewidth() > box["RGP"].get_linewidth()
        assert box["CFP"].get_hatch()
        assert not box["RGP"].get_hatch() and not box["PRP"].get_hatch()
    finally:
        plt.close(figure)


def test_colors_set_text_and_node_boxes():
    """Catches `colors` ignored for text, node fills, or node outlines."""
    model = _model((Strap("Rod", ("PRP", "RGP")),))
    colors = GraphColors(ink="#ffffff", fill="#000000", edge="#808080", hatch="#404040", rod="#c0c0c0")
    figure = plot_graph(model, _LAYOUT, colors=colors)
    try:
        box = _patches(figure, "RGP")[0]
        assert to_hex(box.get_facecolor()) == colors.fill
        assert to_hex(box.get_edgecolor()) == colors.edge
        assert {to_hex(text.get_color()) for text in _texts(figure)} == {colors.ink}
    finally:
        plt.close(figure)


def test_without_layout_nodes_form_a_column_in_model_order():
    """Catches the default layout scrambling or spreading nodes out of one column."""
    model = _model((Strap("Rod", ("PRP", "RGP")), Strap("Bar", ("RGP", "CFP"))))
    figure = plot_graph(model)
    figure.canvas.draw()
    try:
        renderer = figure.canvas.get_renderer()
        extents = [_patches(figure, name)[0].get_window_extent(renderer) for name in ("PRP", "RGP", "CFP")]
        centers = np.array([[(e.x0 + e.x1) / 2, (e.y0 + e.y1) / 2] for e in extents])
        np.testing.assert_allclose(centers[:, 0], centers[0, 0])
        assert np.all(np.diff(centers[:, 1]) < 0.0)
    finally:
        plt.close(figure)


def test_legend_names_link_types_in_sentence_case_and_heads(crowded, head):
    """Catches a legend missing a link type (registered or not), showing class names raw, or omitting the head."""
    names = _names(crowded)

    assert {"Strap", "Conduction", "Radiation", "Gas gap switch"} <= names
    assert any(head.name in name for name in names - set(_LAYOUT))
