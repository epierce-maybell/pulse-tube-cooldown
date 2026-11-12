"""Schematic drawing of a Graph."""

import re
from collections.abc import Callable, Mapping, Sequence

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.artist import Artist
from matplotlib.colors import to_rgba
from matplotlib.figure import Figure
from matplotlib.font_manager import FontProperties
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, Patch, Rectangle
from matplotlib.transforms import Bbox
from numpy.typing import NDArray

from pulse_tube_cooldown.core.types import Graph
from pulse_tube_cooldown.library.links.pulse_tube.head import PulseTubeHead

__all__ = ["plot_graph"]


def plot_graph(graph: Graph, layout: Mapping[str, tuple[int, int]] | None = None) -> Figure:
    """Draws a graph as a schematic of its nodes and links.

    Parallel links share a bundle, pulse-tube heads are rods through their stages, and links on more than two nodes
    are regions.

    Args:
        graph: The graph to draw.
        layout: Column and row of each node, from the top left. Nodes internal to a link are drawn within it. If
            omitted, nodes form one column in graph order.

    Returns:
        A figure sized to fit the drawing, with each link's artists given the link's name as their `gid`.
    """
    internal = {node.name for link in graph.links for node in getattr(link, "internal_nodes", ())}
    if layout is None:
        drawn = [node.name for node in graph.nodes if node.name not in internal]
        layout = {name: (0, row) for row, name in enumerate(drawn)}
    members = [tuple(name for name in link.nodes if name not in internal) for link in graph.links]
    rods, bundles, regions = classify(graph, members, layout)
    routes = {pair: route(pair, layout) for pair in bundles}

    class_names = [type(link).__name__ for link in graph.links]
    palette = ("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948")
    colors = {name: palette[i % len(palette)] for i, name in enumerate(dict.fromkeys(class_names))}
    rank = {name: i for i, name in enumerate(colors)}
    ink, muted, fill, edge = "#0b0b0b", "#52514e", "#f0efec", "#8a8984"

    def in_order(indices: Sequence[int]) -> list[int]:
        return sorted(indices, key=lambda i: (rank[class_names[i]], i))

    figure = plt.figure(figsize=(1.0, 1.0))
    ax = figure.add_axes((0.0, 0.0, 1.0, 1.0))
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.set_axis_off()
    renderer = figure.canvas.get_renderer()
    name_font = FontProperties(size="medium", weight="bold")
    label_font = FontProperties(size="small")
    title_font = FontProperties(size="large", weight="bold")

    def width(text: str, font: FontProperties) -> float:
        return renderer.get_text_width_height_descent(text, font, ismath=False)[0] / figure.dpi

    em = label_font.get_size_in_points() / 72.0
    line = 1.4 * em
    spacing = 0.5 * em
    swatch = 1.2 * em
    slab_height = 2.4 * em
    rod_width = 1.1 * em
    rod_pitch = 2.0 * em
    rod_gap = 0.35 * em

    def list_width(indices: Sequence[int]) -> float:
        return swatch + 0.6 * em + max(width(graph.links[i].name, label_font) for i in indices)

    rod_columns = [layout[members[i][0]][0] for i in rods]
    rod_slots = [rod_columns[:k].count(column) for k, column in enumerate(rod_columns)]
    rod_band = em + (max(rod_slots) + 1) * rod_pitch if rods else 0.0
    name_band = max(width(name, name_font) for name in layout) + 3.0 * em
    slab_width = rod_band + name_band

    def spine(name: str) -> float:
        return rod_band + name_band / 2.0 if layout[name][0] in rod_columns else slab_width / 2.0

    brackets = [pair for pair, kind in routes.items() if kind == "bracket"]
    lanes = dict(zip(brackets, bracket_lanes([tuple(layout[name] for name in pair) for pair in brackets]), strict=True))
    lane_pitch = 1.5 * em + max((len(bundles[pair]) * spacing for pair in brackets), default=0.0)
    column_lanes = {
        layout[pair[0]][0]: max(lane + 1 for other, lane in lanes.items() if layout[other[0]][0] == layout[pair[0]][0])
        for pair in brackets
    }

    columns = max(column for column, _ in layout.values()) + 1
    right, left = [0.0] * columns, [0.0] * columns
    gap = 3.0 * em
    for pair, indices in bundles.items():
        (first_column, _), (second_column, _) = layout[pair[0]], layout[pair[1]]
        match routes[pair]:
            case "down":
                extent = spine(pair[0]) + len(indices) * spacing / 2.0 + em + list_width(indices) - slab_width
                right[first_column] = max(right[first_column], extent)
                gap = max(gap, len(indices) * line + 2.0 * em)
            case "bracket":
                extent = (column_lanes[first_column] + 0.5) * lane_pitch + em + list_width(indices)
                left[first_column] = max(left[first_column], extent)
            case _:
                column = max(first_column, second_column)
                left[column] = max(left[column], em + list_width(indices))
                gap = max(gap, len(indices) * line + 3.0 * em)
    for i in rods:
        gap = max(gap, width(graph.links[i].name, label_font) + 2.0 * em)
    pitch = slab_height + gap
    column_x = [0.0]
    for column in range(1, columns):
        column_x.append(column_x[-1] + slab_width + max(right[column - 1] + left[column] + 2.0 * em, 4.0 * em))

    def x_of(name: str) -> float:
        return column_x[layout[name][0]]

    def y_of(name: str) -> float:
        return -layout[name][1] * pitch

    node_artists: dict[str, list[Artist]] = {}
    link_artists: dict[int, list[Artist]] = {i: [] for i in range(len(graph.links))}

    for name in layout:
        x, y = x_of(name), y_of(name)
        box = FancyBboxPatch(
            (x, y - slab_height / 2.0),
            slab_width,
            slab_height,
            boxstyle=f"round,pad=0,rounding_size={0.4 * em}",
            facecolor=fill,
            edgecolor=edge,
            linewidth=1.0,
            zorder=3,
        )
        ax.add_patch(box)
        text = ax.text(
            x + spine(name), y, name, ha="center", va="center", fontproperties=name_font, color=ink, zorder=5
        )
        node_artists[name] = [box, text]

    def strand(i: int, xs: Sequence[float], ys: Sequence[float]) -> None:
        (artist,) = ax.plot(
            xs,
            ys,
            color=colors[class_names[i]],
            linewidth=1.5,
            solid_capstyle="butt",
            solid_joinstyle="miter",
            gid=graph.links[i].name,
            zorder=1,
        )
        link_artists[i].append(artist)

    def listing(indices: Sequence[int], x: float, y: float, align: str, anchor: str) -> None:
        rows = in_order(indices)
        top = y + (len(rows) - 1) * line / (2.0 if anchor == "center" else 1.0)
        sign = 1.0 if align == "left" else -1.0
        for n, i in enumerate(rows):
            row_y = top - n * line
            (mark,) = ax.plot(
                [x, x + sign * swatch],
                [row_y, row_y],
                color=colors[class_names[i]],
                linewidth=2.2,
                solid_capstyle="butt",
                gid=graph.links[i].name,
            )
            text = ax.text(
                x + sign * (swatch + 0.6 * em),
                row_y,
                graph.links[i].name,
                ha=align,
                va="center",
                fontproperties=label_font,
                color=ink,
            )
            link_artists[i] += [mark, text]

    bracket_y: dict[tuple[str, str], float] = {}
    lowest: dict[int, float] = {}
    for pair in sorted(brackets, key=lambda pair: -(y_of(pair[0]) + y_of(pair[1]))):
        column, height = layout[pair[0]][0], len(bundles[pair]) * line
        center = min((y_of(pair[0]) + y_of(pair[1])) / 2.0, lowest.get(column, np.inf) - height / 2.0)
        bracket_y[pair] = center
        lowest[column] = center - height / 2.0

    for pair, indices in bundles.items():
        upper, lower = pair
        rows = in_order(indices)
        offsets = (np.arange(len(rows)) - (len(rows) - 1) / 2.0) * spacing
        half = len(rows) * spacing / 2.0
        top, bottom = y_of(upper), y_of(lower)
        match routes[pair]:
            case "down":
                x = x_of(upper) + spine(upper)
                for i, offset in zip(rows, offsets, strict=True):
                    strand(i, [x + offset] * 2, [top - slab_height / 2.0, bottom + slab_height / 2.0])
                listing(rows, x + half + em, (top + bottom) / 2.0, "left", "center")
            case "bracket":
                x = x_of(upper)
                lane_x = x - (lanes[pair] + 1) * lane_pitch
                for i, offset in zip(rows, offsets, strict=True):
                    xs = [x, lane_x + offset, lane_x + offset, x]
                    strand(i, xs, [top - offset, top - offset, bottom + offset, bottom + offset])
                label_x = x - (column_lanes[layout[upper][0]] + 0.5) * lane_pitch
                listing(rows, label_x, bracket_y[pair], "right", "center")
            case "across":
                west, east = sorted(pair, key=lambda name: layout[name][0])
                for i, offset in zip(rows, offsets, strict=True):
                    strand(i, [x_of(west) + slab_width, x_of(east)], [top - offset] * 2)
                listing(rows, x_of(east) - em, top + half + 0.6 * em, "right", "bottom")
            case "elbow":
                direction = np.sign(layout[lower][0] - layout[upper][0])
                start = x_of(upper) + (slab_width if direction > 0 else 0.0)
                corner = x_of(lower) + spine(lower)
                for i, offset in zip(rows, offsets, strict=True):
                    xs = [start, corner - direction * offset, corner - direction * offset]
                    strand(i, xs, [top - offset, top - offset, bottom + slab_height / 2.0])
                label_x = max(x_of(upper), x_of(lower))
                listing(rows, label_x - em, top + half + 0.6 * em, "right", "bottom")

    for i, slot in zip(rods, rod_slots, strict=True):
        head = graph.links[i]
        color = colors[class_names[i]]
        x = x_of(members[i][0]) + em + slot * rod_pitch + (rod_pitch - rod_width) / 2.0
        warm = y_of(head.stages[0]) + slab_height / 2.0 + gap - em
        for low, high in head_spans(head, warm, y_of, slab_height):
            segment = FancyBboxPatch(
                (x, low + rod_gap / 2.0),
                rod_width,
                high - low - rod_gap,
                boxstyle=f"round,pad=0,rounding_size={rod_width / 2.0}",
                facecolor=color,
                linewidth=0,
                zorder=2,
            )
            ax.add_patch(segment)
            link_artists[i].append(segment)
        for name in members[i]:
            collar = Rectangle(
                (x, y_of(name) - slab_height / 2.0), rod_width, slab_height, facecolor=color, linewidth=0, zorder=4
            )
            ax.add_patch(collar)
            link_artists[i].append(collar)
        bottom = min(y_of(name) for name in members[i]) - slab_height / 2.0 - em
        label = ax.text(
            x + rod_width / 2.0,
            bottom,
            head.name,
            rotation=90,
            rotation_mode="anchor",
            ha="right",
            va="center",
            fontproperties=label_font,
            color=ink,
        )
        link_artists[i].append(label)
        for artist in link_artists[i]:
            artist.set_gid(head.name)

    def bounds(artists: Sequence[Artist]) -> Bbox:
        to_data = ax.transData.inverted()
        return Bbox.union([artist.get_window_extent(renderer).transformed(to_data) for artist in artists])

    for i in sorted(regions, key=lambda i: len(members[i])):
        enclosed = set(members[i])
        contents = [artist for name in enclosed for artist in node_artists[name]]
        contents += [
            artist
            for j, artists in link_artists.items()
            if j != i and set(members[j]) <= enclosed
            for artist in artists
        ]
        box = bounds(contents).padded(1.5 * em)
        title = 3.2 * em
        color = colors[class_names[i]]
        for style in (
            {"facecolor": to_rgba(color, 0.07), "edgecolor": "none"},
            {"facecolor": "none", "edgecolor": color, "linewidth": 1.0, "linestyle": (0, (4, 3))},
        ):
            region = FancyBboxPatch(
                (box.x0, box.y0),
                box.width,
                box.height + title,
                boxstyle=f"round,pad=0,rounding_size={2.0 * em}",
                gid=graph.links[i].name,
                zorder=0,
                **style,
            )
            ax.add_patch(region)
            link_artists[i].append(region)
        top = box.y1 + title - em
        link_artists[i].append(
            ax.text(box.x0 + 1.5 * em, top, graph.links[i].name, va="top", fontproperties=title_font, color=ink)
        )
        if words(class_names[i]).casefold() != graph.links[i].name.casefold():
            link_artists[i].append(
                ax.text(
                    box.x0 + 1.5 * em,
                    top - 1.6 * em,
                    words(class_names[i]),
                    va="top",
                    fontproperties=label_font,
                    color=muted,
                )
            )

    handles: list[Artist] = []
    for name, color in colors.items():
        first = class_names.index(name)
        if first in regions:
            handles.append(
                Patch(facecolor=to_rgba(color, 0.07), edgecolor=color, linestyle=(0, (4, 3)), label=words(name))
            )
        else:
            thickness = 7.0 if first in rods else 2.2
            handles.append(Line2D([], [], color=color, linewidth=thickness, label=words(name)))

    content = bounds(ax.patches + ax.lines + ax.texts)
    legend = ax.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=((content.x0 + content.x1) / 2.0, content.y0 - em),
        bbox_transform=ax.transData,
        ncols=len(handles),
        frameon=False,
        prop=label_font,
        handlelength=2.5,
    )
    fitted = Bbox.union([content, bounds([legend])]).padded(em)
    ax.set_xlim(fitted.x0, fitted.x1)
    ax.set_ylim(fitted.y0, fitted.y1)
    figure.set_size_inches(fitted.width, fitted.height)
    return figure


def classify(
    graph: Graph, members: Sequence[tuple[str, ...]], layout: Mapping[str, tuple[int, int]]
) -> tuple[list[int], dict[tuple[str, str], list[int]], list[int]]:
    """Splits link indices into rods, bundles by node pair, and regions."""
    rods: list[int] = []
    bundles: dict[tuple[str, str], list[int]] = {}
    regions: list[int] = []
    for i, (link, names) in enumerate(zip(graph.links, members, strict=True)):
        if isinstance(link, PulseTubeHead) and len({layout[name][0] for name in names}) == 1:
            rods.append(i)
        elif len(names) == 2:
            pair = tuple(sorted(names, key=lambda name: layout[name][::-1]))
            bundles.setdefault(pair, []).append(i)
        else:
            regions.append(i)
    return rods, bundles, regions


def route(pair: tuple[str, str], layout: Mapping[str, tuple[int, int]]) -> str:
    """Returns how a bundle runs: across, down, bracket, or elbow."""
    (first_column, first_row), (second_column, second_row) = layout[pair[0]], layout[pair[1]]
    if first_row == second_row:
        return "across"
    if first_column != second_column:
        return "elbow"
    if any(column == first_column and first_row < row < second_row for column, row in layout.values()):
        return "bracket"
    return "down"


def bracket_lanes(cells: Sequence[tuple[tuple[int, int], tuple[int, int]]]) -> list[int]:
    """Assigns each bracket the innermost lane free of overlapping brackets in its column."""
    lanes = [0] * len(cells)
    placed: list[int] = []
    for i in sorted(range(len(cells)), key=lambda i: cells[i][1][1] - cells[i][0][1]):
        (column, top), (_, bottom) = cells[i]
        taken = {
            lanes[j]
            for j in placed
            if cells[j][0][0] == column and min(bottom, cells[j][1][1]) >= max(top, cells[j][0][1])
        }
        lanes[i] = min(set(range(len(taken) + 1)) - taken)
        placed.append(i)
    return lanes


def head_spans(
    head: PulseTubeHead, warm: float, y_of: Callable[[str], float], slab_height: float
) -> list[tuple[float, float]]:
    """Returns the bottom and top height of each regenerator node along a head's rod, warm to cold."""
    first, second = head.stages
    half = slab_height / 2.0
    pieces_2 = [(0.0, 1.0, y_of(first) - half, y_of(second) + half)]
    if head.clamp is not None and y_of(second) < y_of(head.clamp) < y_of(first):
        clamp = tiled_edges(head.regenerator_2.midpoint)[head.memory_node]
        pieces_2 = [
            (0.0, clamp, y_of(first) - half, y_of(head.clamp) + half),
            (clamp, 1.0, y_of(head.clamp) - half, y_of(second) + half),
        ]
    return [
        *node_spans(head.regenerator_1.midpoint, [(0.0, 1.0, warm, y_of(first) + half)]),
        *node_spans(head.regenerator_2.midpoint, pieces_2),
    ]


def node_spans(
    midpoint: NDArray[np.float64], pieces: Sequence[tuple[float, float, float, float]]
) -> list[tuple[float, float]]:
    """Returns the bottom and top height of each regenerator node, mapped through `pieces`."""
    edges = tiled_edges(midpoint)
    spans = []
    for low, high, middle in zip(edges[:-1], edges[1:], midpoint, strict=True):
        start, end, start_height, end_height = next(piece for piece in pieces if piece[0] <= middle <= piece[1])
        heights = [start_height + (edge - start) / (end - start) * (end_height - start_height) for edge in (low, high)]
        spans.append((min(heights), max(heights)))
    return spans


def tiled_edges(midpoint: NDArray[np.float64]) -> NDArray[np.float64]:
    """Returns the edges of nodes tiling [0, 1], from their midpoints."""
    edges = [0.0]
    for middle in midpoint:
        edges.append(2.0 * middle - edges[-1])
    return np.array(edges)


def words(class_name: str) -> str:
    """Converts a CamelCase class name to sentence case."""
    return re.sub(r"(?<=[a-z])(?=[A-Z])", " ", class_name).capitalize()
