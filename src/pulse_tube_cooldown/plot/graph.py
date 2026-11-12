"""Schematic drawing of a Model's thermal network."""

import inspect
import re
from collections.abc import Callable, Collection, Mapping, Sequence

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.figure import Figure
from matplotlib.font_manager import FontProperties
from matplotlib.patches import Circle, FancyBboxPatch, Rectangle
from matplotlib.transforms import Bbox
from numpy.typing import NDArray

from pulse_tube_cooldown.config import LINK_TYPES
from pulse_tube_cooldown.core.types import Link, Model
from pulse_tube_cooldown.library.links.pulse_tube.head import PulseTubeHead

__all__ = ["plot_graph"]


def plot_graph(
    model: Model, layout: Mapping[str, tuple[int, int]] | None = None, sensors: Collection[str] = ()
) -> Figure:
    """Draws a model's nodes and links as a schematic.

    Each pair of nodes joined by two-node links has one line per link type, colored by type in the order types are
    registered in `LINK_TYPES`, followed by unregistered types in the order they first appear. A link on more than two
    nodes is a pipe with a joint at each of its nodes if they share a column, and a dashed region around them
    otherwise. Pulse-tube heads are rods through their stages with one cell per regenerator node, named in the legend
    from left to right. Fixed-temperature nodes are hatched and nodes with thermometers have a heavy outline.

    Args:
        model: The model to draw, with the links of every phase.
        layout: Column and row of each node, from the top left. Nodes internal to a link are drawn within it. If
            omitted, nodes form one column in model order.
        sensors: Names of nodes with thermometers.

    Returns:
        A figure sized to fit the drawing. Each node's box has the node's name as its `gid`. Each head's regenerator
        node cells, each pipe's line and joints, and each region's outline have the link's name as their `gid`.
    """
    links: dict[str, Link] = {}
    for phase in model.procedure.phases:
        for link in phase.links:
            links.setdefault(link.name, link)
    internal = {node.name for link in links.values() for node in getattr(link, "internal_nodes", ())}
    if layout is None:
        drawn = [node.name for node in model.nodes if node.name not in internal]
        layout = {name: (0, row) for row, name in enumerate(drawn)}
    boundary = {model.nodes[i].name for i in model.fixed_boundary_nodes}
    heads = [link for link in links.values() if isinstance(link, PulseTubeHead)]
    shared = [link for link in links.values() if len(link.nodes) > 2 and not isinstance(link, PulseTubeHead)]
    pipes = [link for link in shared if len({layout[name][0] for name in link.nodes}) == 1]
    regions = [link for link in shared if link not in pipes]
    order = link_type_order()
    order += [kind for kind in dict.fromkeys(type(link).__name__ for link in links.values()) if kind not in order]

    pairs: dict[tuple[str, str], list[str]] = {}
    for link in links.values():
        if len(link.nodes) == 2:
            kinds = pairs.setdefault(tuple(sorted(link.nodes, key=lambda name: layout[name][::-1])), [])
            if type(link).__name__ not in kinds:
                kinds.append(type(link).__name__)
    for kinds in pairs.values():
        kinds.sort(key=order.index)
    routes = {pair: route(pair, layout) for pair in pairs}

    palette = ("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948")
    colors = {kind: palette[i % len(palette)] for i, kind in enumerate(order)}
    ink, fill, edge, hatch, rod = "#1f1e1c", "#f4f3f0", "#a3a29c", "#cfcec7", "#46453f"

    figure = plt.figure(figsize=(1.0, 1.0))
    ax = figure.add_axes((0.0, 0.0, 1.0, 1.0))
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.set_axis_off()
    renderer = figure.canvas.get_renderer()
    name_font = FontProperties(size="medium")
    label_font = FontProperties(size="small")
    em = label_font.get_size_in_points() / 72.0
    point = 1.0 / 72.0

    def width(text: str, font: FontProperties) -> float:
        return renderer.get_text_width_height_descent(text, font, ismath=False)[0] / figure.dpi

    line_width, line_gap = 1.4 * point, 2.2 * point

    def thickness(pair: tuple[str, str]) -> float:
        count = len(pairs[pair])
        return count * line_width + (count - 1) * line_gap

    def offsets(pair: tuple[str, str]) -> NDArray[np.float64]:
        return (np.arange(len(pairs[pair])) - (len(pairs[pair]) - 1) / 2.0) * (line_width + line_gap)

    rod_width, rod_pitch = 1.1 * em, 1.8 * em
    rod_columns = {layout[head.stages[0]][0] for head in heads}
    rod_band = 0.8 * em + len(heads) * rod_pitch if heads else 0.0
    pipe_width, pipe_pitch, joint_radius = 2.0 * point, 1.4 * em, 0.38 * em
    columns = max(column for column, _ in layout.values()) + 1
    pipe_column = {link.name: layout[link.nodes[0]][0] for link in pipes}
    pipe_slot = {
        link.name: [pipe_column[other.name] for other in pipes[:k]].count(pipe_column[link.name])
        for k, link in enumerate(pipes)
    }
    pipe_band = [list(pipe_column.values()).count(column) * pipe_pitch for column in range(columns)]
    bar_width = [0.0] * columns
    for name, (column, _) in layout.items():
        band = rod_band if column in rod_columns else 0.0
        bar_width[column] = max(bar_width[column], band + width(name, name_font) + pipe_band[column] + 3.0 * em)

    def spine_offset(column: int) -> float:
        band = rod_band if column in rod_columns else 0.0
        return band + (bar_width[column] - band) / 2.0

    brackets = [pair for pair, kind in routes.items() if kind == "bracket"]
    spans = [(layout[upper][0], layout[upper][1], layout[lower][1]) for upper, lower in brackets]
    lanes = dict(zip(brackets, interval_lanes(spans), strict=True))
    lane_pitch = max((thickness(pair) for pair in brackets), default=0.0) + 1.2 * em
    elbows = {pair: elbow_ends(pair, layout) for pair, kind in routes.items() if kind == "elbow"}

    terminals: dict[tuple[str, str], list[tuple[float, tuple[str, str]]]] = {}
    for pair, kind in routes.items():
        upper, lower = pair
        match kind:
            case "bracket":
                terminals.setdefault((upper, "left"), []).append((lanes[pair] - len(brackets), pair))
                terminals.setdefault((lower, "left"), []).append((len(brackets) - lanes[pair], pair))
            case "across":
                west, east = sorted(pair, key=lambda name: layout[name][0])
                terminals.setdefault((west, "right"), []).append((0.0, pair))
                terminals.setdefault((east, "left"), []).append((0.0, pair))
            case "elbow":
                entered, start = elbows[pair]
                side = "right" if layout[start][0] < layout[entered][0] else "left"
                rank = 1.0 if layout[entered][1] < layout[start][1] else -1.0
                terminals.setdefault((start, side), []).append((rank, pair))

    stack_gap = 0.6 * em
    terminal_y: dict[tuple[str, tuple[str, str]], float] = {}
    edge_height: dict[str, float] = {}
    for (name, _), entries in terminals.items():
        entries.sort(key=lambda entry: -entry[0])
        heights = [thickness(pair) for _, pair in entries]
        total = sum(heights) + stack_gap * (len(heights) - 1)
        edge_height[name] = max(edge_height.get(name, 0.0), total)
        y = total / 2.0
        for height, (_, pair) in zip(heights, entries, strict=True):
            terminal_y[(name, pair)] = y - height / 2.0
            y -= height + stack_gap

    minimum = 2.2 * em
    bar_height = {name: max(minimum, edge_height.get(name, 0.0) + em) for name in layout}
    rows = max(row for _, row in layout.values()) + 1
    row_height = [
        max((bar_height[name] for name, (_, r) in layout.items() if r == row), default=minimum) for row in range(rows)
    ]
    gap = 2.4 * em
    row_y = [0.0]
    for row in range(1, rows):
        row_y.append(row_y[-1] - row_height[row - 1] / 2.0 - gap - row_height[row] / 2.0)

    left_room = [0.0] * columns
    for pair in brackets:
        column = layout[pair[0]][0]
        left_room[column] = max(left_room[column], (lanes[pair] + 1) * lane_pitch)
    column_x = [left_room[0]]
    for column in range(1, columns):
        column_x.append(column_x[-1] + bar_width[column - 1] + max(left_room[column], 7.0 * em))

    def x_of(name: str) -> float:
        return column_x[layout[name][0]]

    def y_of(name: str) -> float:
        return row_y[layout[name][1]]

    def right(name: str) -> float:
        return x_of(name) + bar_width[layout[name][0]]

    def top(name: str) -> float:
        return y_of(name) + bar_height[name] / 2.0

    def bottom(name: str) -> float:
        return y_of(name) - bar_height[name] / 2.0

    node_boxes: dict[str, FancyBboxPatch] = {}
    for name, (column, _) in layout.items():
        x, y, height = x_of(name), y_of(name), bar_height[name]
        hatching = {"hatch": "////", "hatchcolor": hatch} if name in boundary else {}
        outline = {"edgecolor": ink, "linewidth": 1.6} if name in sensors else {"edgecolor": edge, "linewidth": 0.8}
        node_boxes[name] = FancyBboxPatch(
            (x, y - height / 2.0),
            bar_width[column],
            height,
            boxstyle=f"round,pad=0,rounding_size={0.3 * em}",
            facecolor=fill,
            gid=name,
            zorder=3,
            **outline,
            **hatching,
        )
        ax.add_patch(node_boxes[name])
        band = rod_band if column in rod_columns else 0.0
        ax.text(x + band + em, y, name, ha="left", va="center", fontproperties=name_font, color=ink, zorder=5)

    for link in pipes:
        color = colors[type(link).__name__]
        x = right(link.nodes[0]) - 1.2 * em - pipe_slot[link.name] * pipe_pitch
        heights = [y_of(name) for name in link.nodes]
        ax.plot(
            [x, x],
            [max(heights), min(heights)],
            color=color,
            linewidth=pipe_width * 72.0,
            solid_capstyle="butt",
            gid=link.name,
            zorder=2,
        )
        for y in heights:
            ax.add_patch(Circle((x, y), joint_radius, facecolor=color, linewidth=0, gid=link.name, zorder=5))

    def strand(kind: str, xs: Sequence[float], ys: Sequence[float]) -> None:
        ax.plot(
            xs,
            ys,
            color=colors[kind],
            linewidth=line_width * 72.0,
            solid_capstyle="butt",
            solid_joinstyle="miter",
            zorder=1,
        )

    for pair, kind in routes.items():
        upper, lower = pair
        strands = list(zip(pairs[pair], offsets(pair), strict=True))
        match kind:
            case "down":
                x = x_of(upper) + spine_offset(layout[upper][0])
                for link_type, offset in strands:
                    strand(link_type, [x + offset] * 2, [bottom(upper), top(lower)])
            case "bracket":
                lane_x = x_of(upper) - (lanes[pair] + 0.5) * lane_pitch
                start, end = y_of(upper) + terminal_y[(upper, pair)], y_of(lower) + terminal_y[(lower, pair)]
                for link_type, offset in strands:
                    xs = [x_of(upper), lane_x + offset, lane_x + offset, x_of(lower)]
                    strand(link_type, xs, [start - offset, start - offset, end + offset, end + offset])
            case "across":
                west, east = sorted(pair, key=lambda name: layout[name][0])
                start, end = y_of(west) + terminal_y[(west, pair)], y_of(east) + terminal_y[(east, pair)]
                middle = (right(west) + x_of(east)) / 2.0
                for link_type, offset in strands:
                    xs = [right(west), middle, middle, x_of(east)]
                    strand(link_type, xs, [start - offset, start - offset, end - offset, end - offset])
            case "elbow":
                entered, start = elbows[pair]
                from_west = layout[start][0] < layout[entered][0]
                riser = x_of(entered) + 1.5 * em if from_west else right(entered) - 1.5 * em
                origin_x = right(start) if from_west else x_of(start)
                origin_y = y_of(start) + terminal_y[(start, pair)]
                rising = y_of(entered) > y_of(start)
                end = bottom(entered) if rising else top(entered)
                sign = (1.0 if from_west else -1.0) * (1.0 if rising else -1.0)
                for link_type, offset in strands:
                    corner = riser - offset * sign
                    strand(link_type, [origin_x, corner, corner], [origin_y - offset, origin_y - offset, end])

    for k, head in enumerate(heads):
        x = x_of(head.stages[0]) + 0.8 * em + k * rod_pitch + (rod_pitch - rod_width) / 2.0
        warm = top(head.stages[0]) + gap
        for low, high in head_spans(head, warm, y_of, top, bottom):
            ax.add_patch(
                FancyBboxPatch(
                    (x, low + 0.18 * em),
                    rod_width,
                    high - low - 0.36 * em,
                    boxstyle=f"round,pad=0,rounding_size={0.25 * em}",
                    facecolor=rod,
                    linewidth=0,
                    gid=head.name,
                    zorder=2,
                )
            )
        for name in (*head.stages, *(() if head.clamp is None else (head.clamp,))):
            ax.add_patch(
                Rectangle((x, bottom(name)), rod_width, bar_height[name], facecolor=rod, linewidth=0, zorder=4)
            )

    def bounds() -> Bbox:
        to_data = ax.transData.inverted()
        artists = ax.patches + ax.lines + ax.texts
        return Bbox.union([artist.get_window_extent(renderer).transformed(to_data) for artist in artists])

    to_data = ax.transData.inverted()
    for link in sorted(regions, key=lambda link: len(link.nodes)):
        members = [node_boxes[name].get_window_extent(renderer).transformed(to_data) for name in link.nodes]
        if set(layout) <= set(link.nodes):
            members.append(bounds())
        box = Bbox.union(members).padded(1.2 * em)
        label_height = 2.0 * em
        ax.add_patch(
            FancyBboxPatch(
                (box.x0, box.y0),
                box.width,
                box.height + label_height,
                boxstyle=f"round,pad=0,rounding_size={em}",
                facecolor="none",
                edgecolor=colors[type(link).__name__],
                linewidth=1.2,
                linestyle=(0, (5, 3)),
                gid=link.name,
                zorder=0,
            )
        )
        ax.text(box.x0 + em, box.y1 + label_height / 2.0, link.name, va="center", fontproperties=label_font, color=ink)

    content = bounds()
    entries: list[tuple[Callable[[float, float], None], str]] = []
    if heads:

        def head_swatch(x: float, y: float) -> None:
            for k in range(3):
                ax.add_patch(
                    FancyBboxPatch(
                        (x + 0.25 * em, y - em + k * 0.72 * em),
                        rod_width,
                        0.55 * em,
                        boxstyle=f"round,pad=0,rounding_size={0.2 * em}",
                        facecolor=rod,
                        linewidth=0,
                    )
                )

        entries.append((head_swatch, ", ".join(head.name for head in heads) + " regenerator nodes"))
    for kind in order:
        if any(kind in kinds for kinds in pairs.values()):

            def line_swatch(x: float, y: float, color: str = colors[kind]) -> None:
                ax.plot([x, x + 1.6 * em], [y, y], color=color, linewidth=2.0, solid_capstyle="butt")

            entries.append((line_swatch, words(kind)))
    for link in pipes:

        def pipe_swatch(x: float, y: float, color: str = colors[type(link).__name__]) -> None:
            ax.plot([x + 0.8 * em] * 2, [y - 0.9 * em, y + 0.9 * em], color=color, linewidth=pipe_width * 72.0)
            ax.add_patch(Circle((x + 0.8 * em, y), joint_radius, facecolor=color, linewidth=0))

        entries.append((pipe_swatch, link.name))
    for link in regions:

        def region_swatch(x: float, y: float, color: str = colors[type(link).__name__]) -> None:
            ax.add_patch(
                Rectangle(
                    (x, y - 0.5 * em),
                    1.6 * em,
                    em,
                    facecolor="none",
                    edgecolor=color,
                    linewidth=1.2,
                    linestyle=(0, (3, 2)),
                )
            )

        entries.append((region_swatch, link.name))
    if boundary:

        def boundary_swatch(x: float, y: float) -> None:
            ax.add_patch(
                Rectangle(
                    (x, y - 0.5 * em),
                    1.6 * em,
                    em,
                    facecolor=fill,
                    edgecolor=edge,
                    hatch="////",
                    hatchcolor=hatch,
                    linewidth=0.8,
                )
            )

        entries.append((boundary_swatch, "Fixed temperature"))
    if set(sensors) & set(layout):

        def sensor_swatch(x: float, y: float) -> None:
            ax.add_patch(Rectangle((x, y - 0.5 * em), 1.6 * em, em, facecolor=fill, edgecolor=ink, linewidth=1.6))

        entries.append((sensor_swatch, "Thermometer"))
    x = content.x1 + 2.5 * em
    y = (content.y0 + content.y1) / 2.0 + (len(entries) - 1) * 1.15 * em
    for draw, label in entries:
        draw(x, y)
        ax.text(x + 2.4 * em, y, label, va="center", fontproperties=label_font, color=ink)
        y -= 2.3 * em

    fitted = bounds().padded(1.2 * em)
    fitted = Bbox.from_extents(fitted.x0, fitted.y0, fitted.x1 + 2.0 * em, fitted.y1)
    ax.set_xlim(fitted.x0, fitted.x1)
    ax.set_ylim(fitted.y0, fitted.y1)
    figure.set_size_inches(fitted.width, fitted.height)
    return figure


def link_type_order() -> list[str]:
    """Returns link class names in the order their types are registered."""
    return [inspect.signature(loader).return_annotation.__name__ for loader in LINK_TYPES.values()]


def route(pair: tuple[str, str], layout: Mapping[str, tuple[int, int]]) -> str:
    """Returns how a connector runs: across, down, bracket, or elbow."""
    (first_column, first_row), (second_column, second_row) = layout[pair[0]], layout[pair[1]]
    if first_row == second_row:
        return "across"
    if first_column != second_column:
        return "elbow"
    if any(column == first_column and first_row < row < second_row for column, row in layout.values()):
        return "bracket"
    return "down"


def elbow_ends(pair: tuple[str, str], layout: Mapping[str, tuple[int, int]]) -> tuple[str, str]:
    """Returns the node an elbow enters vertically and the node it leaves horizontally.

    The vertical leg runs in a column with no node between the pair's rows, preferring the column with fewer nodes.
    """
    upper, lower = pair
    upper_row, lower_row = layout[upper][1], layout[lower][1]

    def cost(name: str) -> tuple[bool, int]:
        column = layout[name][0]
        blocked = any(c == column and upper_row < r < lower_row for c, r in layout.values())
        return blocked, [c for c, _ in layout.values()].count(column)

    entered = min((lower, upper), key=cost)
    return entered, upper if entered == lower else lower


def interval_lanes(spans: Sequence[tuple[int, int, int]]) -> list[int]:
    """Assigns each (column, top row, bottom row) span the innermost lane free of overlapping spans in its column."""
    lanes = [0] * len(spans)
    placed: list[int] = []
    for i in sorted(range(len(spans)), key=lambda i: spans[i][2] - spans[i][1]):
        column, start, end = spans[i]
        taken = {lanes[j] for j in placed if spans[j][0] == column and min(end, spans[j][2]) >= max(start, spans[j][1])}
        lanes[i] = min(set(range(len(taken) + 1)) - taken)
        placed.append(i)
    return lanes


def head_spans(
    head: PulseTubeHead,
    warm: float,
    y_of: Callable[[str], float],
    top: Callable[[str], float],
    bottom: Callable[[str], float],
) -> list[tuple[float, float]]:
    """Returns the bottom and top height of each regenerator node along a head's rod, warm to cold."""
    first, second = head.stages
    pieces_2 = [(0.0, 1.0, bottom(first), top(second))]
    if head.clamp is not None and y_of(second) < y_of(head.clamp) < y_of(first):
        clamp = tiled_edges(head.regenerator_2.midpoint)[head.memory_node]
        pieces_2 = [(0.0, clamp, bottom(first), top(head.clamp)), (clamp, 1.0, bottom(head.clamp), top(second))]
    return [
        *node_spans(head.regenerator_1.midpoint, [(0.0, 1.0, warm, top(first))]),
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
