"""Drawing of a Graph's nodes and links."""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from matplotlib.text import Text
from matplotlib.transforms import blended_transform_factory

from pulse_tube_cooldown.core.types import Graph

__all__ = ["plot_graph"]


def plot_graph(graph: Graph, ax: Axes | None = None, lane_gap: float = 12.0, edge_spacing: float = 6.0) -> Axes:
    """Draws each node and each link as edges to its nodes.

    Nodes form a column in graph order, first on top. Links run in vertical lanes beside the column, with right-angled
    edges to their nodes; each edge line has its link's name as its `gid`.

    Args:
        graph: The graph to draw.
        ax: Axes to draw on. If omitted, draws on a new figure sized to fit.
        lane_gap: Horizontal gap in points between the column and the lanes that hold links, and between lanes.
        edge_spacing: Vertical distance in points between edges where they meet a node.

    Returns:
        The Axes drawn on, with nodes labelled by name and links by name and class name. Horizontal positions are in
        axes fractions and vertical positions in node spacings, with the first node at 0.
    """
    if ax is None:
        fig = plt.figure(figsize=(6.4, max(2.0, 0.8 * len(graph.nodes))))
        ax = fig.add_axes((0.0, 0.0, 1.0, 1.0))
        fit_figure = True
    else:
        fit_figure = False
    transform = blended_transform_factory(ax.transAxes, ax.transData)
    renderer = ax.figure.canvas.get_renderer()

    def label(text: str, y: float, color: str, facecolor: str, pad: float, fontsize: str) -> tuple[Text, float]:
        """Adds a boxed label and returns it with its width in points."""
        artist = ax.text(
            0.0,
            y,
            text,
            transform=transform,
            ha="center",
            va="center",
            color="#0b0b0b",
            fontsize=fontsize,
            bbox={"boxstyle": f"round,pad={pad}", "facecolor": facecolor, "edgecolor": color},
            zorder=5 if fontsize == "medium" else 4,
        )
        text_points = artist.get_window_extent(renderer).width * 72.0 / ax.figure.dpi
        return artist, text_points + 2.0 * pad * artist.get_fontsize()

    node_y = -np.arange(len(graph.nodes), dtype=np.float64)
    node_labels = [
        label(node.name, y, "#52514e", "#f0efec", 0.4, "medium") for node, y in zip(graph.nodes, node_y, strict=True)
    ]
    half_node = max((width for _, width in node_labels), default=0.0) / 2.0

    palette = ("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948")
    class_names = [type(link).__name__ for link in graph.links]
    colors = {name: palette[i % len(palette)] for i, name in enumerate(dict.fromkeys(class_names))}
    lanes = _assign_lanes(graph.connectivity)
    link_labels = [
        label(f"{link.name}\n{name}", _label_y(indices), colors[name], "white", 0.2, "small")
        for link, name, indices in zip(graph.links, class_names, graph.connectivity, strict=True)
    ]
    sides = [-1.0 if len(indices) == 1 else 1.0 for indices in graph.connectivity]

    lane_widths: dict[tuple[float, int], float] = {}
    for side, lane, (_, width) in zip(sides, lanes, link_labels, strict=True):
        lane_widths[side, lane] = max(lane_widths.get((side, lane), 0.0), width)

    def lane_center(side: float, lane: int) -> float:
        """Distance in points from the column to a lane's center, one gap past each lane inside it."""
        inside = sum(lane_widths[side, inner] + lane_gap for inner in range(lane))
        return half_node + lane_gap + inside + lane_widths[side, lane] / 2.0

    def extent(side: float) -> float:
        outer = max((lane for lane_side, lane in lane_widths if lane_side == side), default=None)
        return half_node if outer is None else lane_center(side, outer) + lane_widths[side, outer] / 2.0

    column_x = lane_gap + extent(-1.0)
    if fit_figure:
        ax.figure.set_figwidth((column_x + extent(1.0) + lane_gap) / 72.0)

    ax.set_ylim(min(node_y, default=0.0) - 0.6, 0.6)
    axes_box = ax.get_window_extent(renderer)
    points_to_fraction = ax.figure.dpi / 72.0 / axes_box.width
    points_to_y = np.ptp(ax.get_ylim()) * ax.figure.dpi / 72.0 / axes_box.height
    lane_x = [
        (column_x + side * lane_center(side, lane)) * points_to_fraction
        for side, lane in zip(sides, lanes, strict=True)
    ]
    column_x *= points_to_fraction
    offsets = _attachment_offsets(graph.connectivity, lane_x, edge_spacing * points_to_y)

    for artist, _ in node_labels:
        artist.set_x(column_x)

    stub_heights: dict[int, list[tuple[int, float]]] = {}
    for indices, lane, (artist, _) in zip(graph.connectivity, lanes, link_labels, strict=True):
        if len(indices) == 1:
            box_points = artist.get_window_extent(renderer).height * 72.0 / ax.figure.dpi + 0.4 * artist.get_fontsize()
            stub_heights.setdefault(indices[0], []).append((lane, box_points / 2.0))

    for i, (link, indices, (artist, _)) in enumerate(zip(graph.links, graph.connectivity, link_labels, strict=True)):
        color = colors[class_names[i]]
        label_xy = (lane_x[i], artist.get_position()[1])
        artist.set_x(lane_x[i])
        if len(indices) == 1:
            # Outer one-node links pass above the labels of the inner ones on the same node, then drop into their own.
            y = node_y[indices[0]]
            inner = [half_height for lane, half_height in stub_heights[indices[0]] if lane < lanes[i]]
            rise = (max(inner) + edge_spacing * (lanes[i] + 1)) * points_to_y if inner else 0.0
            ax.plot(
                [column_x, column_x, lane_x[i], lane_x[i]],
                [y, y + rise, y + rise, y],
                color=color,
                linewidth=1.2,
                gid=link.name,
                transform=transform,
            )
            continue

        routes = []
        for node in indices:
            attach = node_y[node] + offsets[i, node]
            routes.append(([lane_x[i], lane_x[i], column_x, column_x], [label_xy[1], attach, attach, node_y[node]]))
        if len(indices) == 2:
            (first_x, first_y), (second_x, second_y) = routes
            routes = [(first_x[::-1] + second_x[1:], first_y[::-1] + second_y[1:])]
        for xs, ys in routes:
            ax.plot(xs, ys, color=color, linewidth=1.2, gid=link.name, transform=transform)
        if len(indices) > 2:
            junctions = [node_y[node] + offsets[i, node] for node in indices]
            ax.plot([lane_x[i]] * len(junctions), junctions, "o", color=color, markersize=4, transform=transform)

    ax.set_axis_off()
    return ax


def _label_y(indices: tuple[int, ...]) -> float:
    """Height of a link's label: midway between two nodes of its span, so no node's edges pass behind it."""
    low, high = min(indices), max(indices)
    if low == high:
        return -float(low)
    return -(np.floor((low + high) / 2.0) + 0.5)


def _attachment_offsets(
    connectivity: tuple[tuple[int, ...], ...], lane_x: list[float], spacing: float
) -> dict[tuple[int, int], float]:
    """Vertical offset of each link's edge at each of its nodes, keyed by link and node index.

    Edges meeting at a node are stacked so they neither overlap nor cross. From bottom to top: edges running down to
    their lane, inner lanes lowest; then edges running up, inner lanes highest.
    """
    at_node: dict[int, list[int]] = {}
    for i, indices in enumerate(connectivity):
        if len(indices) > 1:
            for node in indices:
                at_node.setdefault(node, []).append(i)

    offsets = {}
    for node, links in at_node.items():
        runs_up = {i: _label_y(connectivity[i]) > -node for i in links}
        links.sort(key=lambda i: (1, -lane_x[i]) if runs_up[i] else (0, lane_x[i]))
        for position, i in enumerate(links):
            offsets[i, node] = (position - (len(links) - 1) / 2.0) * spacing
    return offsets


def _assign_lanes(connectivity: tuple[tuple[int, ...], ...]) -> list[int]:
    """Gives each link the innermost lane free of the links it would collide with, shortest span first.

    Links on two or more nodes collide when their node ranges overlap; ranges that only share an end node do not.
    Links on one node collide only with other one-node links on the same node.
    """
    spans = [(min(indices), max(indices)) for indices in connectivity]

    def collide(a: tuple[int, int], b: tuple[int, int]) -> bool:
        if a[0] == a[1] or b[0] == b[1]:
            return a == b
        return min(a[1], b[1]) > max(a[0], b[0])

    lanes = [0] * len(spans)
    placed: list[int] = []
    for i in sorted(range(len(spans)), key=lambda i: spans[i][1] - spans[i][0]):
        taken = {lanes[j] for j in placed if collide(spans[i], spans[j])}
        lanes[i] = min(set(range(len(taken) + 1)) - taken)
        placed.append(i)
    return lanes
