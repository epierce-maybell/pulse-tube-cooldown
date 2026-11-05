"""Digitizes the Münster load-map figures listed in `digitize.toml` into `muenster_maps.csv`.

Each axis is calibrated against the figure's own gridlines, starting from two reference gridlines and refined by a
least-squares fit to every gridline found. Markers are found by color, and each marker's position is the center of
its bounding box. Grid markers are assigned heater powers by order: markers split into one column per stage-1 power
at the widest gaps in T1, and within a column they take the stage-2 powers in order of T2. An overlay of every
figure's markers and assignments is written to `digitized/` for inspection.
"""

import csv
import tomllib
from dataclasses import dataclass
from pathlib import Path

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import NDArray
from scipy import ndimage

HERE = Path(__file__).parent
DATA = HERE.resolve().parents[2] / "data/test_stand"


@dataclass(frozen=True, slots=True)
class AxisMap:
    """A linear map from pixel position to data value along one axis.

    Attributes:
        intercept: Data value at pixel 0.
        slope: Data units per pixel.
        residual: Root-mean-square distance in pixels of the fitted gridlines from their values.
        gridlines: Number of gridlines in the fit.
    """

    intercept: float
    slope: float
    residual: float
    gridlines: int

    def __call__(self, pixel: NDArray[np.float64]) -> NDArray[np.float64]:
        """Returns the data values at pixel positions."""
        return self.intercept + self.slope * pixel


def read_image(path: Path, crop: list[int] | None = None) -> NDArray[np.float64]:
    """Reads a figure as RGB values from 0 to 1, keeping pixel columns [start, stop) if cropped."""
    image = mpimg.imread(path)[..., :3].astype(np.float64)
    return image[:, crop[0] : crop[1]] if crop else image


def find_gridlines(image: NDArray[np.float64], axis: int, coverage: float = 0.6) -> NDArray[np.float64]:
    """Finds the pixel positions of gray lines spanning the figure.

    Args:
        image: RGB figure.
        axis: 1 for vertical lines, positions along x; 0 for horizontal lines, positions along y.
        coverage: Fraction of the figure's extent a line must cover.

    Returns:
        The center of each run of adjacent line pixels.
    """
    brightest, dimmest = image.max(axis=2), image.min(axis=2)
    gray = (brightest < 0.93) & (brightest - dimmest < 0.06)
    pixels = np.flatnonzero(gray.mean(axis=1 - axis) > coverage)
    runs = np.split(pixels, np.flatnonzero(np.diff(pixels) > 2) + 1)
    return np.array([run.mean() for run in runs if run.size])


def calibrate(gridlines: NDArray[np.float64], reference: list[list[float]], step: float) -> AxisMap:
    """Fits an axis map to the gridlines that fall on multiples of the gridline step.

    Args:
        gridlines: Pixel positions of the figure's gridlines along the axis.
        reference: Two gridlines as [pixel, value], which set the starting map.
        step: Spacing of the finest gridlines in data units.

    Returns:
        The map fitted to every gridline within a sixth of a step of a multiple of the step under the starting map.
    """
    (pixel_a, value_a), (pixel_b, value_b) = reference
    slope = (value_b - value_a) / (pixel_b - pixel_a)
    guess = value_a + slope * (gridlines - pixel_a)
    nearest = np.round(guess / step) * step
    on_grid = np.abs(guess - nearest) < step / 6
    slope, intercept = np.polyfit(gridlines[on_grid], nearest[on_grid], 1)
    residual = np.sqrt(np.mean(((intercept + slope * gridlines[on_grid] - nearest[on_grid]) / slope) ** 2))
    return AxisMap(float(intercept), float(slope), float(residual), int(on_grid.sum()))


def find_markers(image: NDArray[np.float64], color: str, minimum_area: int = 12) -> NDArray[np.float64]:
    """Finds filled markers of one color.

    Args:
        image: RGB figure.
        color: "red", "blue", or "yellow".
        minimum_area: Smallest marker in pixels, after lines thinner than three pixels are removed.

    Returns:
        The bounding-box center of each marker as (x, y) pixel positions.
    """
    red, green, blue = image[..., 0], image[..., 1], image[..., 2]
    mask = {
        "red": (red > 0.75) & (green < 0.35) & (blue < 0.35),
        "blue": (blue > 0.75) & (red < 0.35) & (green < 0.35),
        "yellow": (red > 0.8) & (green > 0.8) & (blue < 0.35),
    }[color]
    labels, _ = ndimage.label(ndimage.binary_opening(mask, structure=np.ones((3, 3))))
    centers = []
    for index, region in enumerate(ndimage.find_objects(labels), start=1):
        if (labels[region] == index).sum() >= minimum_area:
            rows, columns = region
            centers.append(((columns.start + columns.stop - 1) / 2, (rows.start + rows.stop - 1) / 2))
    return np.array(centers)


def assign_grid(t1: NDArray[np.float64], t2: NDArray[np.float64], q1: list[float], q2: list[float]) -> list[tuple]:
    """Assigns grid markers their heater powers.

    Args:
        t1: Stage-1 temperature of each marker.
        t2: Stage-2 temperature of each marker.
        q1: Stage-1 powers, low to high.
        q2: Stage-2 powers, low to high.

    Returns:
        (Q1, Q2, T1, T2) for every marker.
    """
    order = np.argsort(t1)
    gaps = np.sort(np.argsort(np.diff(t1[order]))[::-1][: len(q1) - 1])
    columns = np.split(order, gaps + 1)
    points = []
    for power_1, column in zip(q1, columns, strict=True):
        if column.size != len(q2):
            raise ValueError(f"column at {power_1} W has {column.size} markers, expected {len(q2)}")
        for power_2, index in zip(q2, column[np.argsort(t2[column])], strict=True):
            points.append((power_1, power_2, float(t1[index]), float(t2[index])))
    return points


def digitize_figure(figure: dict, overlay: Path) -> tuple[list[tuple], AxisMap, AxisMap]:
    """Digitizes one figure and draws its overlay.

    Args:
        figure: One `[[figure]]` entry of `digitize.toml`.
        overlay: Path for the overlay image.

    Returns:
        The (Q1, Q2, T1, T2) points, with T1 NaN for a sweep, and the x and y axis maps.
    """
    image = read_image(HERE / "muenster_figures" / figure["image"], figure.get("crop"))
    x_map = calibrate(find_gridlines(image, axis=1), figure["x_reference"], figure["x_step"])
    y_map = calibrate(find_gridlines(image, axis=0), figure["y_reference"], figure["y_step"])
    markers = find_markers(image, figure["marker"])
    x, y = x_map(markers[:, 0]), y_map(markers[:, 1])
    keep = np.ones(x.size, dtype=bool)
    for x_low, x_high, y_low, y_high in figure.get("ignore", []):
        keep &= ~((x >= x_low) & (x <= x_high) & (y >= y_low) & (y <= y_high))
    x, y, markers = x[keep], y[keep], markers[keep]
    if figure.get("kind", "grid") == "sweep":
        order = np.argsort(x)
        if order.size != len(figure["q2"]):
            raise ValueError(f"sweep has {order.size} markers, expected {len(figure['q2'])}")
        points = [
            (figure["q1"][0], power_2, np.nan, float(y[i])) for power_2, i in zip(figure["q2"], order, strict=True)
        ]
    else:
        points = assign_grid(x, y, figure["q1"], figure["q2"])

    fig, ax = plt.subplots(figsize=(image.shape[1] / 100, image.shape[0] / 100), dpi=100)
    ax.imshow(image)
    ax.scatter(markers[:, 0], markers[:, 1], s=120, facecolors="none", edgecolors="#1baf7a", linewidths=1.5)
    for power_1, power_2, t1, t2 in points:
        pixel_x = (
            x[np.argmin(np.abs(y - t2) + (0 if np.isnan(t1) else np.abs(x - t1)))] - x_map.intercept
        ) / x_map.slope
        pixel_y = (t2 - y_map.intercept) / y_map.slope
        ax.annotate(
            f"{power_1:g}/{power_2:g}",
            (pixel_x, pixel_y),
            xytext=(6, -10),
            textcoords="offset points",
            fontsize=8,
            color="#e34948",
        )
    ax.set_axis_off()
    fig.subplots_adjust(0, 0, 1, 1)
    fig.savefig(overlay)
    plt.close(fig)
    return points, x_map, y_map


def main() -> None:
    """Digitizes every figure, writes `muenster_maps.csv`, and prints each calibration."""
    figures = tomllib.loads((HERE / "digitize.toml").read_text(encoding="utf-8"))["figure"]
    (HERE / "digitized").mkdir(exist_ok=True)
    rows = []
    for figure in figures:
        kind = figure.get("kind", "grid")
        name = f"{figure['head']}_{figure['date']}_{kind}"
        points, x_map, y_map = digitize_figure(figure, HERE / "digitized" / f"{name}.png")
        page = int(figure["image"].split("_")[1])
        print(
            f"{name:28s} page {page:2d}  x: {x_map.gridlines:2d} lines, {x_map.residual:.2f} px rms  "
            f"y: {y_map.gridlines:2d} lines, {y_map.residual:.2f} px rms  "
            f"resolution {x_map.slope:.3f} K/px, {-y_map.slope:.4f} K/px"
        )
        for power_1, power_2, t1, t2 in points:
            rows.append(
                {
                    "head": figure["head"],
                    "date": str(figure["date"]),
                    "kind": kind,
                    "page": page,
                    "image": figure["image"],
                    "Q1": power_1,
                    "Q2": power_2,
                    "T1": "" if np.isnan(t1) else f"{t1:.3f}",
                    "T2": f"{t2:.4f}",
                    "inferred_levels": bool(figure.get("inferred")),
                }
            )
    with open(DATA / "muenster_maps.csv", "w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
