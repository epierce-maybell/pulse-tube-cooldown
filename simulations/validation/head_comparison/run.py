"""Writes a single-page population comparison as PNG and SVG.

Run from the repository root:
    python -m simulations.validation.head_comparison.run

Overlay a candidate with:
    python -m simulations.validation.head_comparison.run --candidate my_model:predict --label "My model"

The candidate function takes the complete tuple of HeadData records and returns
a mapping from head serial to Prediction (types in comparison.py). It may use a
population fit performed beforehand or fit shared parameters from these records.
The candidate owns simulation and equilibrium calculations; this script plots
its returned temperatures. Run conditions, trimming, and sources are recorded
in records.toml. Only plots are written.

Styles combine SciencePlots science, no-latex, and nature, with page.mplstyle
providing page size, panel typography, and the six-color cycle.
"""

import argparse
import importlib
from pathlib import Path

import matplotlib.pyplot as plt

from simulations.validation.head_comparison.comparison import comparison_style, load_population, plot_comparison


def main() -> None:
    """Reads the selected measurements and optional candidate, then saves the page."""
    plt.switch_backend("Agg")
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--records", type=Path, default=here / "records.toml", help="Run selection TOML")
    parser.add_argument("--candidate", help="Population prediction function as module:function")
    parser.add_argument("--label", default="Candidate", help="Candidate name printed on the page")
    parser.add_argument("--output", type=Path, default=here / "comparison", help="Output path stem")
    args = parser.parse_args()

    heads = load_population(args.records, root=here.parents[2])
    predictions = None
    if args.candidate:
        module, function = args.candidate.split(":", maxsplit=1)
        predictions = getattr(importlib.import_module(module), function)(heads)
    figure = plot_comparison(heads, predictions=predictions, candidate_name=args.label)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with comparison_style():
        for extension in ("png", "svg"):
            path = args.output.with_suffix(f".{extension}")
            figure.savefig(path, dpi=180)
            print(path)
    plt.close(figure)


if __name__ == "__main__":
    main()
