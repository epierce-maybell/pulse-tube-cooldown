"""Extracts the Münster load-map figures from the manufacturer's overview into `muenster_figures/`.

Each figure is saved as `page_<page>_<index>.png`, its page in the PDF and its order on that page. Run with pypdf
available, for example `uv run --with pypdf python extract_figures.py`.
"""

import argparse
from pathlib import Path

import pypdf

HERE = Path(__file__).parent
OVERVIEW = HERE.parents[2] / "docs" / "pt_performance_overview_2026-10-30.pdf"


def extract_figures(pdf: Path, pages: list[int], directory: Path) -> list[Path]:
    """Writes every image on the given pages of a PDF to a directory.

    Args:
        pdf: The PDF to read.
        pages: Page numbers, counting from 1.
        directory: Where the images go.

    Returns:
        The paths written.
    """
    directory.mkdir(exist_ok=True)
    reader = pypdf.PdfReader(pdf)
    written = []
    for page in pages:
        for index, image in enumerate(reader.pages[page - 1].images):
            path = directory / f"page_{page:02d}_{index}{Path(image.name).suffix}"
            path.write_bytes(image.data)
            written.append(path)
    return written


def main() -> None:
    """Extracts the figures that `digitize.toml` digitizes."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "pages", nargs="*", type=int, default=[8, 9, 11, 12, 15, 18, 21, 24, 26, 29, 32, 36, 39, 42, 45, 48, 51, 54]
    )
    for path in extract_figures(OVERVIEW, parser.parse_args().pages, HERE / "muenster_figures"):
        print(path.name)


if __name__ == "__main__":
    main()
