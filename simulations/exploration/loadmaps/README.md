# Load-map exploration

This experiment characterizes measured maps and screens candidate functional
forms. Its conclusions are in `characterization.md`; plots, digitizing overlays,
screen inputs, and reports live beside the code.

Input tables live in `data/test_stand/`: `runs.csv`, `muenster_maps.csv`,
`muenster_key_points.csv`, and `valve_settings.csv`. The supplied draft test index
is a separate source record, `data/test_stand/load_map_index_2026-11-05.csv`.
Numeric maps referenced as `archive/` in `runs.csv` come from
`../../data/test_stand/loadmap_data/` relative to the repository.

Run from the repository root:

```text
uv run --with pypdf python simulations/exploration/loadmaps/extract_figures.py
uv run python simulations/exploration/loadmaps/digitize.py
uv run python simulations/exploration/loadmaps/characterize.py
uv run python simulations/exploration/loadmaps/functional_forms.py
uv run python simulations/exploration/loadmaps/alpha_screen.py
```

`extract_figures.py` extracts images from
`docs/pt_performance_overview_2026-10-30.pdf` into ignored `muenster_figures/`.
`digitize.py` reads those images using `digitize.toml`, updates the reference
table `data/test_stand/muenster_maps.csv`, and saves inspection overlays in
`digitized/`. `characterize.py` and `functional_forms.py` write plots to
`figures/` and summaries beside their scripts. `alpha_screen.py` reads
`alpha_screen.toml` and prints JSON.
