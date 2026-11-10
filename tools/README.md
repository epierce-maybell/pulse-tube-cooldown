# Data preparation

Run these tools from the repository root. Their inputs and outputs live in
`data/`, except for packaged head-design files under `src/`.

| Tool | Inputs and outputs |
|---|---|
| `reduce_bom.py` | Reads each system's BoM export; writes `bom_by_node.csv` and `bom_line_assignments.csv` beside it. |
| `step_geometry.py` | Extracts a STEP assembly into an inventory and BREP parts, then prints a geometry report for transcription into `data/<system>/geometry.csv`. Requires CadQuery. |
| `derive_fp1_cooldown.py` | Bins an InfluxDB thermometer CSV export into 60 s means; writes `data/fp1/cooldown_<date>.dat` or a supplied output path. Requires pandas. |
| `derive_head_model.py` | Reads `data/heads/pt1_1/geometry.csv` and `data/heads/pt1_5/geometry.csv`; writes `pt1_1.toml` and `pt1_5.toml` in `src/pulse_tube_cooldown/library/heads/data/`. |
| `derive_he4_adsorption.py` | Reads `data/gases/martins2011.toml`; prints fitted Langmuir constants and residuals. |

Examples:

```text
uv run python tools/reduce_bom.py
uv run --with cadquery python tools/step_geometry.py extract SN09 <assembly.step> <work_directory>
uv run --with cadquery python tools/step_geometry.py analyze SN09 <work_directory>
uv run --with pandas python tools/derive_fp1_cooldown.py <export.csv> [<output.dat>]
uv run python tools/derive_head_model.py
uv run python tools/derive_he4_adsorption.py
```

CAD extraction inventories and BREP files belong with their geometry
investigation under `../pulse-tube-research/`. Curated geometry values and their
sources belong in the system's `geometry.csv`.
