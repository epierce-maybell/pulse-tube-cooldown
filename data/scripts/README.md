# Data preparation

Run these tools from the repository root. Their inputs and outputs live in
`data/`, except for the head-model file in the installed library.

| Tool | Inputs and outputs |
|---|---|
| `reduce_bom.py` | Reads each system's BoM export; writes `bom_by_node.csv` and `bom_line_assignments.csv` beside it. |
| `step_geometry.py` | Extracts a STEP assembly into an inventory and BREP parts, then prints a geometry report for transcription into `data/<system>/geometry.csv`. Requires CadQuery. |
| `derive_fp1_cooldown.py` | Bins an InfluxDB thermometer CSV export into 60 s means; writes `data/fp1/cooldown_<date>.dat` or a supplied output path. Requires pandas. |
| `derive_head_model.py` | Reads `data/sn09/pulse_tube/geometry.csv`; writes `src/pulse_tube_cooldown/library/heads/data/pt1_1.toml`. |
| `derive_he4_adsorption.py` | Reads `data/gases/martins2011.toml`; prints fitted Langmuir constants and residuals. |

Examples:

```text
uv run python data/scripts/reduce_bom.py
uv run --with cadquery python data/scripts/step_geometry.py extract SN09 <assembly.step> <work_directory>
uv run --with cadquery python data/scripts/step_geometry.py analyze SN09 <work_directory>
uv run --with pandas python data/scripts/derive_fp1_cooldown.py <export.csv> [<output.dat>]
uv run python data/scripts/derive_head_model.py
uv run python data/scripts/derive_he4_adsorption.py
```

CAD extraction inventories and BREP files belong in the analysis directory for
the geometry investigation. Curated geometry values and their sources belong in
the system's `geometry.csv`.
