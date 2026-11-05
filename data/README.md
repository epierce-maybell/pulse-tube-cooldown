# Data

This directory holds source measurements, reference tables, and inputs to the
simulator. Preparation tools live in [`data/scripts/`](scripts/README.md).
Analysis code, figures, fitted candidates, and reports live with their experiment
under `simulations/`. Published documents and diagrams live under `docs/`.

| Location | Contents |
|---|---|
| `sn09/`, `fp1/` | Refrigerator cooldown records and comparison conventions, BoM exports and node assignments, geometry, and fridge configuration templates. |
| `heads/pt1_1/geometry.csv`, `heads/pt1_5/geometry.csv` | Pulse-tube head geometry by model, from CAD and stated estimates. PT1.5 records what its CAD shows; its other quantities, including the matrix filling, are PT1.1's. |
| `test_stand/<serial>/` | Measured cooldowns and load maps, grouped by cold-head serial. |
| `test_stand/fit_set.csv` | The canonical head's load-map fit set, one map per head, with each run's compressors, valve frequency, and MLI state. |
| `test_stand/load_map_index_2026-11-05.csv` | Supplied November 5, 2026 draft index of load-map tests and setup changes. |
| `test_stand/runs.csv` | Catalog of the datasets used in the load-map characterization, with provenance and assumptions. |
| `test_stand/muenster_maps.csv` | Digitized manufacturer load maps, with source page and image names. |
| `test_stand/muenster_key_points.csv`, `test_stand/valve_settings.csv` | Transcribed manufacturer reference tables. |
| `test_stand/sn0117_load_map.png` | Supplied PTD 420Q (SN117) chart with two F-70H compressors; original filename `load_map.png`. |
| `gases/` | Published adsorption inputs used to derive the library's gas properties. |
| `systems.md`, `test_stand.md`, `compressors.md` | Measurement formats, hardware records, operating conditions, and sources. |

SN09 is a refrigerator; S0009 is a cold-head serial. SN09's installed head serial
is unknown. The test-stand measurements are separate from the
refrigerator cooldown records.

The full test-stand archive is `../../data/test_stand/` relative to the repository.
The characterization and digitizing workflow is in
[`simulations/exploration/loadmaps/`](../simulations/exploration/loadmaps/README.md).
Existing simulation snapshots retain the paths and hashes recorded when they ran.
