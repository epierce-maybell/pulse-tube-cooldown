# Data

Source measurements, hardware records, and curated physical inputs live here.
Preparation tools live in [`tools/`](../tools/README.md). Maintained validation
lives in [`simulations/`](../simulations/README.md); exploratory work and fitted
candidates live in [`../pulse-tube-research/`](../../pulse-tube-research/README.md).

| Location | Contents |
|---|---|
| `sn09/`, `fp1/` | Refrigerator cooldown records and comparison conventions, BoM exports and node assignments, and geometry. |
| `heads/pt1_1/geometry.csv`, `heads/pt1_5/geometry.csv` | Head geometry from CAD and stated estimates. PT1.5 records what its CAD shows; other quantities, including matrix filling, come from PT1.1. |
| `test_stand/<serial>/` | Measured cooldowns and load maps grouped by cold-head serial. |
| `test_stand/fit_set.csv`, `test_stand/cooldown_set.csv` | Curated stand-record selections with compressor counts, valve frequencies, source filenames, and available MLI information. |
| `test_stand/load_map_index_2026-11-05.csv` | Supplied draft index of load-map tests and setup changes. |
| `test_stand/pt_performance_overview_2026-10-30.pdf` | Supplied manufacturer overview: source for digitized maps and transcribed settings. |
| `test_stand/muenster_maps.csv` | Digitized manufacturer maps with source page and image names. |
| `test_stand/muenster_key_points.csv`, `test_stand/valve_settings.csv` | Transcribed manufacturer reference tables. |
| `test_stand/sn0117_load_map.png` | Supplied PTD 420Q (SN117) chart with two F-70H compressors; originally `load_map.png`. |
| `gases/` | Published adsorption inputs used to derive packaged gas properties. |
| `systems.md`, `test_stand.md`, `compressors.md` | Measurement formats, hardware records, operating conditions, and sources. |

SN09 is a refrigerator; S0009 is a cold-head serial. SN09's installed head serial
is unknown. Stand measurements are separate from refrigerator cooldown records.
The full stand archive is `../../data/test_stand/` relative to the repository.

The research snapshot preserves retired fitted head parameters, drive files,
the old stand configuration, and the characterization's `runs.csv` catalog.
Its `simulations/exploration/loadmaps/` contains the digitization workflow,
inspection overlays, and characterization. Historical snapshots retain their
original recorded paths and hashes.

The replacement head formulation is [in progress](../docs/DESIGN.md#open-decisions).
Full SN09 and FP1 configurations remain tasks 4.1 and 4.2.
