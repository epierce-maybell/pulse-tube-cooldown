# Pulse-tube test-stand records

The test stand runs one pulse-tube head at a time and logs three Lake Shore Cernox thermometers: T1 on stage 1, T2 on stage 2, and T3 on the pulse tube itself. The raw files label T3's column `TPR`; the readers retain that column name. Each `test_stand/<serial>/` folder holds that head's cooldown traces (`cooldown_<date>.dat`) and load maps (`load_map_<date>.dat`). The stand's full archive, 43 continuous logs and 22 load maps across eleven heads from March to October 2026, is in `../../data/test_stand/` relative to the repository.

The supplied draft test index is `test_stand/load_map_index_2026-11-05.csv`. Manufacturer reference tables are alongside it. The characterization's dataset catalog, code, and reports are preserved under `../pulse-tube-research/` in their original `data/test_stand/runs.csv` and `simulations/exploration/loadmaps/` locations.

The human confirmed T3's mounting on 2026-11-02: on the main pulse tube reaching stage 2, approximately 2 inches (0.0508 m) below stage 1. This approximate position is recorded in `heads/pt1_1/geometry.csv`. Its attachment details and sensor serial remain unspecified. The model's regenerator clamp temperature is a separate quantity, tested by the SN09 system's RGP trace.

## Inventory

Supplied by the human. Stage 1: the head's upper cold-end heat exchanger, 392.8 g of copper (43.82 cm^3); the flange, 453.6 g of aluminum, taken as 6061-T6; the radiation shield, 1.191 kg of copper; and 70.9 g of mounting hardware, taken as stainless fasteners. Stage 2: the head's lower cold-end heat exchanger, 421.8 g of 304 stainless (52.55 cm^3, as assigned in CAD 350-000162); and 39.7 g of copper mounting hardware. T3's mass is negligible.

The corrected inventory is recorded in `test_stand/cold_end_inventory.toml`. The modeled stand has 4.070932 kg of copper, 0.4536 kg of aluminum, and 0.0709 kg of stainless on stage 1; stage 2 has 0.946218 kg of copper and 0.4218 kg of stainless. The assembly STEP identifies additional lower blocks and lids with volumes of 2.77582e-4 m^3 on stage 1 and 1.01174e-4 m^3 on stage 2. They are provisionally assigned OFHC copper, adding 2.487132 kg and 0.906518 kg respectively. The original upper-component material assignments are retained. Each installed head includes 2.879932 kg of stage-1 copper and 0.906518 kg of stage-2 copper plus 0.4218 kg of stage-2 stainless. These cold-end parts are assumed to share their attached stage temperatures, with the same inventory provisionally used for PT1.1 and PT1.5. The stand flange, shield, and mounting hardware remain stand loads.

MLI wraps the two regenerators, pulse tube 1, and both stages in 4 to 12 layers each, varying per run. S0009's Münster run had 7 on each regenerator and pulse tube 1, 4 on stage 1, and 5 on stage 2; the Denver runs are assumed the same. Layer density is unmeasured and assumed 2000 layers/m, the low end of the 20 to 30 layers/cm typical of installed blankets. Across 10 to 30 layers/cm, the flux through 4 layers from 300 K to 45 K changes by about a factor of 3.5. The wrapping is looser than a true MLI blanket; per the human on 2026-11-07, it has changed the stage loads at base by only about 0.1 W.

## Format

A cooldown file opens with the calibration curve for each channel, T2 through Z1512 and T1 and T3 (`T3/Taux` in the header) both through Z1518, then a column line. These are the header's recorded assignments; matching curves to the installed sensor serials remains unverified. Each row holds the timestamp, elapsed minutes, T1, TPR, T2, two copies of T1 (`Taux` and `Taux2`), and unused fields. Rows come every 2 s, and T3 updates about every 15 s. The elapsed column can restart within a file, so the timestamp is the clock. `test_stand/s0009/cooldown_2026-05-14.dat` has a second short run appended from 14:19:17. The 2026-03-27 10:01:05 and 2026-04-13 09:12:23 logs write the timestamp as ISO 8601 (`2026-04-13T09:12:23`) under a `Date TimeElapsed(m) ...` column line. The archive's day splits (`..._Day_1.dat`) and S0019's 2026-08-06 log have no header.

A load map opens with the valve settings, the valve frequency (`f:1.1`), and the applied powers, then holds one row per steady state: stage-1 and stage-2 heater power in W, elapsed minutes, T1, TPR, T2, `Taux`, and `Taux2`. The operator types in the `f:` field: S0007's header matches its shipped setting, and S0009's and S0034's do not. A run's valve frequency is the one the index states for it, otherwise its head's shipped setting in `data/test_stand/pt_performance_overview_2026-10-30.pdf`, tabulated in `test_stand/valve_settings.csv`.

Some load maps exist only as steady plateaus in a continuous log, about 2 h per point: S0035's on 2026-09-25, and S0038's on 2026-10-02, -06, -08, -09, and -12. A continuous log's trailing fields do not hold the heater powers; they read 0 throughout the Sep 21 log, which contains the plateaus of `dataLoadMap Mon Sep 21 16_36_05 2026.dat`. That file's powers are 25, 20, and 0 W on stage 1 against 0.5 and 0 W on stage 2, and its points match S0035's Sep 25 plateaus within 0.2 K, so it is taken as S0035's.

## Load-map provenance

The maintained head comparison selects maps and records their conditions in
`simulations/validation/head_comparison/records.toml` relative to the repository.
The earlier population-fit catalog is preserved in the research archive.

Two archive files carry other names than the index's. `LoadMap-SN023-06-20-2026.dat` is S0023's 2026-06-19 map, `dataLoadMap Fri Jun 19 10_52_14 2026.dat`: its first point matches the plateau at 13:16 in that day's continuous log. `dataLoadMap Mon Sep 21 16_36_05 2026.dat` is S0035's 2026-09-21 map. The index lists no raw data for S0019's, S0034's, and S0035's selected maps; the archive holds each one's stand file.

## Cooldown comparison windows

The maintained comparison records each selected cooldown and its window in
`simulations/validation/head_comparison/records.toml`. The original cooldown
catalog is preserved in the research archive.

Each copy holds its archive log's header and one contiguous run of its rows, unchanged. The compressor starts at the first row where T2 falls more than 1 K from the row before. A copy starts at the log's first row, except S0038's, which starts 10 min before its compressor starts since its log holds earlier runs. It ends at the last row before the first heater step, or 4 h after the compressor starts, whichever comes first. S0009's copy is the whole log and so ends with 1 min of the next load map's first heater step, from 14:20:39.

| Serial | Compressor start | Copy ends | Archive lines |
|---|---|---|---|
| S0007 | 09:12:52 | 11:56:33, the log's end | all 4919 |
| S0009 | 10:25:21 | 14:21:51, the log's end | all 7349 |
| S0019 | 16:17:15 | 20:17:15, at 4 h | 1 to 7232 |
| S0028 | 17:12:03 | 21:12:02, at 4 h | 1 to 8213 |
| S0035 | 11:24:06 | 14:00:02, before the heater step at 14:00:04 | 1 to 5073 |
| S0038 | 14:44:33 | 17:00:02, before the heater step at 17:00:04 | 1 to 4, and 127186 to 131535 |

## Channel behavior

Found across twelve cooldowns of at least six heads, all logged with the same curves.

- **T1** reads NaN between 70 and 76 K and between 219 and 230 K. Above 230 K it latches at a stale value, 259 to 292 K at rest.
- **T2** reads NaN between 70 and 76 K, and between 20 and 21 K on some heads, SN09 included. At rest it reads 306 to 310 K while TPR reads 292 to 296 K.
- **T3 (logged as TPR)** reads 33 to 35 K at base on every head, and SN09's two load maps agree within 0.7 K at matching loads. During cooldown it rises in a pattern that repeats on a given head and differs between heads: above room temperature on SN09 (325 and 351 K), SN038 (365 and 376 K), and SN035 (385 K), and below it on SN007, SN019, and SN028. It then drops to base within about 35 min once T1 reaches 45 to 60 K. This is the pulse-tube sensor's trace; its thermal behavior and calibration assignment remain to be investigated.

## Use

Load-map fits use T1 and T2 at each point. Cooldown comparisons use T1 below 219 K and T2 throughout. T3 is an unscored pulse-tube diagnostic, plotted separately from the model's clamp temperature. The present network has no thermal node corresponding to T3's mounting point. The first-pass experiment uses T3's initial reading as an assumed room-temperature proxy for the initial state and warm flange; that boundary temperature has no direct measurement in these records.

SN09's May 14 runs are on one F-70 at 1.0 Hz: Münster recorded S0009 that way, within a kelvin of the May points. The August runs' compressor count is unrecorded and taken as two, since their lower T1 and higher zero-load T2 match S0028's change from one compressor to two.
