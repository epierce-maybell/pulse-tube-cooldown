# Validation-system records

SN09 and FP1 share the six-plate architecture: PRP, RGP, CFP, STP, ICP, and MXP. SN09 has one pulse-tube head, run on two compressors. FP1 has two PT1.5 heads on the same plates, run by two compressors as is usual for FP systems. Neither system's head serials are recorded. The run procedures are cold_control's `cooldown` routine for SN09 and `fp_cooldown` for FP1, in `../../dev/cc_routines` relative to the repository.

## Pulse-tube heads

Neither system's heads have test-stand records: FP1's predate the stand, and SN09's serial is unrecorded. Both systems use the canonical head of DESIGN "Coupled capacity maps", fitted to the Denver test-stand records of S0007, S0008, S0009, S0019, S0023, S0028, S0034, S0035, and S0038 (`test_stand.md`). The head-to-head spread among those records is accepted as error. SN09's head was recorded as a PT1.1, but no PT1.1 records exist. `test_stand/s0009/` holds S0009's test-stand records; `heads/pt1_1/geometry.csv` holds the CAD-derived geometry of SN09's PT1.1 head, and `heads/pt1_5/geometry.csv` that of FP1's PT1.5 heads.

## Cooldown records

- **SN09**: `sn09/cooldown_2026-06-04.dat`, 2026-06-04 00:01 to 2026-06-06 23:00 UTC at 60 s. Channels PRP, RGP, CFP, STP, ICP, MXP, and DHX. The compressors start at 00:21 UTC, the last sample before CFP and PRP begin to fall. Condensing begins at 2026-06-06 14:15 UTC, 61.9 h after compressor start, when STP steps down by about 2 K and ICP and DHX follow; MXP reads 14.5 K. RGP reads a constant 294.94 K until 00:29 UTC.
- **FP1**: `fp1/cooldown_2026-10-29.dat`, 2026-10-29 16:15 to 2026-11-01 06:00 UTC at 60 s, derived by `data/scripts/derive_fp1_cooldown.py` (from the repository root) as 60 s means of the InfluxDB export `FP1-10-29.csv` (the `thermometer` measurement's `temperature` field at 1 s). Channels PRP, RGP, CFP, STP, MXP, top DHX, and CMN. FP1 has no ICP thermometer. The compressors start at 16:31 UTC, the last sample before PRP and CFP begin to fall. Condensing begins at 2026-10-31 10:42 UTC, 42.2 h after compressor start, when STP steps down by about 2 K; the MXP channel reads 16.1 K. For this run the MXP channel reads a sensor on the DHX. It logs only from 2026-10-31 04:39 UTC (44.9 K). CMN is a sub-kelvin thermometer, below 0.1 K throughout, with occasional -2.0 sentinels.

Both routines begin condensing once their MXP channel falls below 14 K (`exit_temperature` in the `gas_pushing` subroutine). Compressor start and the onset of condensing are read from the traces above.

The routine logs, `monitor_ovc` P1 diagnostics, and P2 and P3 pressure traces for these runs are not pulled. The model takes its exchange-gas pressure and gas-pushing thresholds from the routines' settings.

## Inventories and geometry

`<system>/bom_by_node.csv` gives node masses by material. `<system>/geometry.csv` is read from STEP exports of the CAD assemblies by `data/scripts/step_geometry.py` (`uv run --with cadquery`, from the repository root), one row per value with its source. Supports are the legacy model's values.

## Gas-gap switches

Supplied by the human. Each system has eleven switches: PRP–RGP 2 × N2, RGP–CFP 2 × N2, CFP–STP 3 × He-4, STP–ICP 2 × He-4, and ICP–MXP 2 × He-4. All use the same core (500-000182), outer core (500-000183), and shell (500-000185). Each gap is 0.5 mm wide with 40 cm^2 of facing area. The shell is 0.50 mm 304 SS sheet, 25.0 mm OD and 56.4 mm long, giving a 39 mm^2 conduction cross-section over a 70 mm developed length, and is its own conduction link. Switches are filled to 1 atm at room temperature and sealed with a crush tube. Each holds one 0.197 cm^3 charge of granulated coconut charcoal (part 200-000021), about 0.1 g at 0.5 g/cm^3.

## Surfaces

The copper plates, covers, and shields are gold-plated (confirmed by the human) and use `Gold-plated Cu`. The aluminum PRP shields, top plate, and OVC use `Machined Al`; their finish is assumed.
