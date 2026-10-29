# Tasks

Build the v1 simulator in `docs/DESIGN.md`: the six-plate systems SN09 and FP1 from room temperature to the start of condensing, matching their measured cooldowns at least as well as the legacy model does. "As well" means the least-squares residual in log temperature across the system thermometers other than MXP.

## Reading this list

DESIGN.md is the specification and `docs/adrs/` holds the decisions behind it. AGENTS.md says how we work. The legacy code at `../pulse_tube_cooldown_simulation`, branch `Advanced_Feature_SIM`, is the reference for physics and data. Only the human edits DESIGN; when a discussion settles something, propose the wording.

Pick any unchecked task whose prerequisites are done. Prerequisites are the order within a section and whatever the entry names. Mark a task `[o]` while you have it and `[x]` when it lands. Refer to tasks by number; a new task takes the next free number in its section. Each entry says what the piece is and what finished looks like. How to build and test it is up to you. If work turns up that isn't here, add it.

Three kinds of entry:

- **Auto**: enough is settled to build it and hand it back for review.
- **Discuss**: start by talking it through with the human. A short discussion rolls straight into implementation; a longer one ends in a DESIGN update and the task becomes Auto.
- **Human**: information only the human can get.

## Done

The first vertical slice is complete, plus a few links beyond it:

- Core types, graph, procedure, and model; the engine, from heat balances through phase switching to a Result as an xarray Dataset.
- Material property fits and loader with OFHC copper RRR100 and 304 stainless; the Surface type and loader; the gas loader with He-4 and N2.
- Configuration loading, `simulate_cooldown`, `save_output` to netCDF, the `run-ptc` CLI, and the toy run in `simulations/validation/toy/` with its integration test.
- Conduction and radiation links, the gas-pushing-begins and condensing-begins conditions, and the temperature plot.
- Decisions: validate on SN09 and FP1; the cooler model and the gas-pushing model (ADRs); gas-gap switch and radiation physics (DESIGN); exchange gas (recorded in 1.3, not yet in DESIGN).

Where DESIGN is silent, the code and its tests are the record.

## 1. Links

- [ ] **1.1 Gas-pushing link** · Auto
  The `gas_pushing_pm5` model from the gas-pushing ADR, in `library/links/gas_pushing.py`: conductances between adjacent nodes in the chain from the downstream volumes and temperatures, summed into signed node flows that conserve energy. Register it so the DESIGN configuration example loads.

- [ ] **1.2 Gas-gap switch link** · Auto
  The reduced-order switch in DESIGN "Gas-gap heat switch": pressure from the adsorption balance, capped at the cold node's vapor pressure, with the brackets in series. Finished when a configured run shows the switch carrying heat early and shutting itself off as it cools, with no phase change. The He-4 adsorption constants are a pending input (5.6); build against placeholders and swap them in when they arrive. The real hardware is described at the end of section 5.

- [ ] **1.3 Exchange-gas link** · Auto
  Decided in discussion, not yet in DESIGN. One multi-node link that sums the switch's gap conductance over facing surface pairs, each pair giving `nodes`, `area`, and `gap`, including vacuum can to PRP, with N2 standing in for air. Pressure is the `monitor_ovc` P1(T_CFP) cubic from `../../dev/cc_routines`, in mbar rather than legacy's torr and clamped at zero, so the link shuts off near 68 K on its own. Share the gap conductance with the switch rather than duplicating it. Flows across the pairs sum to zero in both the viscous and free-molecular limits. Propose the DESIGN wording while you're at it.

- [ ] **1.4 Multi-layer insulation** · Discuss
  How MLI modifies radiative transfer, what it's parameterized by, and whether it's a separate link or an option on radiation. Legacy has MLI functions to start from. Then build it.

- [ ] **1.5 Joint contact** · Discuss
  Contact conductance at bolted joints and the parameters or measurements it needs. Legacy has `thermal_conductance_calc.py` and `conductance_systems.json`. Then build it.

## 2. Cooler

The model is `docs/adrs/coupled-capacity-maps.md`. The milestone is reproducing the bare head's cooldown on the test stand. Registering the cooler as a configurable link waits on that.

- [ ] **2.1 Helium properties** · Auto
  Add `he34pak` (`../../dev/he34pak`; pure He-4 is its `hepak` subpackage) as a dependency and expose He-4 density, isobaric specific heat, and isobaric expansion coefficient as array functions of temperature and pressure over 2 K to 300 K at the head's pressures. Check the near-ideal limit at 300 K and the real-gas regime at 4 K.

- [ ] **2.2 Head model** · Auto
  The ADR's governing equations in `library/cooler/model.py`, for one compressor driving one or more heads. Solve the drive constraint for the pressure swing under its relief cap, then the mass flows, utilizations, losses, enthalpy flows, and tube conduction. It returns the cooling at stage 1, the regenerator clamp, and stage 2 for each head, given those three temperatures plus ambient, the head's design values, its four fitted numbers, and the compressor values. Check a hand-computed single head, the ideal-gas limit, the relief cap, and that two heads on one compressor share the swing.

- [ ] **2.3 Test-stand data and load-map fit** · Auto
  Readers for the cooldown traces (time, T1, TPR, T2) and load maps (Q1, Q2, T1, TPR, T2) in `data/sn09/pulse_tube/`, then a least-squares fit of a head's four numbers to a load map, taking zero clamp load and the measured clamp temperature at each point, written as a per-unit parameter file beside the data. A synthetic load map generated by the model should recover its own parameters. Flag a packing constant far from one.

- [ ] **2.4 Cooler link** · Auto
  `library/links/cooler.py`: a link over the stage-1, clamp, and stage-2 nodes plus a named ambient node, holding one compressor and one or more heads and returning the model's flows into each node.

- [ ] **2.5 Bare-head cooldown test** · Auto
  Fit SN09's two load maps, then simulate the bare head on the test stand as a three-node network and compare T1, TPR, and T2 against both measured cooldowns. Give the clamp node the thermometer block's mass so zero clamp load emerges instead of being imposed. Needs 5.1, 5.3, 5.4, and 5.5, plus material files for the regenerator matrix and the flange's aluminum (6061-T6 unless the inventory says otherwise). Deliver a reproducible report saying whether the four numbers moved between load maps and, for any miss, which term the ADR blames and what follows.

- [ ] **2.6 Cooler configuration** · Auto, after 2.5 passes
  Register `cooler` and load its nodes, head parameter files, design values, and compressor values from TOML. FP1's two heads go on one link when they share a compressor and on two otherwise.

## 3. Plots

- [ ] **3.1 Heat-flow plot** · Auto
  Companion to `plot_temperature`: flows by link and node, with the sign convention, units, and phase-switch markers.

- [ ] **3.2 Measured-temperature overlay** · Auto
  Measured thermometer traces drawn over their simulated nodes, once 4.1 settles the mapping.

- [ ] **3.3 Graph plot** · Auto
  Draw a Graph in `plot/graph.py`: nodes and links labelled by name and type, multi-node links included.

## 4. Validation

SN09 and FP1 share the six-plate architecture. SN09 has one PT1.1 pulse tube and FP1 has two PT1.5s on the same plates. The run procedures are cold_control's `cooldown` (SN09) and `fp_cooldown` (FP1) routines. SN09's cooldown is `data/sn09/cooldown_2026-06-04.dat`: 2026-06-04 00:01 to 2026-06-06 23:00 UTC at 60 s, channels PRP, RGP, CFP, STP, ICP, MXP, and DHX, with MXP reaching 10 K at 63.1 h. Node masses by material for both systems are in `data/<system>/bom_by_node.csv`.

- [ ] **4.1 Measurement comparison conventions** · Discuss
  For each system: which thermometer maps to which node, the cooldown's time origin, the comparison window, how to treat invalid readings, and how to sample and aggregate the log-temperature residual so the new and legacy models are scored identically. Needs 5.8.

- [ ] **4.2 Measurements and residual** · Auto
  Prepare each system's mapped, aligned measurements with units and provenance, and compute the least-squares log-temperature residual over the chosen thermometers and window, MXP excluded. Identical traces score zero.

- [ ] **4.3 Legacy baseline** · Auto
  Run the legacy model on each system, archive its version, inputs, and traces under `simulations/validation/`, and score it with the same residual.

- [ ] **4.4 Validation-system configurations** · Auto
  Each system's TOML and runner from its inventory, procedure, and cooler data. Port whatever materials the inventories need that the library lacks, with sources. Needs every link in section 1, 2.6, and 5.9.

- [ ] **4.5 Full cooldowns and comparison** · Auto
  Run each system through its cooldown, archive inputs, version, output, and solver outcome, then score against the measurements and the legacy baseline, with overlays. The report says whether each system meets the v1 target and where any excess residual comes from.

- [ ] **4.6 Resolve the outcome** · Discuss
  Accept when both systems qualify, or pick the next bounded investigation and repeat. Accepted runs become regression tests under `tests/regression/`.

## 5. Inputs

Things only a human can find. Check what's already in the repository and the related resources first, and record what's adequate by reference.

- [o] **5.1 Test-stand inventory** · Human
  Masses and materials on each stage of the pulse-tube test stand: flange, shield, mounting hardware, the cold-end heat exchangers inside the head, and the thermometer block at the clamp, with what each bolts to. So far: mounting hardware 110.6 g, of which 39.7 g is copper; radiation shield 1.191 kg of copper; flange 453.6 g of aluminum. Still missing: which stage each part sits on, the material of the other 70.9 g of hardware, and the cold-head internal masses.

- [ ] **5.2 FP1 pulse-tube records** · Human
  The test-stand cooldown and load map for each of FP1's two PT1.5 heads, under `data/fp1/pulse_tube/` like SN09's.

- [ ] **5.3 Head geometry and regenerator matrix** · Human
  For PT1.1 and PT1.5, per the ADR's "Data requirements": the gas volume of each gas-filled component and whose temperature it takes, the clamp position as a fraction of the regenerator length, the matrix masses and materials per section, and the pulse-tube wall dimensions and material. Expected from the heads' designers.

- [ ] **5.4 Compressor values** · Human
  Per compressor model: displacement rate, waveform factor, relief swing, mean pressure, and valve frequency, with sources, or a note that no specification exists. Also which compressor drives each head on the test stand and in each fridge, and whether FP1's two heads share one.

- [ ] **5.5 Thermometer calibrations** · Human
  In `data/sn09/pulse_tube/`, T1 and TPR both read through Z1518, TPR climbs to 350 K during cooldown, and at rest T1, TPR, and T2 read 276, 294, and 309 K. Which sensor and calibration each reads through, where TPR is mounted, and corrected traces.

- [ ] **5.6 He-4 adsorption and accommodation data** · Human
  Langmuir constants for He-4 on the switches' charcoal, and accommodation coefficients across temperature for He-4 and N2. Corruccini, Vacuum 7–8, 19 (1959), also reproduced in Ekin, *Experimental Techniques for Low-Temperature Measurements*. Add them to the gas files and drop the xfail from the He-4 tests.

- [ ] **5.7 Surface emissivities** · Human
  Emissivity over temperature for each radiating finish in SN09 and FP1, as surface files in `library/surfaces/data/`.

- [ ] **5.8 Validation measurements** · Human
  For each system: the cooldown record, the cold_control routine logs for that run, the `monitor_ovc` P1 diagnostics, the P2 and P3 pressure traces through gas pushing, and channel metadata. FP1's cooldown is the main gap.

- [ ] **5.9 FP1 geometry from CAD** · Human
  Into `data/fp1/geometry.csv`, one row per value with its unit and source part: shield and vacuum-can dimensions and panel thicknesses; plate spacings and the cross-section of each support type; switch core cross-sections and lengths along the heat path, gap, facing area, and sealed gas volume; internal gas volumes of the still, ICP heat exchangers, each DHX, and the mixing chamber; each pulse tube's clamp fraction. SN09's geometry follows from these and its BoM shield masses.

Already supplied, for the switch link and the system configurations: each system has eleven gas-gap switches, PRP–RGP 2 × N2, RGP–CFP 2 × N2, CFP–STP 3 × He-4, STP–ICP 2 × He-4, and ICP–MXP 2 × He-4. All use the same core (500-000182), outer core (500-000183), and shell (500-000185). Each gap is 0.5 mm wide with 40 cm^2 of facing area. The shell is 0.50 mm 304 SS sheet, 25.0 mm OD and 56.4 mm long, giving a 39 mm^2 conduction cross-section over a 70 mm developed length, and is its own conduction link. Switches are filled to 1 atm at room temperature and sealed with a crush tube. Each holds one 0.197 cm^3 charge of granulated coconut charcoal (part 200-000021), about 0.1 g at 0.5 g/cm^3.
