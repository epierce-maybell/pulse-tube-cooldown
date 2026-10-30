# Tasks

Build the v1 simulator in `docs/DESIGN.md`: the six-plate systems SN09 and FP1 from room temperature to the start of condensing, matching their measured cooldowns at least as well as the legacy model does. "As well" means the least-squares residual in log temperature across the system thermometers other than MXP.

## Reading this list

DESIGN.md is the specification and `docs/adrs/` holds the decisions behind it. AGENTS.md says how we work. The legacy code at `../pulse_tube_cooldown_simulation`, branch `Advanced_Feature_SIM`, is the reference for physics and data. The manufacturer's record of every shipped head, with its tuning, valve frequency, Münster load map, and test configuration, is `docs/pt_performance_overview_2026-10-30.pdf`; heads are named by serial there (S0009 is SN09's head, S0019 is `data/sn019`). Only the human edits DESIGN; when a discussion settles something, propose the wording.

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

- [x] **1.4 Multi-layer insulation** · Discuss
  How MLI modifies radiative transfer, what it's parameterized by, and whether it's a separate link or an option on radiation. Legacy has MLI functions to start from. Then build it. Built as `multi_layer_insulation`: the Lockheed correlation with the blanket's coefficients and source in the link's table. DESIGN wording pending.

- [ ] **1.5 Joint contact** · Discuss
  Contact conductance at bolted joints and the parameters or measurements it needs. Legacy has `thermal_conductance_calc.py` and `conductance_systems.json`. Then build it.

## 2. Cooler

The model is `docs/adrs/coupled-capacity-maps.md`. The milestone is reproducing the bare head's cooldown on the test stand. Registering the cooler as a configurable link waits on that.

- [x] **2.1 Helium properties** · Auto
  Add `he34pak` (`../../dev/he34pak`; pure He-4 is its `hepak` subpackage) as a dependency and expose He-4 density, isobaric specific heat, and isobaric expansion coefficient as array functions of temperature and pressure over 2 K to 300 K at the head's pressures. Check the near-ideal limit at 300 K and the real-gas regime at 4 K.

- [ ] **2.2 Section physics** · Auto
  The ADR's network equations for one head at given supply and return pressures, as plain functions of arrays and floats with no new types. The mass flow across each column boundary from the gas volumes beyond it, with each matrix node's pore gas at its own temperature and each pulse tube along a linear profile from the warm flange to its stage; the boundary resistances; and the stage enthalpy sinks from the swing and the waveform factor κ, with the real-gas correction at stage 2. Heat flows into the warm flange, column nodes, and stages sum to zero. The tube walls are separate conduction links. Check a hand-computed head, the ideal-gas limit, that one node per section reproduces the ADR's section equations, and that splitting a layer into more nodes leaves the steady state unchanged. Report the largest utilization per section as a diagnostic rather than capping it.

- [ ] **2.8 Drive solve** · Auto
  The mean pressure from the gas charge over every volume at its temperature, the manifold's warm volume included, then the pressure swing from the drive constraint for any number of heads and compressors on one manifold, under the smallest relief cap. Plain functions over 2.2's mass flows. Check the relief cap, that two heads on one compressor share the swing, and that two compressors on one head stay at the relief cap to a colder temperature than one.

- [ ] **2.9 Head description** · Discuss
  How a configuration describes a head: its two columns and their layers (matrix material, matrix mass, pore gas volume), its pulse tubes, and the clamp as the stage-2 column node whose span contains `x_R`. How SN09's `data/sn09/pulse_tube/geometry.csv` reaches a configuration with one node per material layer, and whether any types are needed. This settles how a configuration spells a head's columns for 2.4 and 2.6. Ends in a DESIGN update.

- [ ] **2.3 Test-stand data and load-map fit** · Auto
  Readers for the cooldown traces (time, T1, TPR, T2) and load maps (Q1, Q2, T1, TPR, T2) in `data/sn09/pulse_tube/`, then a least-squares fit of a head's four numbers to a load map, including its tube walls' conduction, written as a per-unit parameter file beside the data. The residuals are the stage flows against the heater powers at the measured T1 and T2, with each column solved to steady state and zero clamp load at every point. TPR is not a fit input; report the predicted clamp temperature beside it as a diagnostic. A synthetic load map generated by the model should recover its own parameters.
  Fit SN09 on `load_map_2026-05-14.dat` alone, as one F-70 at 1.0 Hz: Münster recorded S0009 that way within a kelvin of the May points, and the August map is a prediction target in 2.5, not a fit input. The valve frequency is a per-unit number taken from the overview's shipped setting; the `f:` field in the load-map headers is typed in by the operator (S0007's header reproduces its shipped settings, S0009's and S0034's do not match theirs). Report the packing constants rather than flagging them: shipped double-inlet openings range from 0 to 3 turns on DIV1B and 3 to 9.75 on DIV2B, and the bypassed regenerator flow lands in the packing constants, so they will differ between heads.

- [ ] **2.4 Cooler link** · Auto
  `library/links/cooler.py`: one link per manifold over the warm flange, every head's column nodes, and its stage nodes, holding the compressors and heads on that manifold and returning the model's flows into each node. The Nodes for a column and the link's description of it come from the one description settled in 2.9, so matrix mass has one source.

- [ ] **2.5 Bare-head cooldown test** · Auto
  With the May fit from 2.3, simulate the bare head on the test stand as a network of the warm flange, the two stages, and the column nodes, with its tube walls as conduction links and stage-1 radiation from the 300 K chamber through the stand's MLI (5.1), and compare T1 and T2 against `cooldown_2026-05-14.dat` over the ranges in `data/test_stand.md`. The column nodes carry only their matrix (5.3), so zero clamp load emerges instead of being imposed. Plot TPR against the clamp node without scoring it.
  Then, with the same four numbers and two compressors, predict `load_map_2026-08-25.dat` and `cooldown_2026-08-24.dat`. The August run's compressor count is unrecorded; its lower T1 and higher zero-load T2 against May match S0028's change from one compressor to two in the overview, so two is the stated assumption. A matching prediction confirms the assumption and the ADR's drive independence together; a miss says which one failed. Needs 1.4, 2.7, 5.1, 5.3, 5.4, and 5.5, plus a material file for the flange's aluminum (6061-T6 unless the inventory says otherwise). 

- [x] **2.7 Regenerator matrix materials** · Auto
  Material files in `library/materials` for lead, Er3Ni, HoCu2, and Gd2O2S: specific heat from 2 K to 300 K and density, with sources. The rare-earth heat capacities have peaks below 10 K that the fits must keep. Geometry and layer fractions come from 5.3.
  Done: lead from Arblaster 2012; the rare earths as `measured_with_lattice` fits. Estimates: the Gd2O2S Debye temperature (325 K), the Er3Ni density, and Er3Ni above 30 K, which recovers 75% of its magnetic entropy.

- [ ] **2.6 Cooler configuration** · Auto, after 2.5 passes
  Register `cooler` and load its nodes, head parameter files, design values, and compressor values from TOML. Spell a head's columns as 2.9 settles, so that their matrix Nodes and the link are built from one description. Each manifold is one link, so FP1's two heads share a link when they share a manifold.

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

- [x] **5.1 Test-stand inventory** · Human
  Masses and materials on each stage of the pulse-tube test stand. Stage 1: the head's upper cold-end heat exchanger, 392.8 g of copper (43.82 cm^3); the flange, 453.6 g of aluminum; the radiation shield, 1.191 kg of copper; and 70.9 g of mounting hardware, taken as stainless fasteners. Stage 2: the head's lower cold-end heat exchanger, 421.8 g of 304 stainless (52.55 cm^3, material as assigned in CAD 350-000162); and 39.7 g of copper mounting hardware. Clamp: a Lake Shore Cernox only, of negligible mass.
  MLI: the overview's test configurations wrap the two regenerators, pulse tube 1, and both stages in 4 to 12 layers each, varying per run; S0009's Münster run had 7 on each regenerator and pulse tube 1, 4 on stage 1, and 5 on stage 2. Assume those counts for the Denver runs. Layer density is unmeasured; assume 2000 layers/m (20 layers/cm), the low end of the 20 to 30 layers/cm typical of installed blankets. Across 10 to 30 layers/cm, the flux through 4 layers from 300 K to 45 K changes by about a factor of 3.5.

- [ ] **5.2 FP1 pulse-tube records** · Human
  The test-stand cooldown and load map for each of FP1's two PT1.5 heads, under `data/fp1/pulse_tube/` like SN09's. Each head's Münster key points, valve settings, and shipped frequency are in the overview by serial, so the first step is the two serials.

- [x] **5.3 SN09 head geometry and regenerator matrix** · Human
  For the SN09 head, per the ADR's "Data requirements": the gas volume of each gas-filled component and whose temperature it takes, the clamp position as a fraction of the regenerator length, the matrix materials and layer fractions per section, and the tube dimensions and material. Recorded in `data/sn09/pulse_tube/geometry.csv`: tube radii, lengths, and the clamp position from CAD; matrix layers from the TransMIT PTD 420RX paper in `docs/papers/Falter2024.pdf`, the same head family; wall thickness, porosity, and the gas volumes derived from them marked as estimates. FP1's heads are in 5.9.

- [o] **5.4 Compressor values** · Human
  Every compressor is a Sumitomo F-70H. The valve frequency is set per head at tuning, 0.9 to 1.1 Hz in the overview, S0009 at 1.0 Hz, and is a per-unit number rather than a drive value; the fridges run the shipped setting unless retuned. Needed: displacement rate, waveform factor, relief swing, gas charge, and manifold volume, with sources, or a note that no specification exists.
  So far: the test stand runs one compressor on the head, SN09 runs two compressors on its one head, and FP1 is believed to run each of its two heads on its own compressor. Every head in the overview is filled to 17.0 to 17.1 bar at rest. The running mean pressure is not constant: S0023's Münster cooldown trace (overview, p. 33) shows the return pressure rising from about 6.8 bar to 8.7 bar an hour in and settling near 6.5 bar at base, and a base swing at that return would exceed the relief cap unless the mean has fallen as gas moved into the cold volumes. The model takes the mean from the charge and the gas inventory over every volume at its temperature, so the manifold's warm volume is a drive value. On the stand it is two 20 m DN25 flexlines and two 5 L buffers; that is the default for the fridges until measured. The F-70H/F-70L/F-70LP Operating Manual (SHI 267472A, rev. E, 2011) gives, for a 4K pulse-tube head: equalization pressure at 20 C of 13.6, 14.3, or 16.6 bar gauge nominal, with the table not saying which configuration takes which; a factory charge of 1430 to 1460 kPa at 20 C; operating supply pressure of 18 to 25 bar gauge, varying with head and load; shutdown on return pressure above 140 psig or below 15 psig; and power draw of 6.7 to 7.2 kW on 50 Hz mains and 8.0 to 8.5 kW on 60 Hz. Neither that manual nor the F-70H/F-70L Service Manual gives a displacement rate or internal-relief setpoint. The waveform factor appears only squared and multiplied into each valve conductance, so it is fixed at 1 and absorbed into the conductances.
  Front panel of one of FP1's F-70H compressors during cooldown, with the head at 260 to 280 K, where the swing sits at its relief cap: supply 360 to 380 psig on the analog gauge and return 90 psig on the LCD. That puts the relief swing at 270 to 290 psi (18.6 to 20.0 bar) and the running mean pressure at 230 psig (15.9 bar gauge, 16.9 bar absolute), closest to the manual's 16.6 bar gauge equalization option. The F-70's scroll capsule is a Copeland ZCH72C3G/C4G, rated 3500 rpm on 60 Hz; Copeland does not publish its displacement, but the ZR72KC on what appears to be the same frame sweeps 5.98 in^3/rev, or 727 ft^3/h (20.6 m^3/h) at 3500 rpm. That swept rate is an upper bound on the displacement rate, since helium leakage, warm suction gas, and injected oil all reduce delivery. Until a reading at base replaces it, the displacement rate is the swept rate, 5.72e-3 m^3/s, which delivers 6.6 g/s at the 90 psig return and 300 K. Still missing: supply and return pressures with a system at base, with the stage temperatures at that moment.

- [x] **5.5 Thermometer calibrations** · Human
  Which test-stand readings hold, recorded in `data/test_stand.md` from twelve cooldowns of at least six heads in `../../data/test_stand/`. T1 and T2 hold apart from fixed dropout bands and T1's latched warm readings. TPR is repeatable but is not a clamp measurement: during cooldown it rises in a pattern that depends on the head, reaching 351 K on SN09, and at base it sits 8% to 38% of the way down the T1-to-T2 drop on a clamp 66% of the way along the column, while SN09's own RGP sits 88% of the way down at base despite carrying loads the stand's clamp lacks. Every head reads through the same curves, so this is not a calibration error. Where TPR is mounted, and whether its sensor moves between heads, is unknown.

- [ ] **5.6 He-4 adsorption and accommodation data** · Human
  Langmuir constants for He-4 on the switches' charcoal, and accommodation coefficients across temperature for He-4 and N2. Corruccini, Vacuum 7–8, 19 (1959), also reproduced in Ekin, *Experimental Techniques for Low-Temperature Measurements*. Add them to the gas files and drop the xfail from the He-4 tests.

- [ ] **5.7 Surface emissivities** · Human
  Emissivity over temperature for each radiating finish in SN09 and FP1, as surface files in `library/surfaces/data/`.

- [ ] **5.8 Validation measurements** · Human
  For each system: the cooldown record, the cold_control routine logs for that run, the `monitor_ovc` P1 diagnostics, the P2 and P3 pressure traces through gas pushing, and channel metadata. FP1's cooldown is the main gap.

- [ ] **5.9 FP1 geometry from CAD** · Human
  Into `data/fp1/geometry.csv`, one row per value with its unit and source part: shield and vacuum-can dimensions and panel thicknesses; plate spacings and the cross-section of each support type; switch core cross-sections and lengths along the heat path, gap, facing area, and sealed gas volume; internal gas volumes of the still, ICP heat exchangers, each DHX, and the mixing chamber; for each pulse tube, the clamp distances to PRP and CFP and the outer radius and length of its four tubes, as in `data/sn09/pulse_tube/geometry.csv`, once 2.5 has tested the SN09 values. SN09's geometry follows from these and its BoM shield masses.

Already supplied, for the switch link and the system configurations: each system has eleven gas-gap switches, PRP–RGP 2 × N2, RGP–CFP 2 × N2, CFP–STP 3 × He-4, STP–ICP 2 × He-4, and ICP–MXP 2 × He-4. All use the same core (500-000182), outer core (500-000183), and shell (500-000185). Each gap is 0.5 mm wide with 40 cm^2 of facing area. The shell is 0.50 mm 304 SS sheet, 25.0 mm OD and 56.4 mm long, giving a 39 mm^2 conduction cross-section over a 70 mm developed length, and is its own conduction link. Switches are filled to 1 atm at room temperature and sealed with a crush tube. Each holds one 0.197 cm^3 charge of granulated coconut charcoal (part 200-000021), about 0.1 g at 0.5 g/cm^3.
