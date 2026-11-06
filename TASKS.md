# Tasks

Build the v1 simulator in `docs/DESIGN.md`: the six-plate systems SN09 and FP1 from room temperature to the start of condensing, matching their measured cooldowns at least as well as the legacy model does. "As well" means the least-squares residual in log temperature across the system thermometers other than MXP.

## Reading this list
Only the human edits DESIGN and ADRs; when a discussion settles something, propose the wording.

Read your own entry, not the whole list. Stay focused on the single task you've been assigned. Mark a task `[o]` while you have it and `[x]` when it lands. Refer to tasks by number; a new task takes the next free number in its section. If work turns up that isn't here, surface that as a recommended addition to the list. An entry names the piece and the context needed to start it; findings and progress belong in the code, data, or docs they concern.

Three kinds of entry:
1. **Build**: first, ensure alignment between the agent and the user. Discuss any outstanding points. Then build and hand back the work for review. Iterate as needed.
2. **Discuss**: the output of the task is a decision, more detailed follow-up tasks, an ADR, or other documentation. Ensure that the result is produced in collaboration with the user; a decision that the user cannot defend and take ownership over is a failure state.
3. **Human**: work that agents cannot do with what they have access to on this machine. Help the user identify where to track down the information. Remember that doing so takes time and effort; focus on finding the minimum set of necessary outside input, rather than comprehensive details that may not be needed.

A completed entry stays while its numbered section has unfinished work. When every task in a section is done and its decisions live in ADRs, DESIGN, or self-documenting code, delete the section. Before deleting, update any remaining task that refers to it.

## 2. Cooler

- [x] **2.1 Cooler ADR** · Discuss
  Settle `docs/adrs/03-coupled-capacity-maps.md` with the human, superseding ADR 01.

- [x] **2.2 Cooler tasks** · Discuss, after 2.1
  Break the settled ADR into section 2's tasks with the human, using the `discretize` skill.

Work the rest of the section in list order. The first slice runs from 2.5 through 2.3: the canonical head fitted to the test-stand load maps.

- [x] **2.5 Head design data** · Build
  The head-model schema, loader, and TOMLs for what ADR 03 reads per design: regenerator section lengths and walls, pulse-tube walls, reservoir and dead volumes, and κ, λ, and σ from Snodgrass 2024. PT1.1 comes from `data/heads/pt1_1/geometry.csv` and PT1.5 from `data/heads/pt1_5/geometry.csv`, both with PT1.1's matrix filling.
  Context: ADR 03 "Parameters"; DESIGN "Cooler configuration"; `library/heads/`, `library/links/cooler/load.py`, `data/scripts/derive_head_model.py`, `docs/papers/pulse_tube/`.
  Verified by: both models load, with each value traced to its source.

- [x] **2.6 Manifold pressures** · Build, after 2.5
  ADR 03's gas inventory and drive: each head's mass swing, the charge and mean pressure, and the supply equation with a smooth relief, with the helium properties these need.
  Context: ADR 03 "Governing equations"; `library/links/cooler/drive.py`, `head.py`; `library/gases/helium.py`; `data/compressors.md`.
  Verified by: mean pressure and swing vary smoothly with temperature and follow the readings in `data/compressors.md`.

- [x] **2.7 Cooling power** · Build, after 2.6
  The head's heat flows of ADR 03, without losses: acoustic power, regenerator energy flow in both sections, pulse-tube wall conduction, and Q1, Q2, and QR. The fitted parameters become a_i and R_i.
  Context: ADR 03 "Governing equations" and "Regenerator network"; `library/links/cooler/`.
  Verified by: heat flows with continuous derivatives that reach a cold steady state at plausible parameters.

- [x] **2.8 Breakthrough loss** · Build, after 2.7
  Each pulse tube's breakthrough loss, with helium enthalpy and a smooth positive part.
  Context: ADR 03 "Governing equations"; `docs/papers/pulse_tube/Snodgrass2024.pdf`.
  Verified by: no loss while the stroke stays within the buffer, rising smoothly past it.

- [x] **2.12 Smooth measured specific heats** · Build
  Tabulated properties interpolate with continuous derivatives, such as a shape-preserving cubic, including where a measured table joins its lattice continuation, so the cooler's heat flows meet ADR 03 requirement 7 with the library's regenerator materials. Gd2O2S's peak at 5.24 K is the sharpest case.
  Context: `library/property_fits.py` (the `table` and `measured_with_lattice` forms); `library/materials/data/`, `library/gases/data/`.
  Verified by: unchanged values at the table points and continuous first derivatives, in the fit forms and the library's regenerator materials.

- [x] **2.9 Cooler configuration** · Build, after 2.7
  The cooler as a configured link: `"cooler"` in `LINK_TYPES`, with its drive and per-head parameter files read from TOML. Propose matching DESIGN "Cooler configuration" text.
  Context: DESIGN "Cooler configuration" and "Public API"; `library/links/load.py`, `library/links/cooler/load.py`.
  Verified by: a fridge TOML with a cooler builds and runs.

- [x] **2.10 Test-stand model** · Build, after 2.9
  The test stand as a fridge template, `data/test_stand/fridge.toml`: one head on one or two compressors, with the stage masses and radiation shield of `data/test_stand.md`; foil on the tubes is omitted for lack of geometry. An engine steady-state solver, with heater powers as its input, solves it at a load map's points.
  Context: `data/test_stand.md`, `data/compressors.md`; `simulations/validation/sn09_bare_head/` for the earlier version.
  Verified by: the solver's analytic steady states, and the template building.

- [x] **2.11 Load-map fit set** · Build
  The canonical fit's inputs from `data/test_stand/load_map_index_2026-11-05.csv`: one Denver map per PT1.5 head with raw data on file, in the standard setup, which the index notes when it uses two compressors. Among a head's repeats, the more standard run by the notes, otherwise the more recent. Each point carries its run's compressor count, valve frequency, and MLI state. The selection replaces the record list in `data/systems.md`.
  Context: the index; `data/test_stand.md`; `simulations/exploration/loadmaps/`; the archive at `../../data/test_stand/`.
  Verified by: each selected map traces to its index row and raw file.

- [ ] **2.3 Canonical pulse-tube head** · Build, after 2.10 and 2.11
  SN09's and FP1's head model: one PT1.5 head's a_i and R_i fitted to 2.11's set, replacing the S0009-fitted parameters. Ends the first slice.
  Context: ADR 03 "Requirements" and "Parameters"; `simulations/exploration/loadmaps/characterization.md`.
  Verified by: each run's residual, and Q1's coupling to T2 against ADR 03 requirement 4.

- [ ] **2.4 Test the canonical pulse-tube head** · Build, after 2.3
  Use the test stand cooldown data to validate the pulse tube model in a controlled setting.
  Context: `data/test_stand.md` "Use"; 2.10's model.

## 3. Coupled head and canonical fit

This section builds 2.3 and 2.4. ADR 03 as written cannot meet requirement 4, so the head gains the cold heat exchangers' gas voids and a stage-1 coupling term, `c E_2 x (1 - x)` with `x = beta T_2 / (1 + epsilon)` at the stage-2 boundary. 3.1 records the evidence. Out of scope: PT1.1's values (4.4), and `simulations/validation/test_stand_steady_state/` with 2.10's open sub-item.

- [x] **3.1 Findings and ADR 03 proposal** · Discuss
  `simulations/validation/test_stand_head/findings.md`: the energy-budget argument, the void and parameter scans, the coupling gut check, and proposed ADR 03 wording for the voids, the coupling term, and requirement 4. The human edits ADR 03 from it.
  Context: ADR 03; `simulations/exploration/loadmaps/summary.txt`; `data/heads/pt1_5/geometry.csv`.
  Verified by: the human's ADR 03 edit.

- [x] **3.2 Specific-heat integral** · Build
  Each material's specific-heat integral over temperature, alongside `conductivity_integral`, for every fit form.
  Context: `library/materials/load.py`, `library/property_fits.py`.
  Verified by: agreement with quadrature for each fit form.

- [x] **3.3 Matrix heat capacity at node boundaries** · Build, after 3.2
  In the regenerator energy flow, the matrix heat capacity at a node boundary is the mean specific heat between the neighbouring temperatures, `(H(Ta) - H(Tb)) / (Ta - Tb)`.
  Context: `library/links/cooler/head.py` (`_energy_flow`); `simulations/validation/test_stand_steady_state/report.txt`.
  Verified by: unchanged results for a constant specific heat, and a stable cold root where `report.txt` found an unstable one.

- [x] **3.4 Heat-exchanger voids** · Build, after 3.1
  Each stage's cold heat-exchanger void, read from `data/heads/pt1_5/geometry.csv` through `data/scripts/derive_head_model.py` into `library/heads/data/pt1_5.toml`. The void's gas is isothermal at its stage temperature. It enters the gas inventory, the mass swing, the boundary mass flows, and the stage's pore heat capacity. A design file without voids loads them as zero.
  Context: ADR 03 as edited in 3.1; `library/heads/`, `library/links/cooler/head.py`, `drive.py`.
  Verified by: zero voids reproduce today's heat flows, and a hand-computed void mass swing.

- [x] **3.5 Stage-coupling term** · Build, after 3.1
  The design's coupling coefficient `c` in `HeadModel` and its TOML, zero when absent, and `c E_2 x (1 - x)` in stage 1's heat flow.
  Context: ADR 03 as edited in 3.1; `library/heads/`, `library/links/cooler/head.py`.
  Verified by: `c = 0` reproduces today's heat flows, and a hand-computed term.

- [x] **3.6 Stage-2 gate** · Discuss, after 3.3 and 3.4
  With the voids and 3.3, at T1 = 50 K and T2 = 3.4 K over a2 and R2: can dQ2/dT2 reach -0.4 to -0.7 W/K, with heat into the node positive, on a stable state? If not, the result goes to the human before 3.8.
  Context: 3.1's findings; `simulations/validation/test_stand_steady_state/coupling.py` for the derivative method.
  Also: `simulations/validation/test_stand_head/boundary_capacity_stability.txt`. After 3.3 the template's no-load cold root is stable, but stage 2 sits at the 2.1 K helium floor, some Gd2O2S nodes fall below it, and their temperatures alternate node to node. The point-sampled model before 3.3 also relaxed to a stable root, at T2 = 2.32 K.
  Result, `simulations/validation/test_stand_head/stage2_gate.txt`: yes, on the cold matrix branch, which is stable with the stages held and with heater power held. dQ2/dT2 reaches -0.70 W/K at a2 = 0.038, R2 = 4.3e10 and -0.63 at a2 = 0.034, R2 = 4.8e10, with a stage-2 heater of 0.57 to 0.67 W; a1, R1, and `c` move it by under 0.02 W/K. The branch ends where the warmest Gd2O2S node reaches about 6.2 K, beyond which the matrix settles warm and stage 2 cannot hold 3.4 K; along it, -0.55 W/K at a2 = 0.03 was the closest to the measured -0.46.

- [x] **3.7 Robust steady state** · Build
  `solve_steady_state` relaxes in log temperature and finishes with Newton. Its temperature floor is an argument defaulting to 2.1 K, the helium-property limit, and a state pinned there raises promptly. Callers continue between points by passing the previous solution.
  Context: `engine/steady_state.py`, `tests/engine/test_steady_state.py`.
  Verified by: the analytic steady states, a prompt failure where no state exists above the floor, and the template converging at a load-map point.

- [ ] **3.8 Canonical PT1.5 fit** · Build, after 3.5, 3.6, and 3.7
  a1, a2, R1, R2, and `c` fitted jointly to every map in `data/test_stand/fit_set.csv`, with equal weight per map and each run at its own compressor count and valve frequency. At each point the stages are held at the measured T1 and T2, the matrix is solved, and the stage heat flows are compared with the applied heater powers. Every run uses the template's bare-shield radiation link. a, R, and the frequency go to `data/heads/pt1_5/parameters.toml`, `c` to `pt1_5.toml`, and the test-stand template points at both. `data/scripts/derive_head_model.py` writes `c` and its source into `pt1_5.toml`.
  Context: `simulations/validation/test_stand_head/`; `simulations/exploration/loadmaps/characterization.md`.
  Verified by: each run's residual, with S0028 reported separately, and the four couplings against figure 07's grid cells.

- [ ] **3.9 Stand cooldown set** · Build
  `read_cooldown` accepts ISO timestamps, as in S0007's April log. `data/test_stand/cooldown_set.csv` lists the stand cooldowns, with trimmed verbatim copies under `data/test_stand/<serial>/`: S0007, S0009, S0019, S0035, and S0038 on one compressor, and S0028 on two. The experiment reads the list with `csv`; there is no library reader.
  Context: `measurements/`; `data/test_stand.md`; the archive at `../../data/test_stand/`.
  Verified by: an ISO-timestamp test, and each copy traced to its archive log.

- [ ] **3.10 Stand cooldowns** · Build, after 3.8 and 3.9
  The fitted template's cooldown on one and on two compressors, scored on T1 below 219 K and on T2 throughout. If the transient solver fails, the remedy goes to the human.
  Context: `data/test_stand.md` "Use"; `simulations/validation/test_stand_head/`.
  Verified by: each simulated cooldown's score against the six measured cooldowns.

- [ ] **3.11 Fit and cooldown figures** · Build, after 3.8 and 3.10
  PNG and PDF in `simulations/validation/test_stand_head/`: the model's standard points against every PT1.5 run's spread on one and two compressors; the four couplings over figure 07's cells; each run's heater-power residuals; and T1 and T2 cooldowns over the six stand cooldowns.
  Context: `simulations/exploration/loadmaps/figures/02_standard_points.png` and `07_capacity_partials.png`.

## 4. Validation

Each system's records, procedure, and inventory are in `data/systems.md`.

- [x] **4.1 Measurement comparison conventions** · Discuss
  Times are in `data/systems.md`.
  - Mapping: each plate thermometer maps to the node of the same name, and FP1's top DHX maps to DHX. Scored: SN09's PRP, RGP, CFP, STP, ICP, and DHX; FP1's PRP, RGP, CFP, STP, and top DHX. SN09's MXP is overlaid, not scored. FP1's MXP channel, a DHX sensor reading up to 25 K below top DHX, and its CMN are unused.
  - Time origin: compressor start. Each node starts at its measured temperature there.
  - Window: from compressor start to the last sample before condensing begins. Validation runs end at the window's length (`max_time`), not at a condensing condition.
  - Invalid readings: NaN, zero or negative, or within a run of 10 or more identical consecutive samples. They are excluded from the residual and the overlays.

- [x] **4.7 Residual scoring** · Discuss, after 4.1
  How to sample and aggregate the log-temperature residual over 4.1's window and thermometers, so the new and legacy models are scored identically.
  - Samples: every valid measured sample of each scored thermometer in the window. The simulated node temperature is interpolated there, linearly in time in log temperature.
  - Score: one per system, the RMS of `ln(T_simulated / T_measured)` pooled over all of the system's samples. Each thermometer's RMS is reported alongside.
  - Coverage: a simulation that ends before the window's last sample is not scored.

- [x] **4.2 Measurements and residual** · Build, after 4.7
  Each system's mapped, aligned measurements, and the log-temperature residual over the chosen thermometers.

- [x] **4.3 Legacy baseline** · Build, after 4.2
  The legacy model run on each system and scored with the same residual, archived under `simulations/validation/`. Each system's baseline is its best legacy score over committed legacy inputs: any committed legacy version, cooler data, and mix of committed saved-system pieces, with no parameter fit to the measured cooldown and no change to legacy code. The runner runs legacy at the commit in its inputs; moving to a later legacy commit is recorded here.
  Runs `sn09_tuned`, `sn09_corrected`, `fp1`, and `fp1_generous` are in `simulations/validation/legacy_baseline/`; the baselines are `sn09_tuned` and `fp1_generous`.
  - [x] Plot FP1's archived baseline with all stages on one panel: `simulations/validation/legacy_baseline/fp1_baseline_all_stages.png` and `.pdf`.
  - [x] Plot FP1's measured cooldown stages on one panel, matching the baseline figure: `simulations/validation/legacy_baseline/fp1_data_all_stages.png` and `.pdf`.

- [ ] **4.4 Validation-system configurations** · Build, after section 2
  Each system's fridge template, `data/<system>/fridge.toml`, from its inventory, geometry, `data/systems.md`, procedure, and cooler data, porting any materials the library lacks. Its validation run directory, such as `simulations/validation/sn09_cooldown/`, holds a copy of the template as its config and a runner that sets initial temperatures and `max_time` from `simulations/validation/measurements.py`.
  PT1.1 needs its own coupling coefficient `c` and cold heat-exchanger void volumes. The head-model loader treats them as zero until `library/heads/data/pt1_1.toml` states them.

- [ ] **4.5 Full cooldowns and comparison** · Build, after 4.3 and 4.4
  Each system run through its cooldown and scored against the measurements and the legacy baseline, with overlays.

## 7. Repository health (once all other tasks completed)

- [ ] **7.1 Test suite** · Build
  Make sure tests are necessary, relevant, and catch problems that are *likely* and *impactful*. A broad test suite of low-leverage tests is a failure state. We also need to make sure our tests do not take too long to run, because agents are instructed to run them on every single /implement call.

- [ ] **7.2 DESIGN reconciliation** · Discuss
  Reconcile DESIGN.md with code. This is our primary design document, and should serve as the equivalent of an engineering drawing for our codebase. It needs to remain concise, human-readable, and up-to-date with the repo. The "Cooler configuration" example predates the regenerated `pt1_1.toml`.

- [ ] **7.3 ADR reconciliation** · Discuss
  Reconcile ADRs with specific implementations. Clean up prose and logical order so that a reader can easily follow why the decision was made.

- [ ] **7.4 Simulations directory** · Build
  Clean up the entire simulations directory. It is a nest of terrible code and unreviewed decisions. Keep what's important, remove the obvious slop.

- [ ] **7.5 Comments and documentation** · Build
  De-slop the comments and documentation across the codebase. Remove "Memento" style comments that record conversations, old decisions, information recorded elsewhere, etc.

- [ ] **7.6 Module interfaces** · Discuss
  Audit the interfaces and seams between modules. This is cheap to change now, before a v1.0.0 release.

- [ ] **7.7 Findings and model upgrades** · Discuss
  Extract scientific findings and major model upgrades compared to legacy code. Split into 3 categories: (1) new scientific findings, (2) model architecture changes that directly supported those findings, and (3) model architecture changes that will support future work.

- [ ] **7.8 Agent navigation** · Discuss, after 7.2 and 7.3
  Improve agent navigation and interaction economy. Make it fast and cheap for agents to edit the repository. AGENTS rule 1 and `CONTEXT.md` are in place. Remaining: measure the context a fresh agent reads before it can start implementing (plotting, a heat-transfer link, cooler work), split into preloaded and repository reading; target under 30k tokens.

- [x] **7.9 Data directory cleanup** · Build
  Keep source measurements, reference data, and system inputs in `data/`; put processing tools and exploratory work in their owning directories. Update path consumers and document the layout.

- [x] **7.10 Data preparation tools location** · Build
  Keep data preparation tools in `data/scripts/`, alongside the measurements and reference inputs they prepare. Update path consumers and documentation.
