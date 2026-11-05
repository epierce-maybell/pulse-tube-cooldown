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

- [ ] **2.5 Head design data** · Build
  The head-model schema, loader, and TOMLs for what ADR 03 reads per design: regenerator section lengths and walls, pulse-tube walls, reservoir and dead volumes, and κ, λ, and σ from Snodgrass 2024. PT1.1 comes from `data/sn09/pulse_tube/geometry.csv` and PT1.5 from `data/fp1/geometry.csv`, both with PT1.1's matrix filling.
  Context: ADR 03 "Parameters"; DESIGN "Cooler configuration"; `library/heads/`, `library/links/cooler/load.py`, `scripts/data/derive_head_model.py`, `docs/papers/pulse_tube/`.
  Verified by: both models load, with each value traced to its source.

- [ ] **2.6 Manifold pressures** · Build, after 2.5
  ADR 03's gas inventory and drive: each head's mass swing, the charge and mean pressure, and the supply equation with a smooth relief, with the helium properties these need.
  Context: ADR 03 "Governing equations"; `library/links/cooler/drive.py`, `head.py`; `library/gases/helium.py`; `data/compressors.md`.
  Verified by: mean pressure and swing vary smoothly with temperature and follow the readings in `data/compressors.md`.

- [ ] **2.7 Cooling power** · Build, after 2.6
  The head's heat flows of ADR 03, without losses: acoustic power, regenerator energy flow in both sections, pulse-tube wall conduction, and Q1, Q2, and QR. The fitted parameters become a_i and R_i.
  Context: ADR 03 "Governing equations" and "Regenerator network"; `library/links/cooler/`.
  Verified by: heat flows with continuous derivatives that reach a cold steady state at plausible parameters.

- [ ] **2.8 Breakthrough loss** · Build, after 2.7
  Each pulse tube's breakthrough loss, with helium enthalpy and a smooth positive part.
  Context: ADR 03 "Governing equations"; `docs/papers/pulse_tube/Snodgrass2024.pdf`.
  Verified by: no loss while the stroke stays within the buffer, rising smoothly past it.

- [ ] **2.9 Cooler configuration** · Build, after 2.7
  The cooler as a configured link: `"cooler"` in `LINK_TYPES`, with its drive and per-head parameter files read from TOML. Propose matching DESIGN "Cooler configuration" text.
  Context: DESIGN "Cooler configuration" and "Public API"; `library/links/load.py`, `library/links/cooler/load.py`.
  Verified by: a fridge TOML with a cooler builds and runs.

- [ ] **2.10 Test-stand model** · Build, after 2.9
  The test stand as a fridge: one head on one or two compressors, with the stage masses and MLI of `data/test_stand.md`, solved to steady state at a load map's heater powers.
  Context: `data/test_stand.md`, `data/compressors.md`; `simulations/validation/sn09_bare_head/` for the earlier version.
  Verified by: a load-map surface and a cooldown at placeholder parameters.

- [ ] **2.11 Load-map fit set** · Build
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

## 3. Plots

- [x] **3.1 Heat-flow plot** · Build

- [ ] **3.2 Measured-temperature overlay** · Build, after 4.1
  Measured thermometer traces drawn over their simulated nodes.

- [x] **3.3 Graph plot** · Build

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

- [ ] **4.3 Legacy baseline** · Build, after 4.2
  The legacy model run on each system and scored with the same residual, archived under `simulations/validation/`.

- [ ] **4.4 Validation-system configurations** · Build, after section 2
  Each system's fridge template, `data/<system>/fridge.toml`, from its inventory, geometry, `data/systems.md`, procedure, and cooler data, porting any materials the library lacks. Its validation run directory, such as `simulations/validation/sn09_cooldown/`, holds a copy of the template as its config and a runner that sets initial temperatures and `max_time` from `simulations/validation/measurements.py`.

- [ ] **4.5 Full cooldowns and comparison** · Build, after 4.3 and 4.4
  Each system run through its cooldown and scored against the measurements and the legacy baseline, with overlays.

- [ ] **4.6 Resolve the outcome** · Discuss, after 4.5
  Accept when both systems qualify, or pick the next bounded investigation. Accepted runs become regression tests under `tests/regression/`. The aluminum surface finish is assumed, so emissivity scaling is one candidate.

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
