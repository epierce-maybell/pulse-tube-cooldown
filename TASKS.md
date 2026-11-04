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

- [o] **2.1 Cooler ADR** · Discuss
  Settle `docs/adrs/03-coupled-capacity-maps.md` with the human, superseding ADR 01.

- [ ] **2.2 Cooler tasks** · Discuss, after 2.1
  Break the settled ADR into section 2's tasks with the human, using the `discretize` skill.

## 3. Plots

- [x] **3.1 Heat-flow plot** · Build

- [ ] **3.2 Measured-temperature overlay** · Build, after 4.1
  Measured thermometer traces drawn over their simulated nodes.

- [x] **3.3 Graph plot** · Build

## 4. Validation

Each system's records, procedure, and inventory are in `data/systems.md`.

- [ ] **4.1 Measurement comparison conventions** · Discuss
  For each system: which thermometer maps to which node, the cooldown's time origin, the comparison window, how to treat invalid readings, and how to sample and aggregate the log-temperature residual so the new and legacy models are scored identically.

- [ ] **4.2 Measurements and residual** · Build, after 4.1
  Each system's mapped, aligned measurements, and the log-temperature residual over the chosen thermometers.

- [ ] **4.3 Legacy baseline** · Build, after 4.2
  The legacy model run on each system and scored with the same residual, archived under `simulations/validation/`.

- [ ] **4.4 Validation-system configurations** · Build, after section 2
  Each system's TOML and runner from its inventory, procedure, and cooler data, porting any materials the library lacks.

- [ ] **4.5 Full cooldowns and comparison** · Build, after 4.3 and 4.4
  Each system run through its cooldown and scored against the measurements and the legacy baseline, with overlays.

- [ ] **4.6 Resolve the outcome** · Discuss, after 4.5
  Accept when both systems qualify, or pick the next bounded investigation. Accepted runs become regression tests under `tests/regression/`. The aluminum surface finish is assumed, so emissivity scaling is one candidate.

## 5. Inputs

- [x] **5.1 Test-stand inventory** · Human

- [ ] **5.2 FP1 pulse-tube records** · Human
  The test-stand cooldown and load map for each of FP1's two PT1.5 heads, under `data/fp1/pulse_tube/` like SN09's. The overview indexes heads by serial, so the first step is the two serials.

- [x] **5.3 SN09 head geometry and regenerator matrix** · Human

- [o] **5.4 Compressor values** · Human
  The drive values ADR 03 needs for each manifold. What is known is in `data/compressors.md`. Still missing: supply and return pressures with a system at base, and the stage temperatures at that moment.

- [x] **5.5 Thermometer calibrations** · Human

- [x] **5.6 He-4 adsorption and accommodation data** · Human

- [x] **5.7 Surface emissivities** · Human

- [x] **5.8 Validation measurements** · Human

- [x] **5.9 FP1 geometry from CAD** · Human

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
