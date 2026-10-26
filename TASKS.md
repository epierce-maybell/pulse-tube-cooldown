# Tasks

The work from an empty package to v1, as specified in `docs/DESIGN.md`. v1 is a full cooldown simulation whose log-temperature residual is at or below the legacy model's.

## How to use this list

Tasks are grouped by area. Within an area, tasks run in order. Each area marks where the **first vertical slice** ends: the full canonical run (configure, simulate, save, CLI) on a toy configuration, checked against an analytic or hand-computed answer. Slice tasks come first.

The next task is the first unchecked task in its area whose earlier tasks are done and whose **Needs** are met. Every link task after Conduction needs task 6 (config loader), because each link adds its own registry entry and configuration schema.

There are three task types:

- **Auto**: fully scoped. An agent implements it unattended, and the work is reviewed afterwards.
- **Discuss**: starts with a working session with Ethan to settle the listed decisions. A short session goes on to implementation in the same session and updates `docs/DESIGN.md` afterwards. A longer one ends with the decisions written into `docs/DESIGN.md`. The entry then becomes **Auto**, points at that section, and a later agent implements it.
- **Human**: data or information only Ethan can supply.

Each entry lists the context it touches, what to do, and how the work is verified. Unit tests are part of every implementation task. `docs/DESIGN.md` is the authority. Where it and `docs/proposals/` disagree, `docs/DESIGN.md` holds. The brownfield at `../pulse_tube_cooldown_simulation` is a reference for physics and data, not for structure.

## Core

- [ ] **1. Data model and equations** · Auto
  - Context: DESIGN.md "Data model", "Phases and Procedure", "Engine"; `core/types.py`, `core/errors.py`, `engine/equations.py`.
  - Do: Define the data model types, including `Gas`, and the package's error types. Implement heat capacity, link and net heat flows, temperature rate of change with fixed boundary nodes, and the Jacobian.
  - Verified by: Hand-computed heat flows and rates on small graphs; the Jacobian agrees with finite differences.

- [ ] **2. Integration** · Auto
  - Context: DESIGN.md "Engine", "Phases and Procedure"; `engine/solver.py`.
  - Do: Implement `integrate_phase`, `run_procedure` with phase-ending and stop events, `reconstruct_heat_flow` and `build_result`, including duplicate entries at phase switches.
  - Verified by: Analytic multi-node relaxation; a two-phase run switches at the expected time and records both sides of the switch.

*End of first vertical slice.*

## Materials

- [ ] **3. Property fits and material loading** · Auto
  - Context: DESIGN.md "Dilution refrigerator library"; `library/property_fits.py`, `library/materials/`; brownfield `ptc_material_properties.json`.
  - Do: Implement the property-fit forms in `PROPERTY_FITS`, the material TOML format and `load_material`. Port OFHC copper and 304 stainless steel.
  - Verified by: Ported materials reproduce reference values at spot temperatures.

*End of first vertical slice.*

- [ ] **4. Remaining materials** · Auto
  - Context: `library/materials/data/`; brownfield `ptc_material_properties.json`.
  - Do: Port the remaining brownfield materials, recording each source.
  - Verified by: Every material loads and matches brownfield values at spot temperatures.

## Configuration and API

- [ ] **5. Phase semantics** · Discuss
  - Context: DESIGN.md "Configuration files", "Phases and Procedure"; brownfield phase handling in `ptc_cooldown_simulator.py`.
  - Decide: How `enable` and `disable` become a per-phase graph, and how the procedure's stop condition is configured.

- [ ] **6. Configuration loading** · Auto
  - Needs: 3, 5.
  - Context: DESIGN.md "Configuration files", "Dilution refrigerator library"; `config.py`, `library/conditions.py`, `library/links/__init__.py`.
  - Do: Load a TOML configuration into a `Model` through `configure_model`, with the `LINK_TYPES` and `CONDITIONS` registries, `load_link`, and the `exchange_gas_cutoff` and `condensing_start` conditions.
  - Verified by: A toy configuration loads into the expected model; malformed configurations raise the package's errors.

- [ ] **7. Public API and CLI** · Auto
  - Needs: 2, 8.
  - Context: DESIGN.md "Public functions", "Command-Line Interface"; `__init__.py`, `__main__.py`, `output.py`, `pyproject.toml`.
  - Do: Implement `simulate_cooldown` and `save_output`, export the public functions, and add `__main__.py` and a `run-ptc` console script.
  - Verified by: The canonical run and `run-ptc` on the toy configuration produce saved output that matches an analytic or hand-computed answer.

*End of first vertical slice.*

## Links

- [ ] **8. Conduction** · Auto
  - Needs: 6.
  - Context: DESIGN.md "Dilution refrigerator library"; `library/links/conduction.py`.
  - Do: Implement conduction through a member of given material and geometry, integrating thermal conductivity over temperature.
  - Verified by: Constant-conductivity case matches the analytic heat flow; temperature-dependent case matches numerical integration.

*End of first vertical slice.*

- [ ] **9. Gas pushing** · Auto
  - Context: DESIGN.md "Gas pushing"; `library/links/gas_pushing.py`, `library/gases/`.
  - Do: Implement `gas_pushing_pm5` as specified, with the gas loader and helium data it uses.
  - Verified by: Hand-computed conductance and heat flow for a chain of nodes.

- [ ] **10. Gas-gap switch** · Discuss
  - Context: `library/links/gas_gap_switch.py`; brownfield `ptc_cooldown_simulator.py`, `thermal_link_library.json`.
  - Decide: Physics and parameters of the switch in its on and off states, and how its state follows the phase.

- [ ] **11. Exchange gas** · Discuss
  - Context: `library/links/exchange_gas.py`, `library/gases/`; brownfield `ptc_cooldown_simulator.py`, `user_gas_links.json`.
  - Decide: Physics of heat transfer through exchange gas shared by several nodes, and its parameters.

- [ ] **12. Radiation and MLI** · Discuss
  - Context: `library/links/radiation.py`, `library/links/multi_layer_insulation.py`; brownfield `ptc_cooldown_simulator.py`.
  - Decide: Radiation physics between plates and shields, how MLI modifies it, and their parameters.

- [ ] **13. Joint contact** · Discuss
  - Context: `library/links/joint_contact.py`; brownfield `ptc_cooldown_simulator.py`, `conductance_systems.json`, `thermal_conductance_calc.py`.
  - Decide: Physics and parameters of contact conductance at bolted joints.

## Cooler

- [ ] **14. Test-stand stage masses** · Human
  - Context: `docs/proposals/cooler-capacity-from-existing-curves.md` "Data needed"; `data/`.
  - Do: Add masses and materials of stage 1 and stage 2 on the test stand, and the test-stand cooldown trace they belong to.

- [ ] **15. Path-derived capacity curves** · Auto
  - Needs: 3.
  - Context: DESIGN.md "Coupled capacity maps"; the proposal's "Capacity along the measured path"; `data/`, `data/LoadMap-SN019-08-07-2026.dat`.
  - Do: Derive `Q1(T1)` and `Q2(T2)` from the test-stand cooldown using the proposal's processing steps, and store them as table-form data with their source.
  - Verified by: Derived capacity agrees with the load map where the path enters it, within heat-capacity uncertainty.

- [ ] **16. Cooler link** · Auto
  - Context: DESIGN.md "Coupled capacity maps"; `library/links/cooler.py`.
  - Do: Implement the uncoupled cooler link: `Q1` and `Q2` from the path-derived curves, and `Q_R` from regenerator conduction as specified in DESIGN.md.
  - Verified by: Heat flows match the curves at the stage nodes and the hand-computed `Q_R` at the RGP node.

- [ ] **17. Coupling bracket** · Auto
  - Needs: 24.
  - Context: the proposal's "Coupling bracket"; `library/links/cooler.py`, `simulations/validation/sn019/`.
  - Do: Implement the coupled-bracket variant of the cooler and compare it with the uncoupled run on SN019 by residual and time to end temperature.
  - Verified by: A short write-up of the comparison in `simulations/` stating which side of the error budget the difference falls on.

## Plot

- [ ] **18. Result plots** · Auto
  - Context: DESIGN.md "Data model"; `plot/result.py`.
  - Do: Plot temperatures and heat flows over time from a `Result`, with phase switches marked and measured traces overlaid when provided.
  - Verified by: Plots render for the toy configuration with and without a measured overlay.

- [ ] **19. Graph plot** · Auto
  - Context: DESIGN.md "Data model"; `plot/graph.py`.
  - Do: Draw a graph's nodes and links, labelled by name and link type.
  - Verified by: The toy configuration's graph renders with every node and link present.

## Validation

- [ ] **20. SN019 fridge inventory** · Human
  - Context: DESIGN.md "Configuration files"; brownfield `ptc_systems.json`.
  - Do: Supply SN019's nodes, masses, materials and links.

- [ ] **21. Measured-data alignment** · Discuss
  - Context: DESIGN.md "Scope"; `data/Cooldown-SN019-08-06-2026.dat`.
  - Decide: Which thermometer maps to which node, and how t = 0 is defined in the measured cooldown.

- [ ] **22. Residual metric** · Auto
  - Needs: 21.
  - Context: DESIGN.md "Scope".
  - Do: Implement the least-squares log-temperature residual between a `Result` and a measured cooldown, across all thermometers except MXP.
  - Verified by: Zero residual for identical traces; hand-computed residual for a known offset.

- [ ] **23. Legacy baseline** · Auto
  - Needs: 20, 22.
  - Context: `../pulse_tube_cooldown_simulation`; `simulations/validation/`.
  - Do: Run the legacy model on SN019 and record its residual as the v1 target.
  - Verified by: The baseline residual and the run that produced it are recorded in `simulations/validation/`.

- [ ] **24. SN019 validation** · Auto
  - Needs: all of Links, 16, 18.
  - Context: DESIGN.md "Validation and testing"; `simulations/validation/sn019/`, `tests/regression/`.
  - Do: Write the SN019 configuration, run it, and compare with the measured cooldown and the legacy baseline. Once its residual is at or below the baseline, promote it to a regression test.
  - Verified by: The validation run's residual is at or below the legacy baseline, with comparison plots, and the regression test passes.
