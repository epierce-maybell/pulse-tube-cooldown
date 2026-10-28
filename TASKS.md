# Tasks — 40/81 done

Build the v1 simulator specified in `docs/DESIGN.md`: the six-plate systems SN09 and FP1 from room temperature until condensing begins, each with least-squares log-temperature residual comparable to or lower than the legacy model's across system thermometers except MXP.

## How to use this list

`docs/DESIGN.md` is authoritative; `AGENTS.md` governs contributions. Proposals supply discussion context. The legacy repository at `../pulse_tube_cooldown_simulation`, branch `Advanced_Feature_SIM`, supplies physics and data references. Ethan records agreed DESIGN changes; agents can propose wording during discussion.

The **first slice** is a toy thermal network through configuration, simulation, saved output, and CLI, checked against a known reference. It includes a phase switch and a fixed boundary node. It ends at 7.4; areas mark their slice membership. Discussions 5.1 and 5.2 gate this slice.

Each task has one focused review outcome and includes relevant tests. Constructed inputs allow checking a task before its production caller exists.

- **Auto**: an agent carries out the assignment unattended for subsequent review.
- **Discuss**: settle the named choices with Ethan and record the outcome; a separate implementation task follows the recorded decision.
- **Human**: supply the specific missing input after checking available sources.

Within each numbered area, work in order. **Needs** states dependencies across areas. Choose an unchecked task whose earlier tasks in its area and explicit dependencies are complete, preferring first-slice work. Mark claimed work `[o]` and completed work `[x]`. Conditional tasks complete when their condition is resolved and the outcome recorded.

Headings show completed/total tasks and active claims. Update the area and overall counts whenever task markers or the task list change.

Identifiers retain the original task number as a prefix. Context paths are relative to `src/pulse_tube_cooldown/` unless another repository location is given. Follow-up tasks receive new identifiers and concrete outcomes before implementation starts, with dependencies checked against the existing order. A discussion that exposes missing human inputs adds explicit Human dependencies to its implementation task.

## 1. [7/7] Core types and equations

- [x] **1.1 Material, gas, node, and link types** · Auto
  - Context: DESIGN "Data model" and "Dilution refrigerator library"; `src/pulse_tube_cooldown/core/types.py`.
  - Do: Define Material, Gas, Node, and the Link protocol with the specified fields and link-local heat-flow convention.
  - Verified by: Constructed components evaluate material properties and return signed flows in their declared node order.

- [x] **1.2 Graph connectivity and lookups** · Auto
  - Context: DESIGN "Data model"; `core/types.py`, `core/errors.py`.
  - Do: Define Graph and name lookups, mapping link-local node order to graph indices. Introduce errors needed by these operations.
  - Verified by: Differently ordered, overlapping links resolve correctly; invalid references produce useful errors.

- [x] **1.3 Procedure and run types** · Auto
  - Context: DESIGN "Data model" and "Phases and Procedure"; `core/types.py`.
  - Do: Define Phase, Procedure, Model, and Result, adding the `xarray` dependency. Phases carry their links; Model carries the shared nodes and builds one Graph per phase. Procedure raises ConfigurationError when it has no phases, phase names repeat, a phase other than the last lacks `until`, the last phase has `until`, or a link name refers to different nodes in different phases. Model raises ConfigurationError, naming the phase, when a phase graph is invalid.
  - Verified by: A constructed multi-phase Procedure and Model are accepted, with one graph per phase over the shared nodes; each invalid Procedure or phase graph produces a useful error; a Result holds its status, model, and Dataset.
- [x] **1.4 Node heat capacities** · Auto
  - Context: DESIGN "Engine"; `engine/equations.py`.
  - Do: Implement `calc_heat_capacity` by summing mass times specific heat over each node's composition.
  - Verified by: Single-material and mixed-material nodes match hand-computed temperature-dependent capacities.

- [x] **1.5 Link flows and node heat balances** · Auto
  - Context: DESIGN "Engine" and Link convention; `engine/equations.py`.
  - Do: Implement `calc_link_heat_flows` and `calc_net_heat_flow`, routing link-local outputs through connectivity.
  - Verified by: Networks with overlapping links, a multi-node link, and a source or sink give the expected signed node totals.

- [x] **1.6 Temperature derivatives and fixed boundaries** · Auto
  - Context: DESIGN "Engine"; `engine/equations.py`.
  - Do: Implement `calc_temperature_rate_of_change` from net heat flows and capacities, with zero derivatives at fixed boundary nodes.
  - Verified by: A small network matches hand-computed rates while retaining heat transfer to fixed boundaries.

- [x] **1.7 Jacobian** · Auto
  - Context: DESIGN "Engine"; `engine/equations.py`.
  - Do: Implement `build_jacobian`, including temperature-dependent capacities and fixed boundaries.
  - Verified by: A nonlinear network agrees with independent finite differences and has the expected fixed-boundary rows.
  - Outcome: no build_jacobian; the solver estimates the Jacobian (DESIGN "Engine").

*End of this area's first-slice work.*

## 2. [4/4] Integration and results

- [x] **2.1 Single-phase integration** · Auto
  - Context: DESIGN "Engine" and "Phases and Procedure"; `engine/solver.py`.
  - Do: Implement `integrate_phase` with solver options, phase and procedure events, time limits, and observable solver failures.
  - Verified by: Analytic relaxation integrates correctly; each termination mechanism and solver failure is observable.
  - Outcome: `integrate_phase` returns a `PhaseSolution` whose `end` names what ended the phase. `Model.solver_options` holds the `[simulation]` keys: `max_time` (end of the whole run, s) and `save_every` (sample spacing from the phase start, s) are required; `solver` defaults to "Radau"; other keys pass to `solve_ivp`.

- [x] **2.2 Procedure execution** · Auto
  - Context: DESIGN "Phases and Procedure"; `engine/solver.py`.
  - Do: Implement `run_procedure`, handing the ending state to each subsequent phase and respecting procedure termination and solver outcomes.
  - Verified by: A two-phase run switches at the expected state and continues from it; a procedure stop prevents further integration.

- [x] **2.3 Heat-flow reconstruction** · Auto
  - Needs: 1.5.
  - Context: DESIGN "Engine"; `engine/solver.py`.
  - Do: Implement `reconstruct_heat_flow` over sampled temperatures using the phase graph, producing flows on the phase graph's (sample, link, node) axes.
  - Verified by: Samples agree with direct evaluations, and unconnected link-node pairs are NaN.

- [x] **2.4 Result assembly** · Auto
  - Context: DESIGN "Data model" and "Engine"; `engine/solver.py`.
  - Do: Implement `build_result`, assembling the DESIGN Result Dataset from solutions, temperatures, heat flows, and status. Place each phase's flows on the link axis of every link used in any phase. Preserve duplicate switch times with each side's phase name.
  - Verified by: A two-phase result retains both switch-side rows with their phase names and distinct flows; a link inactive in one phase is NaN there; coordinates match graph names, `connected` matches connectivity, and attributes record status and units.
  - Outcome: `status` is the final solution's `end`; Dataset attributes hold it and that solution's solver `message`. The link axis lists every link in the procedure in order of first appearance, so links of phases that never ran are NaN. `time`, `temperature`, and `heat_flow` carry `units` attributes.

*End of this area's first-slice work.*

## 3. [7/7] Material loading

- [x] **3.1 Log-polynomial fits** · Auto
  - Context: DESIGN "Dilution refrigerator library"; `library/property_fits.py`; legacy fit definitions.
  - Do: Implement `eval_log_poly` for temperature arrays and register it in `PROPERTY_FITS`.
  - Verified by: Constant and nonconstant examples match hand calculations.

- [x] **3.2 Rational-polynomial fits** · Auto
  - Context: DESIGN property-fit registry; `library/property_fits.py`; legacy rational-polynomial definition.
  - Do: Implement and register `eval_rational_poly`.
  - Verified by: Hand-computed examples and representative copper values match independent evaluation.

- [x] **3.3 Tabulated properties** · Auto
  - Context: DESIGN property-fit registry; `library/property_fits.py`.
  - Do: Implement and register `eval_fit_from_table`, documenting interpolation and range behavior.
  - Verified by: Knots and intermediate values behave as documented; invalid tables receive useful errors.

- [x] **3.4 Material TOML loader** · Auto
  - Needs: 1.1.
  - Context: DESIGN material TOML example; `library/materials/load.py`, `library/materials/data/`.
  - Do: Implement `load_material` using `PROPERTY_FITS`. Retain sources and validity ranges in material data and document range handling.
  - Verified by: Fixtures for supported forms load and evaluate correctly; missing or malformed materials produce useful errors.
  - Outcome: `load_material(name, /, *, directory=...)` matches the `name` field across `*.toml` files. Each property table holds `form`, `range`, and an optional `source`; other keys are the fit's arguments. Outside its range, a property takes its value at the nearest end and issues a RuntimeWarning once per side; an infinite end leaves that side unbounded.

- [x] **3.7 Conductivity integrals** · Auto
  - Needs: 8.1.
  - Context: DESIGN `Material`; 8.1's outcome; `core/types.py`, `library/materials/load.py`.
  - Do: Add `conductivity_integral` to Material and build it in `load_material` as agreed in 8.1.
  - Verified by: Constant, analytic, and tabulated conductivities match their exact integrals; out-of-range values extend linearly and warn once per side.

- [x] **3.5 Copper reference material** · Auto
  - Context: DESIGN material example; `library/materials/data/`; legacy material database.
  - Do: Port OFHC copper RRR100 with property definitions, ranges, and sources.
  - Verified by: Loaded values match references at representative temperatures.

- [x] **3.6 Stainless steel reference material** · Auto
  - Context: DESIGN material example; `library/materials/data/`; legacy material database.
  - Do: Port 304 stainless steel with property definitions, ranges, and sources.
  - Verified by: Loaded values match references at representative temperatures.

*End of this area's first-slice work.*

## 4. [0/1] Validation-system materials

- [ ] **4.1 Enumerate required material ports** · Auto
  - Needs: 20.2.
  - Context: DESIGN material library; each validation system's inventory; legacy material database.
  - Do: Identify additional materials required by each validation system. Add one Auto task per material with its source and verification reference, and a Discuss task for any unsettled property representation.
  - Verified by: Each required material has a usable definition or an individually actionable port task. Add these ports as explicit dependencies of 24.1.

## 5. [3/3] Configuration decisions

- [x] **5.1 Phase activation semantics** · Discuss
  - Context: DESIGN "Configuration files" and "Phases and Procedure"; `enable` and `disable` placeholders.
  - Decide: How activation references links, carries between phases, and produces phase graphs with shared node ordering. Record a concrete configuration example.

- [x] **5.2 Run inputs and termination configuration** · Discuss
  - Context: DESIGN configuration example, Model fields, Procedure stop condition, and condition registry.
  - Decide: How TOML supplies initial temperatures, fixed boundary nodes, and the procedure stop condition. Settle enough to specify a complete toy run.

*End of this area's first-slice work.*

- [x] **5.3 Condensing start condition** · Discuss
  - Context: DESIGN "Scope" and condition registry; the operating procedure for starting condensing.
  - Decide: The physical state that starts condensing and the inputs `condensing_begins` takes.

## 6. [5/5] Configuration loading

- [x] **6.1 Fridge graph loading** · Auto
  - Needs: 1.2, 3.6, 8.2.
  - Context: DESIGN configuration and library loading; `config.py`, `library/links/__init__.py`.
  - Do: Load compositions and links into a graph. Implement `LINK_TYPES` and `load_link`, initially registering conduction and using the shared material loader.
  - Verified by: A toy fridge section yields expected nodes, parameters, and connectivity; bad references and type names produce useful errors.

- [x] **6.2 Condition factories** · Auto
  - Needs: 5.2.
  - Context: DESIGN condition example; `library/conditions.py`.
  - Do: Implement `CONDITIONS` and `gas_pushing_begins`, resolving node names against a graph.
  - Verified by: Conditions have the specified signs on either side of their thresholds and equal zero at the thresholds.

- [x] **6.3 Phase and procedure loading** · Auto
  - Needs: 1.3, 5.1.
  - Context: DESIGN "Phases and Procedure"; `config.py`; area 5 decisions.
  - Do: Build each phase's links and the Procedure from TOML using the agreed activation rules and condition factories.
  - Verified by: A two-phase fixture activates expected links and binds the expected phase and procedure conditions.
  - Outcome: `load_procedure(procedure, graph)` takes the `[procedure]` table and `load_fridge`'s full graph, and binds conditions against that graph through `load_condition` in `library/conditions.py`. A phase may disable only active links and enable only inactive ones; an unknown name fails the same check.

- [x] **6.4 Complete model configuration** · Auto
  - Context: DESIGN `configure_model` and Model fields; `config.py`.
  - Do: Implement `configure_model(path_to_config)`, assembling the `[meta]` table, procedure, initial temperatures, boundary indices, solver options, and configuration-relative data references where used.
  - Verified by: A complete toy TOML loads into the expected Model; malformed run inputs fail with useful context.
  - Outcome: `[meta]` is optional. `[simulation]` requires `initial_temperature`, a scalar or one value per node, and takes optional `fixed_boundary_nodes` names; its other keys become `solver_options`. No first-slice input references configuration-relative data.

*End of this area's first-slice work.*

- [x] **6.5 Condensing start condition** · Auto
  - Needs: 5.3.
  - Context: 5.3's recorded design; `library/conditions.py`.
  - Do: Implement and register `condensing_begins`, resolving node names against a graph.
  - Verified by: The condition has the specified sign on either side of its threshold and equals zero at the threshold.

## 7. [4/4] Public run, persistence, and CLI

- [x] **7.1 Simulation entry point** · Auto
  - Needs: 2.4.
  - Context: DESIGN "Public functions"; `src/pulse_tube_cooldown/__init__.py`.
  - Do: Implement `simulate_cooldown(model)` by executing the procedure and assembling its Result.
  - Verified by: A constructed model returns expected temperatures, switch records, heat flows, and status through the public call.
  - Outcome: `simulate_cooldown` lives in `engine/solver.py`; `__init__.py` holds no code.

- [x] **7.2 Saved output** · Auto
  - Context: DESIGN `save_output`; `output.py`.
  - Do: Implement `save_output`, writing the Result Dataset as netCDF with run provenance, and add the netCDF backend dependency.
  - Verified by: Opening the file with xarray reproduces the Result's data and attributes, including both sides of a phase switch.
  - Outcome: `save_output` lives in `output.py` and writes the Dataset with the netCDF4 engine. Provenance is `Model.meta`, the configuration's `[meta]` table, which `build_result` copies into the Dataset attributes as flat keys beside `status` and `message`.

- [x] **7.3 CLI and public exports** · Auto
  - Needs: 6.4.
  - Context: DESIGN public API and CLI; `__init__.py`, `__main__.py`, existing console-script declaration in `pyproject.toml`.
  - Do: Export the three public functions and implement `run-ptc --config ... --output ...`, including useful failures and exit status.
  - Verified by: An installed CLI configures, simulates, and saves a run; a bad configuration reports a failing exit status.
  - Outcome: `main` lives in `__main__.py`. A `PTCError`, unreadable file, or malformed TOML prints `run-ptc: error: ...` and exits 1 without saving. A `solver_failure` status saves the run, prints the solver message, and exits 1. A missing top-level table, NaN initial temperature, or unknown solver name still surfaces as a traceback.

- [x] **7.4 Canonical toy run** · Auto
  - Context: DESIGN canonical run and testing; `simulations/validation/`, integration tests.
  - Do: Add a reproducible toy configuration and reference calculation exercising conduction, a fixed boundary, and a phase switch. Run through both the three public calls and CLI.
  - Verified by: Both paths save equivalent results agreeing with the reference and preserving switch-side records.
  - Outcome: `simulations/validation/toy/` holds `toy.toml` and `run.py`. The reference in `tests/integration/test_toy_run.py` integrates the Plate's cooling time by quadrature over the library's material properties; the run ends at 21 h, before the Plate settles at the 4 K bath.

*First vertical slice complete.*

## 8. [2/2] Conduction

- [x] **8.1 Conductivity integral** · Discuss
  - Context: DESIGN `Material` and Link library; legacy `conductive_link_heat`.
  - Decide: How conduction obtains the integral of thermal conductivity between node temperatures.
  - Outcome: Each Material provides `conductivity_integral`, in W/m from the low end of the conductivity range. The loader builds it once as a cubic Hermite spline with the conductivity as its slopes; outside the range it extends linearly and warns as the properties do. Conduction models a uniform cross-section, with flow `area / length` times the integral's difference. Gas-gap brackets keep conductivity at the mean temperature.

- [x] **8.2 Conduction Link** · Auto
  - Needs: 1.1, 3.7.
  - Context: DESIGN Link library; 8.1's outcome; `library/links/conduction.py`.
  - Do: Implement conduction for a material, area, and length from the material's conductivity integral.
  - Verified by: Constant conductivity matches analytic flow; temperature-dependent conductivity matches independent integration and conserves energy under reversed temperatures.

*End of this area's first-slice work.*

## 9. [1/3] Gas data and gas pushing

- [x] **9.1 Gas loader and helium data** · Auto
  - Needs: 1.1.
  - Context: DESIGN Gas type and "Gas-gap heat switch"; `library/gases/`; legacy gas data.
  - Do: Implement the gas loader and sourced helium and nitrogen definitions for the Gas fields, using literature values.
  - Verified by: Each gas loads with expected SI properties at representative temperatures; invalid or missing definitions produce useful errors.

- [ ] **9.2 PM5 gas-pushing Link** · Auto
  - Needs: 6.1.
  - Context: DESIGN "Gas pushing"; `library/links/gas_pushing.py`.
  - Do: Implement and register `gas_pushing_pm5`, calculating adjacent-pair conductances from downstream volumes and temperatures and accumulating node flows.
  - Verified by: A hand-computed chain has expected conductances and signed flows, with conserved total energy.

- [ ] **9.3 Helium adsorption and accommodation data** · Human
  - Needs: 9.1.
  - Context: `library/gases/data/`; Corruccini, Vacuum 7–8, 19 (1959), also reproduced in Ekin, *Experimental Techniques for Low-Temperature Measurements*; 20.2's switch charcoal.
  - Do: Supply sourced Langmuir constants for He-4 on the switches' charcoal and accommodation coefficients measured across the temperature ranges of He-4 and N2. Add them to the gas files and remove the xfail marker from the He-4 data tests.
  - Verified by: He-4 loads and its data tests pass; each accommodation table cites values across its range.

## 10. [1/3] Gas-gap switch

- [x] **10.1 Switch physics and parameters** · Discuss
  - Context: DESIGN Link and phase interfaces; legacy `build_gghs_switch`, `gas_link_heat_transfer`, and thermal-link data.
  - Decide: On-state and off-state heat-transfer models, their parameters, and how phase activation selects the intended state.
  - Outcome: DESIGN "Gas-gap heat switch" records a reduced-order model whose shutoff follows from adsorption and condensation, with gas data from literature values. The switch shell is a separate conduction link from each system's inventory.

- [ ] **10.2 Switch heat transfer** · Auto
  - Needs: 9.1.
  - Context: DESIGN "Gas-gap heat switch"; `library/links/gas_gap_switch.py`.
  - Do: Implement the switch Link, solving for gas pressure from adsorption and the cold-node vapor-pressure cap, with brackets in series.
  - Verified by: Independent calculations match conductance where gas is fully free, where adsorption limits pressure, and where condensation limits it; flow direction reverses with temperature order.

- [ ] **10.3 Switch configuration and shutoff** · Auto
  - Needs: 6.4, 9.3.
  - Context: DESIGN "Gas-gap heat switch"; link registry and configuration fixtures.
  - Do: Register the switch and load it from configuration, requiring exactly two nodes.
  - Verified by: A configured run shows the switch carrying heat early and shutting off as it cools, without a phase change; a switch with other than two nodes fails to load.

## 11. [1/3] Exchange gas

- [x] **11.1 Exchange-gas physics and parameters** · Discuss
  - Context: DESIGN multi-node Link; legacy `exchange_gas_link`, pressure handling, and gas-link data; `../../dev/cc_routines` `monitor_ovc`.
  - Decide: The shared-gas heat-transfer model, pressure and geometry inputs, participating nodes, and phase activation.
  - Outcome: One multi-node link sums the gas-gap switch's gap conductance over facing surface pairs (`nodes`, `area`, `gap`), including vacuum can to PRP, with N2 standing in for air. Pressure is `monitor_ovc`'s P1(T_CFP) cubic, in mbar rather than legacy's torr, clamped at zero, so the link shuts itself off near 68 K without a phase change.

- [ ] **11.2 Exchange-gas Link** · Auto
  - Needs: 9.1, 10.2.
  - Context: 11.1's recorded design; `library/links/exchange_gas.py`.
  - Do: Implement the agreed exchange-gas transfer across participating pairs, sharing the gap conductance with the gas-gap switch.
  - Verified by: Independent small-network calculations establish signed flows that sum to zero, in both the viscous and free-molecular limits; flows vanish where the pressure curve reaches zero.

- [ ] **11.3 Exchange-gas configuration** · Auto
  - Needs: 6.4.
  - Context: 11.1's recorded design; link registry and configuration fixtures.
  - Do: Register the Link and load its pairs, gas, reference node, and pressure coefficients from configuration.
  - Verified by: A configured run carries gas heat transfer while the reference node is warm and none once the pressure reaches zero, without a phase change; a pair naming a node outside the link fails to load.

## 12. [2/6] Radiation and insulation

- [x] **12.1 Radiation physics and parameters** · Discuss
  - Context: DESIGN radiation Link; legacy `radiation_parallel_plates` and `radiation_enclosed`.
  - Decide: Supported plate/shield geometries and area, emissivity, and enclosure parameters.
  - Outcome: DESIGN "Radiation" records one gray-body link between two nodes ordered `[inner, outer]`, with per-node areas and named Surfaces whose emissivity is evaluated at each surface's own temperature. Parallel plates have equal areas. An inner area larger than the outer area raises ConfigurationError.

- [x] **12.5 Surface type and loader** · Auto
  - Needs: 1.1.
  - Context: DESIGN `Surface` and "Radiation"; `library/surfaces/load.py`, `library/surfaces/data/`; `make_fit` in `library/property_fits.py`.
  - Do: Define Surface and implement `load_surface`, building emissivity in the material property format.
  - Verified by: Fixture surfaces load and evaluate emissivity, with range handling as for materials; missing or malformed surfaces produce useful errors.
  - Outcome: `load_surface(name, /, *, directory=...)` reads the `[emissivity]` table with `make_fit`. Material, gas, and surface loaders share `read_named_toml` in `library/toml_files.py`, which picks the one `*.toml` file whose `name` field matches. `library/surfaces/data/` arrives with 12.6.

- [ ] **12.2 Radiation Link** · Auto
  - Needs: 6.1, 12.5.
  - Context: 12.1's recorded design; `library/links/radiation.py`.
  - Do: Implement and register the agreed radiation Link.
  - Verified by: Reference calculations cover parallel plates, unequal enclosed areas, temperature-dependent emissivity, equal temperatures, and reversed temperature ordering; an inner area larger than the outer area raises ConfigurationError.

- [ ] **12.6 Surface emissivity data** · Human
  - Needs: 12.5, 20.2.
  - Context: DESIGN "Radiation"; each validation system's inventory; `library/surfaces/data/`.
  - Input: Supply sourced emissivity over temperature for each radiating surface finish in SN09 and FP1, as surface files.

- [ ] **12.3 MLI physics and parameters** · Discuss
  - Context: DESIGN MLI Link; legacy MLI functions; outcome of 12.1.
  - Decide: How insulation modifies transfer, its parameters, and its relationship to the radiation Link.

- [ ] **12.4 MLI Link** · Auto
  - Context: 12.3's recorded design; `library/links/multi_layer_insulation.py`.
  - Do: Implement and register the agreed insulation heat-transfer behavior.
  - Verified by: Reference cases establish the temperature and insulation-parameter dependence, including a loaded configuration.

## 13. [0/2] Joint contact

- [ ] **13.1 Contact physics and parameters** · Discuss
  - Context: DESIGN joint-contact Link; legacy `thermal_conductance_calc.py` and `conductance_systems.json`.
  - Decide: Contact-conductance physics for bolted joints and the parameters or measured data it requires.

- [ ] **13.2 Joint-contact Link** · Auto
  - Needs: 6.1.
  - Context: 13.1's recorded design; `library/links/joint_contact.py`.
  - Do: Implement and register the agreed contact Link.
  - Verified by: Independent cases establish flow and temperature dependence, including configuration loading.

## 14. [1/3, 1 claimed] Cooler characterization inputs

- [x] **14.1 Cooler characterization requirements** · Auto
  - Context: DESIGN "Coupled capacity maps".
  - Do: Determine the inputs needed to derive Q1(T1), Q2(T2), and the regenerator term.
  - Verified by: Each term in the DESIGN capacity and regenerator equations traces to a named input.
  - Outcome: Q1(T1) and Q2(T2) need a measured cooldown path T1(t), T2(t) and each stage's heat capacity C_i(T) from its masses and materials. The regenerator term needs G_R(T) and x_R, which the cooldown path does not determine.

- [o] **14.2 Supply missing characterization inputs** · Human
  - Context: 14.1's outcome and DESIGN stage energy balance.
  - Input: Supply a cooldown path T1(t), T2(t), each stage's masses and materials, and a source for G_R and x_R.
  - Test-stand masses: mounting hardware 3.9 oz (110.6 g), of which 1.4 oz (39.7 g) is copper; radiation shield 42 oz (1.191 kg) of copper; flange 16 oz (453.6 g) of aluminum.
  - Gaps: stage assignment of each part; material of the remaining 2.5 oz (70.9 g) of mounting hardware; cold-head internal masses; cooldown path; G_R and x_R.

- [ ] **14.3 FP1 pulse-tube characterization records** · Human
  - Context: 14.1's outcome; 20.1's outcome; SN09's records in `data/sn09/pulse_tube/`.
  - Input: Locate the test-stand cooldown and load map for each of FP1's two PT1.5 pulse tubes and add them under `data/fp1/pulse_tube/`.

## 15. [0/5] Path-derived capacity curves

- [ ] **15.1 Capacity derivation choices** · Discuss
  - Needs: 14.2.
  - Context: DESIGN "Coupled capacity maps"; characterization inputs. `docs/proposals/cooler-capacity-from-existing-curves.md` supplies candidate methods.
  - Decide: Cleaning, thermometer alignment, differentiation, parasitic-load estimation, and temperature-range treatment for Q1(T1) and Q2(T2). Define evidence for assessing consistency with the matching load map.

- [ ] **15.2 Prepare characterization traces** · Auto
  - Context: 15.1's processing choices; characterization data; `simulations/validation/`.
  - Do: Implement reproducible loading and preparation of the identified stage traces.
  - Verified by: Prepared traces have agreed channels, units, alignment, and sample handling, with raw-data provenance.

- [ ] **15.3 Characterization heat loads** · Auto
  - Needs: 3.4.
  - Context: 15.1's model; test-stand inventory and materials.
  - Do: Calculate stage C_i(T) and agreed Qp_i(t), adding individually scoped material-port dependencies first where needed.
  - Verified by: Known-inventory examples match hand calculations; mass, material, and parasitic-load assumptions are recorded.

- [ ] **15.4 Derive and store stage capacity curves** · Auto
  - Needs: 3.3.
  - Context: DESIGN Q_i = -C_i(T_i) * dT_i/dt + Qp_i; prepared traces and heat loads.
  - Do: Derive Q1(T1) and Q2(T2) and save sourced table data using the agreed processing and range treatment.
  - Verified by: A synthetic trace recovers its known capacity; saved curves reproduce derived stage values.

- [ ] **15.5 Capacity consistency assessment** · Auto
  - Context: Matching load map; 15.1's assessment criteria; derived curves.
  - Do: Compare capacities with load-map evidence in their shared measured region and record uncertainties and discrepancies.
  - Verified by: A reproducible comparison establishes supported temperature ranges and identifies any specific follow-up needed before using the curves.

## 16. [0/4] Cooler Link

- [ ] **16.1 Uncoupled stage cooling** · Auto
  - Needs: 15.5, 6.1.
  - Context: DESIGN "Coupled capacity maps"; `library/links/cooler.py`.
  - Do: Implement stage heat sinks from supported Q1(T1) and Q2(T2) curves, resolving any blockers identified in 15.5 first.
  - Verified by: Each stage's signed flow matches its curve independently of the other stage's temperature.

- [ ] **16.2 Regenerator conduction inputs** · Discuss
  - Context: DESIGN Q_R = G_R(T_R) * (T_RGP - T_R), with T_R = T1 + x_R * (T2 - T1); available cooler data.
  - Decide: The source or parameterization of G_R and x_R; identify any specific human input needed before implementation.

- [ ] **16.3 Regenerator heat transfer** · Auto
  - Context: 16.2's recorded parameters; `library/links/cooler.py`.
  - Do: Add the specified independent regenerator heat sink using the interpolated regenerator temperature.
  - Verified by: Hand calculations establish Q_R and its sign while Q1 and Q2 retain their specified independent behavior.

- [ ] **16.4 Cooler configuration** · Auto
  - Context: DESIGN Link loading; link registry and configuration fixtures.
  - Do: Register the complete cooler Link and load node references, capacity data, and regenerator parameters from configuration.
  - Verified by: A configured graph evaluates expected flows at all three cooled points with identifiable data sources; two cooler links on the same plates add their flows.

## 18. [0/3] Result plots

- [ ] **18.1 Temperature plot** · Auto
  - Needs: 2.4.
  - Context: DESIGN Result; `plot/result.py`.
  - Do: Plot node temperatures with units, node labels, and phase-switch markers.
  - Verified by: A multi-phase result renders expected traces and switch positions.

- [ ] **18.2 Heat-flow plot** · Auto
  - Context: DESIGN Result; `plot/result.py`.
  - Do: Plot heat flows with link/node identities, sign convention, units, and phase-switch markers.
  - Verified by: A multi-phase example displays identifiable flows on both sides of a switch.

- [ ] **18.3 Measured-temperature overlay** · Auto
  - Needs: 22.1.
  - Context: Measurement mapping and prepared traces; `plot/result.py`.
  - Do: Overlay measured temperatures on corresponding simulated node traces.
  - Verified by: A mapped, aligned fixture places each measured series at its expected node and time.

## 19. [0/1] Graph plot

- [ ] **19.1 Graph rendering** · Auto
  - Needs: 1.2.
  - Context: DESIGN Graph; `plot/graph.py`.
  - Do: Draw nodes and links labelled by name and type, including multi-node links.
  - Verified by: A toy graph renders every node, link, and intended connection with readable labels.

## 20. [2/4] Validation case and inventory

- [x] **20.1 Select the v1 validation case** · Discuss
  - Context: DESIGN "Scope"; available repository and legacy measurements and inventories; 14.1's data identification where relevant.
  - Decide: Which six-plate system has suitable paired inventory and measurements for a legacy comparison. Record identity, evidence, data locations, and remaining gaps.
  - Outcome: Validate on SN09 and FP1, which share the six-plate architecture. SN09 has one PT1.1 pulse tube; FP1 has two PT1.5 pulse tubes on the same plates. Run the new and legacy models on both systems; legacy runs are regenerated as needed. Builds are documented for every fridge, and cold_control's `cooldown` (SN09) and `fp_cooldown` (FP1) routines define each run's procedure.
  - SN09 data: `data/sn09/cooldown_2026-06-04.dat` covers 2026-06-04 00:01 to 2026-06-06 23:00 UTC at 60 s, with PRP, RGP, CFP, STP, ICP, MXP, and DHX. MXP reaches 10 K, where condensing starts, at 63.1 h. `data/sn09/pulse_tube/` holds the test-stand cooldowns and load maps for SN09's pulse tube.
  - Gaps: FP1 cooldown record; routine logs for each run; characterization records for FP1's pulse tubes (14.3).

- [x] **20.2 Supply missing system inventory** · Human
  - Context: 20.1's outcome; DESIGN nodes, links, and configuration; each system's BoM and its reduction to node masses, `data/<system>/bom_by_node.csv`.
  - Supplied: Node masses by material for both systems, with eleven gas-gap switches in each: PRP–RGP 2 × N2, RGP–CFP 2 × N2, CFP–STP 3 × He-4, STP–ICP 2 × He-4, and ICP–MXP 2 × He-4. All eleven use the same core (500-000182), outer core (500-000183), and shell (500-000185). Each gap is 0.5 mm wide with 40 cm^2 of facing area. The shell is 0.50 mm 304 SS sheet, 25.0 mm OD and 56.4 mm long, giving a 39 mm^2 conduction cross-section over a 70 mm developed length. Switches are filled to 1 atm at room temperature and sealed with a crush tube. Each switch holds one 0.197 cm^3 (0.012 in^3 in CAD) charge of granulated coconut charcoal (part 200-000021), about 0.1 g at a packed density of 0.5 g/cm^3.

- [ ] **20.3 Supply missing validation measurements** · Human
  - Context: 20.1's outcome and trace metadata.
  - Input: For each validation system, supply the cooldown record, the cold_control routine logs for that run, the `monitor_ovc` P1 diagnostics, the P2 and P3 pressure traces through gas pushing, and missing acquisition/channel metadata. Record existing adequate inputs, including `data/sn09/cooldown_2026-06-04.dat`, by reference.

- [ ] **20.4 Measure FP1 geometry from CAD** · Human
  - Context: FP1 CAD; DESIGN radiation, exchange-gas, conduction, gas-gap switch, gas-pushing, and regenerator inputs; `data/fp1/bom_by_node.csv`.
  - Input: Record in `data/fp1/geometry.csv`, one row per value with its unit and source part:
    - Outer width, depth, height, and panel thickness of the PRP, CFP, STP, and MXP shields, and the vacuum can's inner dimensions.
    - Spacing between adjacent plates, from the vacuum can top through ICP–MXP, and the cross-section of each support type spanning each spacing.
    - Gas-gap switch core and outer-core cross-sections and lengths along the heat path, the gap width and facing area, and the sealed gas volume.
    - Internal gas volumes of the still, the ICP heat exchangers, each DHX, and the mixing chamber.
    - Each pulse tube's regenerator clamp position as a fraction of the distance from first to second stage.
  - SN09 geometry follows from these values and SN09's BoM shield masses.

## 21. [0/2] Measurement interpretation

- [ ] **21.1 Map and align measured temperatures** · Discuss
  - Needs: 20.3.
  - Context: DESIGN "Scope"; each validation system's inventory and measured cooldown.
  - Decide: For each validation system, thermometer-to-node mapping, cooldown time origin, comparison window, and invalid-reading treatment. Identify MXP for its specified residual exclusion.

- [ ] **21.2 Residual sampling convention** · Discuss
  - Context: DESIGN least-squares log-temperature residual; outcome of 21.1.
  - Decide: Sampling/interpolation and aggregation conventions for computing the same residual for the new and legacy models.

## 22. [0/2] Measured data and residual

- [ ] **22.1 Measurement preparation** · Auto
  - Needs: 21.2.
  - Context: Each validation system's measured cooldown and recorded interpretation; `simulations/validation/<identifier>/`.
  - Do: For each validation system, prepare mapped, aligned measurements with units, valid comparison samples, and provenance.
  - Verified by: A fixture with known channel mapping, invalid samples, and time origin yields expected comparison series.

- [ ] **22.2 Log-temperature residual** · Auto
  - Needs: 2.4.
  - Context: DESIGN "Scope"; 21.2's convention.
  - Do: Calculate the least-squares log-temperature residual over specified thermometers and comparison samples.
  - Verified by: Identical traces give zero; a known temperature ratio gives the hand-computed residual; MXP is excluded.

## 23. [0/2] Legacy baseline

- [ ] **23.1 Reproduce the legacy cooldowns** · Auto
  - Needs: 20.2, 20.3.
  - Context: Legacy `Advanced_Feature_SIM`; each validation system's inventory and procedure; `simulations/validation/<identifier>/`.
  - Do: Run the legacy model for each validation system and record its version, inputs, and temperature traces.
  - Verified by: Each archived run matches its system and contains enough inputs to reproduce its traces.

- [ ] **23.2 Score the legacy baseline** · Auto
  - Needs: 22.2.
  - Context: Legacy runs and prepared measurements.
  - Do: Apply the agreed residual calculation to each system's legacy traces and record its baseline and comparison coverage.
  - Verified by: Repeated scoring reproduces each baseline on the same measurements and window used for the new model.

## 24. [0/5] Full-system validation

- [ ] **24.1 Validation-system configurations** · Auto
  - Needs: 4.1 and its required material ports, 6.5, 7.4, 9.2, 10.3, 11.3, 12.4, 12.6, 13.2, 14.3, 16.4, 20.2, 20.4.
  - Context: DESIGN configuration and validation; each validation system's inventory; `simulations/validation/<identifier>/`.
  - Do: Write each validation system's TOML and reproducible runner from its inventory, procedure, and cooler data.
  - Verified by: Loaded compositions, connectivity, activation, and parameters match the supplied inventory.

- [ ] **24.2 Full cooldown execution** · Auto
  - Context: Validation-system configurations and runners.
  - Do: Execute each validation system through its cooldown window and archive inputs, version, output, and solver outcome.
  - Verified by: Each reproducible run reaches intended termination with finite temperatures and inspectable phase/flow records. An execution failure yields a specific follow-up task before progressing to comparison.

- [ ] **24.3 Baseline comparison** · Auto
  - Needs: 18.3, 22.2, 23.2.
  - Context: Completed new-model runs, measured traces, and legacy baselines.
  - Do: Score each system's new run with the shared residual convention and produce measurement overlays and baseline comparisons.
  - Verified by: A report establishes whether each system meets the v1 target and identifies temperatures/time regions dominating excess residual.

- [ ] **24.4 Resolve the validation outcome** · Discuss
  - Context: 24.3's comparison and DESIGN v1 target.
  - Decide: Accept once both systems have qualifying runs, or identify the next bounded investigation or correction for the discrepancy. Add individually scoped follow-up tasks and repeat the comparison after their completion until the target is met.

- [ ] **24.5 Regression coverage for the accepted run** · Auto
  - Context: DESIGN testing; accepted validation runs; `tests/regression/`.
  - Do: Preserve each accepted run as regression coverage with reproducible inputs and justified checks on residual and relevant behavior.
  - Verified by: The tests pass for the accepted runs and detect meaningful degradation in the validated outcomes.
