# Pulse-Tube Cooldown Simulator

Navigation and vocabulary for a newcomer. Specifications and decisions live in `docs/`, and the work list is `TASKS.md`. Where this file disagrees with those, they are current.

## The problem

A dilution refrigerator (DR) starts a cooldown at room temperature. A pulse-tube cooler carries it from about 300 K to the point where condensing begins and the fridge is more or less operational. This project predicts how long that takes for a fridge with given specs, along with every part's temperature and heat flows along the way. It is a design tool for R&D, aimed at a per-unit twin of individual fridges, and it should accommodate future DR architectures.

The project rebuilds an earlier model, `../pulse_tube_cooldown_simulation` (branch `Advanced_Feature_SIM`), called the legacy or brownfield model. It is the reference for physics and data.

## The physical system

The reference system is a 6-plate fridge: PRP, RGP, CFP, STP, ICP and MXP, plus masses linked to those plates, inside a vacuum can at room temperature. Validation systems are SN09 (one PT1.1 head) and FP1 (two PT1.5 heads). The cooler removes heat at its first stage, second stage, and regenerator clamp, each with a capacity that depends on temperature. Heat also moves through gas-gap switches, exchange gas, radiation (reduced by multi-layer insulation, MLI), support struts, bolted joints, and heaters.

Physics lives in `docs/DESIGN.md` and the ADRs, not here.

## Language

- **Material**: a named substance with temperature-dependent specific heat and thermal conductivity.
- **Node**: a thermal mass with a temperature, made of materials, each with a mass. A **boundary node** is held at its initial temperature.
- **Link**: anything that transfers heat between nodes, several nodes sharing a medium, or sources and sinks at nodes. Given node temperatures, it returns the heat flow into each.
- **Graph**: the nodes, the links, and which nodes each link touches.
- **Model**: a graph plus the conditions of a run: initial temperatures, boundary nodes, a stop condition, and solver options.
- **State**: the quantities the solver integrates, chiefly node temperatures.
- **Result**: node temperatures and link heat flows over time.
- **Engine**: the generic numerical layer. It knows components, states, and equations.
- **DR library**: the physics layer: component types and the data that parameterizes them. It builds on the engine.
- **Head**: one pulse-tube cold head. A **manifold** is the set of compressors and heads sharing one gas supply.
- **Stage 1, stage 2, clamp**: the cooler's two cold stages and its regenerator clamp.
- **Load map**: a cooler's measured stage temperatures against stage heater powers.

## Where to work

Paths are under `src/pulse_tube_cooldown/` and `tests/` unless noted. DESIGN names are section headings in `docs/DESIGN.md`.

| Work area | Specification | Implementation | Tests | Data and records |
|---|---|---|---|---|
| Types, graph, procedure, model | DESIGN "Data model", "Phases and Procedure" | `core/`, `library/conditions.py` | `core/`, `library/test_conditions.py` | |
| Engine and solver | DESIGN "Engine" | `engine/` | `engine/` | |
| Configuration, CLI, output | DESIGN "Public API" | `config.py`, `__main__.py`, `output.py` | `test_config.py`, `test_cli.py`, `test_output.py` | `data/<system>/fridge.toml`, `tests/integration/fixtures/`, `simulations/validation/toy/` |
| Materials and fits | DESIGN "Dilution refrigerator library" | `library/materials/`, `library/property_fits.py` | `library/test_materials.py`, `test_material_data.py`, `test_property_fits.py` | `library/materials/data/` |
| Gases, surfaces | DESIGN "Dilution refrigerator library" | `library/gases/`, `library/surfaces/` | `library/test_gases.py`, `test_gas_data.py`, `test_helium.py`, `test_surfaces.py` | `library/gases/data/`, `data/gases/` |
| Conduction, joints, radiation, MLI | DESIGN "Radiation", "Multi-layer insulation", "Joint-contact resistance" | `library/links/` | `library/test_conduction.py`, `test_joint_contact.py`, `test_radiation.py`, `test_multi_layer_insulation.py` | |
| Gas-gap switch, exchange gas | DESIGN "Gas-gap heat switch", "Exchange gas" | `library/links/gas_gap_switch.py`, `exchange_gas.py` | `library/test_gas_gap_switch.py`, `test_exchange_gas.py`, `integration/test_gas_gap_switch_run.py` | `data/systems.md`, `data/gases/` |
| Gas pushing | DESIGN "Gas pushing"; ADR 02 | `library/links/gas_pushing.py` | `library/test_gas_pushing.py` | |
| Cooler (heads, drive, link) | DESIGN "Coupled capacity maps", "Cooler configuration"; ADR 03 | `library/links/cooler/`, `library/heads/` | `library/test_head.py`, `test_drive.py`, `test_cooler.py`, `test_head_data.py` | `library/heads/data/`, `data/heads/`, `data/test_stand/`, `data/test_stand.md`, `data/compressors.md`, `docs/pt_performance_overview_2026-10-30.pdf` |
| Cooler fits and cold-head experiments | ADR 03; TASKS section 2 | `simulations/validation/sn09_load_map/`, `sn09_bare_head/`, `simulations/exploration/loadmaps/` | `integration/test_sn09_load_map_fit.py` | `data/test_stand/`, `data/heads/`, `docs/papers/pulse_tube/` |
| Plots | DESIGN "Public API" | `plot/result.py` (temperature, heat flow), `plot/graph.py` | `plot/` | |
| Measurements and validation | TASKS section 4; DESIGN "Configuration files" | `measurements/` (test-stand records), `simulations/validation/measurements.py` (cooldowns and residual) | `measurements/`, `simulations/validation/test_measurements.py` | `data/systems.md`, `data/test_stand.md`, `data/<system>/` (`fridge.toml` template, `bom_by_node.csv`, cooldown trace and its `cooldown.toml`), `docs/papers/` |

A link's physics is in its DESIGN section; its TOML registration lives in `library/links/load.py`.

## Elsewhere in the repository

- `AGENTS.md`: rules for contributing.
- `docs/adrs/`: numbered decisions. `docs/proposals/`: proposals under consideration. `docs/slides/`: project overview.
- `simulations/`: targeted experiments and validation runs, each in its own directory.
- `data/README.md`: measurement and reference-data layout. `data/scripts/`: preparation tools. `docs/diagrams/`: supplied system diagrams.

## Related resources

- `../../dev/cold_control`: operational code for the refrigerators. `python/src/cold_control/routine_framework/diagnostics.py` is the diagnostics saved from each run.
- `../../dev/cc_routines`: routines for the refrigerators.
- `../../dev/he34pak`: helium properties.
- `../../data/test_stand/`: test-stand archive.
