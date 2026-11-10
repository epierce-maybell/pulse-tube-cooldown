# Pulse-Tube Cooldown Simulator

Navigation for a newcomer. [DESIGN.md](docs/DESIGN.md) describes architecture and
modeling approach; ADRs record significant decisions. [GLOSSARY.md](GLOSSARY.md)
defines shared terminology and its interaction rules. [TASKS.md](TASKS.md) tracks
current work. Those records take precedence over this navigation map.

## The problem

This project predicts cooldown time, node temperatures, and heat flows for a
dilution refrigerator, from room temperature to the start of condensing. It is a
design tool for R&D, aimed at individual fridges and future DR architectures.

The reference systems have six plates: PRP, RGP, CFP, STP, ICP, and MXP, with linked
masses inside a room-temperature vacuum can. SN09 has one PT1.1 head; FP1 has two
PT1.5 heads. The head's first stage, second stage, and regenerator clamp connect
to PRP, CFP, and RGP respectively. Other heat paths include support conduction,
bolted joints, radiation, multi-layer insulation, gas-gap switches, exchange gas,
and gas pushing.

The legacy simulator is `../pulse_tube_cooldown_simulation`, branch
`Advanced_Feature_SIM`. It supplies reference physics, data, and comparison runs.
The replacement head formulation is in progress; its current discussion is in
[the head-model proposal](docs/proposals/pulse-tube-model-direction.md).

## Where to work

Implementation and test paths below are relative to `src/pulse_tube_cooldown/`
and `tests/`. Data and experiment paths are relative to the repository root.
Detailed component explanations live beside their implementations or in the
linked ADRs.

| Work area | Architecture and decisions | Implementation | Tests | Data and records |
|---|---|---|---|---|
| Types, graph, procedure, model | DESIGN "Architecture and responsibilities", "How a run proceeds" | `core/`, `library/conditions.py` | `core/`, `library/test_conditions.py` | |
| Energy balances, transient and steady-state solvers | DESIGN "Architecture and responsibilities", "How a run proceeds" | `engine/` | `engine/` | |
| Configuration, CLI, output | DESIGN "How a run proceeds" | `config.py`, `__main__.py`, `output.py` | `test_config.py`, `test_cli.py`, `test_output.py`, `acceptance/test_interfaces.py` | `simulations/validation/toy/`, `tests/integration/fixtures/`, `data/<system>/fridge.toml` |
| Materials and property fits | DESIGN "Physical modeling approach" | `library/materials/`, `library/property_fits.py` | `library/test_materials.py`, `test_material_data.py`, `test_property_fits.py` | `library/materials/data/` under the package |
| Gases and surfaces | DESIGN "Physical modeling approach" | `library/gases/`, `library/surfaces/` | `library/test_gases.py`, `test_gas_data.py`, `test_helium.py`, `test_surfaces.py` | Packaged gas and surface data; `data/gases/` source inputs |
| Conduction, joint contact, radiation, MLI | DESIGN "Physical modeling approach"; module docstrings | `library/links/` | `library/` component tests; `acceptance/` physical requirements | `data/systems.md` |
| Gas-gap switches and exchange gas | DESIGN "Physical modeling approach"; [ADR-05](docs/adrs/05-exchange-gas-pressure.md) | `library/links/gas_gap_switch.py`, `exchange_gas.py` | `library/test_gas_gap_switch.py`, `test_exchange_gas.py`, `integration/test_gas_gap_switch_run.py` | `data/systems.md`, `data/gases/` |
| Gas pushing | [ADR-02](docs/adrs/02-gas-pushing-routine.md) | `library/links/gas_pushing.py` | `library/test_gas_pushing.py` | `data/systems.md` |
| Head model and shared drive | DESIGN "Open decisions"; [current discussion](docs/proposals/pulse-tube-model-direction.md); existing implementation follows [archived ADR-03](docs/adrs/archive/03-coupled-capacity-maps.md) | `library/links/cooler/`, `library/heads/` | `library/test_head.py`, `test_drive.py`, `test_head_data.py`, `test_cooler.py` | `data/heads/`, `data/test_stand/`, `data/test_stand.md`, `data/compressors.md` |
| Head comparisons and investigations | Head-model discussion and each experiment's records | `simulations/validation/head_comparison/`, `test_stand_head/`, `test_stand_steady_state/`, `simulations/exploration/loadmaps/` | Tests alongside the experiments | `data/test_stand/`, `docs/papers/pulse_tube/` |
| Earlier head experiments | [Archived ADR-03](docs/adrs/archive/03-coupled-capacity-maps.md) | `simulations/validation/sn09_load_map/`, `sn09_bare_head/` | Historical experiment tests and outputs | Saved inputs and findings beside the experiments |
| Plots | DESIGN "Architecture and responsibilities" | `plot/result.py`, `plot/graph.py` | `plot/` | |
| Measurements and refrigerator validation | DESIGN "Repository organization"; system comparison conventions | `measurements/`, `simulations/validation/measurements.py`, `legacy_baseline.py` | `measurements/`, `simulations/validation/test_measurements.py`, `test_legacy_baseline.py` | `data/<system>/cooldown.toml`, cooldown traces, `data/systems.md`, `simulations/validation/legacy_baseline/` |

Link types are registered in `library/links/load.py`; condition factories are
registered in `library/conditions.py`. [Test policy](tests/README.md) governs work
on the test suite. Current completion status belongs in `TASKS.md`.

## Elsewhere in the repository

- `AGENTS.md`: contribution rules and documentation responsibilities.
- `docs/adrs/`: significant decisions; `archive/` retains earlier decisions.
- `docs/proposals/`: ongoing discussions. `docs/slides/`: reporting material.
- `simulations/`: experiments and validation runs with their code and outputs.
- `data/README.md`: measurement and input-data map.
- `data/scripts/README.md`: data-preparation tools and workflows.
- `docs/papers/`, `docs/diagrams/`: reference literature and supplied diagrams.

## Related resources

- `../../dev/cold_control`: operational code for the refrigerators.
  `python/src/cold_control/routine_framework/diagnostics.py` defines saved diagnostics.
- `../../dev/cc_routines`: refrigerator routines.
- `../../dev/he34pak`: helium properties.
- `../../data/test_stand/`: the full test-stand archive.
