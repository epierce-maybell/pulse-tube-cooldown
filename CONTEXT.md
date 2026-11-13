# Pulse-Tube Cooldown Simulator

Navigation for a newcomer. [DESIGN.md](docs/DESIGN.md) describes architecture and
modeling approach; ADRs record significant decisions. [GLOSSARY.md](GLOSSARY.md)
defines shared terminology. [TASKS.md](TASKS.md) tracks current work. Those records
take precedence over this navigation map.

## The problem

The simulator predicts cooldown time, node temperatures, and heat flows for a
dilution refrigerator, from room temperature to the start of condensing. SN09
has one PT1.1 head; FP1 has two PT1.5 heads. Stage 1, stage 2, and the clamp
connect to PRP, CFP, and RGP respectively. Other heat paths include conduction,
bolted joints, radiation, MLI, gas-gap switches, exchange gas, and gas pushing.

The adopted head formulation is recorded in
[ADR-06](docs/adrs/06-reduced-head-dynamics.md), including its fixed power-law
warm response, seven fitted coefficients, and clamp temperature and
cooling-power closure. The head, its fit, and its packaged parameters are in
`library/links/pulse_tube/empirical/`. Head designs are in
`library/links/pulse_tube/data/heads/`. Full refrigerator configurations are
tasks 4.1 and 4.2.

## Where to work

Implementation and test paths are relative to `src/pulse_tube_cooldown/` and
`tests/` respectively. Other paths are relative to the repository root.

| Work area | Architecture and decisions | Implementation | Tests and inputs |
|---|---|---|---|
| Types, graph, procedure | DESIGN "Architecture and responsibilities", "How a run proceeds" | `core/`, `library/conditions.py` | `core/`, `library/test_conditions.py` |
| Energy balances and solvers | DESIGN "How a run proceeds" | `engine/` | `engine/` |
| Configuration, CLI, output | DESIGN "How a run proceeds" | `config.py`, `__main__.py`, `output.py` | `test_config.py`, `test_cli.py`, `test_output.py`, `acceptance/test_interfaces.py`; `integration/fixtures/toy.toml` |
| Materials and property fits | DESIGN "Physical modeling approach" | `library/materials/`, `library/property_fits.py` | Material and fit tests under `library/`; packaged material data |
| Gases and surfaces | DESIGN "Physical modeling approach" | `library/gases/`, `library/surfaces/` | Component and data tests under `library/`; packaged property data and source coefficients |
| Conduction, contacts, radiation, MLI | DESIGN "Physical modeling approach" | `library/links/` | `library/`, `acceptance/`; `data/systems.md` |
| Gas-gap switches and exchange gas | [ADR-05](docs/adrs/05-exchange-gas-pressure.md) | `library/links/gas_gap_switch.py`, `exchange_gas.py` | Component, integration, and acceptance tests; `data/systems.md` |
| Gas pushing | [ADR-02](docs/adrs/02-gas-pushing-routine.md) | `library/links/gas_pushing.py` | Component and acceptance tests; `data/systems.md` |
| Pulse-tube head and its measurements | [ADR-06](docs/adrs/06-reduced-head-dynamics.md), [archived ADR-03](docs/adrs/archive/03-coupled-capacity-maps.md) | `library/links/pulse_tube/empirical/` and its packaged `data/`, head designs in `library/links/pulse_tube/data/heads/`, `measurements/test_stand.py` | `measurements/`; `data/test_stand/`, `data/test_stand.md`, `data/compressors.md` |
| Plots | DESIGN "Architecture and responsibilities" | `plot/` | `plot/` |
| Refrigerator measurements and validation | [Simulation guide](simulations/README.md) | `measurements/refrigerator.py`, `measurements/comparison.py`, `simulations/validation/` | `measurements/test_refrigerator.py`, `validation/`; `data/<system>/cooldown.toml` and measured traces |
| Helium adsorption fit | [Data guide](data/README.md) | `library/gases/prepare.py` | Packaged `library/gases/data/he4.toml` source coefficients |

Link types are registered in `config.py`; condition factories are
registered in `library/conditions.py`. [Test policy](tests/README.md) governs the
test suite. All maintained tests live under `tests/`.

## Research and historical work

Exploratory models, fits, diagnostics, and their outputs live in the local sibling
directory [`../pulse-tube-research/`](../pulse-tube-research/README.md). It is
unversioned. Its `snapshot.json` identifies the source commit and SHA-256 hashes
of the files preserved during tasks 2.3 and 2.4. It retains the original directory
layout, source-data copies, documentation, and package source at transfer.

The former `simulations/exploration/` and head investigations under
`simulations/validation/{sn09_load_map,sn09_bare_head,test_stand_head,test_stand_steady_state}/`
are there. `archive/e40e98a-directory-audit/` preserves the subsequent preparation
scripts, retired catalogs, head definitions, and legacy comparisons, together
with their committed inputs and a SHA-256 manifest. Several older scripts
require a retired cooler implementation from Git history.

The accepted fixed-form head investigation is in
[`parametric_warm_response/`](../pulse-tube-research/simulations/exploration/parametric_warm_response/README.md).
Its `crossing_fits/` directory contains the final per-head coefficients,
six-head results, and comparison figure. Its `diagnostics/` directory records
the S0035 transient-shape investigation. ADR-06 describes the calibration and
the scope of that evidence.

The production repository retains the head-comparison figure and legacy baseline.
Research is promoted into the library against the agreed design and acceptance
criteria. New investigations belong in the research directory.

## Related resources

- `../pulse_tube_cooldown_simulation`, branch `Advanced_Feature_SIM`: legacy
  reference physics, data, and comparison runs.
- `../pulse-tube-model`: earlier exploratory model and its findings.
- `../../dev/cold_control`: refrigerator operational code; saved diagnostics are
  defined in `python/src/cold_control/routine_framework/diagnostics.py`.
- `../../dev/cc_routines`: refrigerator routines.
- `../../dev/he34pak`: helium properties.
- `../../data/test_stand/`: full test-stand archive.
