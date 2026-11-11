# Architecture and design

The pulse-tube cooldown simulator predicts how a dilution refrigerator cools from
room temperature to the start of condensing. This document describes its
architecture, the responsibilities of the repository's parts, and the philosophy
behind the modeling approach. [CONTEXT.md](../CONTEXT.md) maps individual work
areas to code, tests, and data. [GLOSSARY.md](../GLOSSARY.md) defines shared terms.

## Purpose and scope

The simulator is a design tool for understanding cooldown time, node temperatures,
and heat flows in individual refrigerators. Different refrigerator architectures
are represented by different graphs, material inventories, and procedures.
SN09 and FP1 are the reference systems for the v1 release; the release work is
tracked in [TASKS.md](../TASKS.md).

A useful prediction follows from physics we understand. Component behavior is
checked against physical invariants and known solutions, and assembled models
are compared with measurements. Agreement with data needs an explanation in the
modeled energy balances. The legacy simulator supplies reference physics, data,
and comparison runs.

## Architecture and responsibilities

A graph represents the thermal system. Nodes hold material inventories at a
single temperature; links describe heat transfer, including shared media and
sources or sinks. A configured model combines the nodes with initial conditions,
boundary nodes, a procedure, and solver options.

| Part | Responsibility |
|---|---|
| [`core/`](../src/pulse_tube_cooldown/core/) | Shared types for physical properties, nodes, links, graphs, procedures, and results. |
| [`engine/`](../src/pulse_tube_cooldown/engine/) | Assemble node energy balances, integrate a procedure, find steady states, and reconstruct results. |
| [`library/`](../src/pulse_tube_cooldown/library/) | DR component physics, physical properties, conditions, and loaders for their configuration and data. |
| [`config.py`](../src/pulse_tube_cooldown/config.py), [`__main__.py`](../src/pulse_tube_cooldown/__main__.py), [`output.py`](../src/pulse_tube_cooldown/output.py) | Build configured models, provide the command-line entry point, and save results. |
| [`measurements/`](../src/pulse_tube_cooldown/measurements/) | Read test-stand and refrigerator records into labeled datasets and compare measured and simulated cooldowns. |
| [`plot/`](../src/pulse_tube_cooldown/plot/) | Draw graphs and plot simulation histories and comparisons with measurements. |

The engine operates on the shared node and link interfaces. The DR library
supplies the physical laws and data behind those interfaces. This separation
lets a new link participate in an existing graph and solver. Models can also be
built directly in Python with custom materials and links.

A link receives temperatures only for its own nodes, in its declared order. It
returns the heat flow into each node and any additional heat capacity sharing
that node's temperature, such as attached gas. The engine adds these terms to
each node's material heat capacity and energy balance. Names identify nodes and
links in configuration and results; the graph maps link-local order to the
engine's array indices.

The current engine integrates node temperatures. Boundary nodes retain their
initial temperatures while links continue to exchange heat with them. The
adopted head model carries seven material-piece energy entries, specified in
[ADR-06](adrs/06-reduced-head-dynamics.md). Their engine
integration is task 3.4.

## How a run proceeds

The public entry points are `configure_model`, `simulate_cooldown`, and
`save_output`. The `run-ptc` command uses the same sequence.

1. **Configure.** Read the TOML configuration, load material and link data, build
   the graph, and resolve named conditions and boundary nodes. Link loaders can
   also define nodes required by their component model.
2. **Build the procedure.** Each phase selects its active links. Configuration
   applies link enables and disables cumulatively; direct Python construction
   supplies each phase's links explicitly. Each phase has a graph over the same
   nodes.
3. **Integrate.** The engine assembles energy balances from the active graph and
   passes temperature rates to SciPy's ODE solver. A phase-ending condition
   advances to the next phase with the current temperatures. A procedure-ending
   condition, the time limit, or a solver failure ends the run.
4. **Build the result.** Reconstruct link heat flows at the sampled temperatures
   and collect the phase histories in a labeled xarray dataset. Phase switches
   retain samples on both sides so changes in heat flow remain visible.
5. **Save and inspect.** Save the dataset as netCDF and use the plotting functions
   to inspect temperatures, signed heat flows, and comparisons with measurements.

Heat flow is positive into the receiving node. Results identify nodes and links
by name, carry units and supplied run metadata, and distinguish inactive links
from links carrying zero heat. The detailed result schema and interfaces are
documented with the [core types](../src/pulse_tube_cooldown/core/types.py).

Maintained validation runs group their input configuration, additional run inputs, results,
figures, and interpretation in a run directory. Run scripts manage this
arrangement and record the inputs and provenance needed to understand a result.
The output function writes the result dataset to the requested file. The
[toy integration test](../tests/integration/test_toy_run.py) checks this sequence
against a conduction reference using its [fixture](../tests/integration/fixtures/toy.toml).

The engine also provides a
[steady-state solver](../src/pulse_tube_cooldown/engine/steady_state.py) that
balances net heat flow at free nodes while holding boundary temperatures and
applying specified heater powers. It shares the transient solver's energy-balance
assembly and supports component investigations and load-map comparisons.

## Physical modeling approach

Nodes combine the masses and temperature-dependent specific heats of their
materials. Links represent the processes coupling those thermal masses. Keeping
individual contributions explicit makes heat-flow budgets available for
investigation and comparison with the physical system.

| Process | Approach |
|---|---|
| Conduction | Integrate a material's conductivity between endpoint temperatures for a solid heat path. |
| Joint contact | Use a force-dependent contact conductance with a temperature power law and sourced coefficients. |
| Radiation | Exchange gray-body radiation between enclosing surfaces, with each surface's temperature-dependent emissivity. |
| Multi-layer insulation | Use the Lockheed correlation for spacer conduction and inter-layer radiation over the blanket's area. The blanket replaces direct radiation over that area. |
| Gas-gap heat switch | Combine molecular and viscous gas conduction with bracket resistance. Determine sealed-gas pressure from free gas, charcoal adsorption, and saturation at the colder surface. |
| Exchange gas | Couple facing surface pairs through gas at a shared pressure. Use the configured pumping threshold as a pressure proxy for the reference systems, as recorded in [ADR-05](adrs/05-exchange-gas-pressure.md). |
| Gas pushing | Average gas transport over the routine's pressure cycles, accounting for the ordered gas volumes and operational limits described in [ADR-02](adrs/02-gas-pushing-routine.md). |

Short explanations and component interfaces live beside their implementations
in [`library/links/`](../src/pulse_tube_cooldown/library/links/). ADRs record
significant modeling decisions and their assumptions.

Materials, gases, and surfaces have named library data files.
Published source coefficients used to derive a library property live with that
property's data. Head designs, with their source geometry and material filling,
are packaged in `library/links/pulse_tube/data/heads/`.
Property loaders construct temperature-dependent functions and integrals from
sourced fits or tables. Validity ranges accompany the data; values outside them
use the endpoint behavior and warnings described by the loaders. Real-helium
properties for head-model work are supplied through
[`library/gases/helium.py`](../src/pulse_tube_cooldown/library/gases/helium.py).

Physical inputs use SI units. Values that depend on the system, material, or run
belong in configuration or data, with their sources recorded nearby. Separating
shared component data from individual run inputs lets comparisons identify what
changed between configurations.

## Repository organization

| Location | Contents and role |
|---|---|
| `src/pulse_tube_cooldown/` | The reusable simulator, component library, measurement readers, and plotting functions. Packaged data live with their library loaders. |
| `tests/` | Core, engine, component, interface, integration, and acceptance tests. [Test policy](../tests/README.md) describes their responsibilities and authorship rules. |
| `data/` | Source measurements, hardware records, geometry, and curated physical inputs. [The data map](../data/README.md) identifies their locations and provenance. |
| `simulations/` | Maintained validation: the head-comparison target and legacy refrigerator baseline, with their selected inputs, retained results, and findings. |
| `../pulse-tube-research/` | Local, unversioned workspace for exploratory models, fits, diagnostics, and generated research outputs. |
| `docs/` | Architecture, ADRs, proposals, reference papers, diagrams, and reporting material. Archived ADRs retain earlier decisions. |
| `CONTEXT.md`, `GLOSSARY.md`, `TASKS.md`, `AGENTS.md` | Detailed navigation, shared language, current work, and contribution rules. |

Library code reads packaged data by name and receives everything else by path.
Data a library user needs without this repository are packaged in a `data/`
directory beside the module that loads them: material, gas, and surface
properties, archetypal component data such as head designs, the head model's
shared parameters, and an archetypal head's calibration, and model defaults. Data specific to one system, serial-numbered
unit, or run live in the repository's `data/`: measurements, hardware records,
system configurations, and fits to an individual unit. Derived files record
their inputs and the command that regenerates them.

Source measurements and prepared physical inputs live in this repository. Fitted
candidates and exploratory results live in `../pulse-tube-research/`. For example,
test-stand readers live in the package, while the head comparison chooses its records and analyzes
predictions under `simulations/`. Refrigerator comparison windows, channel
mappings, and stuck-reading thresholds are recorded with the system data.
Reusable readers and log-temperature scoring live in the measurements package.
The helium adsorption fitting script reads source coefficients from packaged gas
data; its [command](../data/README.md#helium-adsorption-fit) is recorded in the data guide.
Historical BoM, CAD, and data-export preparation is preserved in
`../pulse-tube-research/archive/`.

Acceptance tests express required physical behavior and public interfaces.
Component tests and assembled simulations provide complementary evidence about
the model. Current coverage and unfinished validation work are tracked in
[TASKS.md](../TASKS.md).

## Open decisions

The adopted pulse-tube head formulation is a reduced nonlinear dynamical system
model with seven internal energy states, stage and clamp temperatures as inputs,
and cooling powers as outputs. [ADR-06](adrs/06-reduced-head-dynamics.md)
specifies its vector field, output map, calibration, and algebraic clamp
relation. Each stage's warm correction has the fixed form
`A_i * (T_i / T0)^b_i`, where `T_i` and `T0` are the stage and warm-flange
temperatures in K, `A_i` is an amplitude in W, and `b_i` is dimensionless.
The four warm coefficients, two regenerator exchange coefficients, and one
internal-energy coupling coefficient are calibrated separately per head and
drive condition, with the functional form held fixed. ADR-06 records their
interpretation and the six-head forward calibration results.
In the reference fridge, stage 1 connects to PRP, stage 2 to CFP, and
the clamp to RGP.
The stand condition sets clamp load to zero and solves its free temperature.
Production implementation and engine integration remain work for task 3.4.

The retired cooler implementation, head-design loader, generated inventories,
and experiments from [archived ADR-03](adrs/archive/03-coupled-capacity-maps.md),
together with subsequent candidate investigations, are preserved in
`../pulse-tube-research/`, with original paths and a hashed snapshot manifest.
The replacement work and the archetypal head are tracked in [TASKS.md](../TASKS.md).
