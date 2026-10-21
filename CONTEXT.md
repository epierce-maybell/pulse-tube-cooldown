# Pulse-Tube Cooldown Simulator

A ground-up rebuild of the pulse-tube cooldown model for dilution refrigerators. It predicts how long a fridge with given specs takes to cool down. The previous implementation (`../pulse_tube_cooldown_simulation`) is the reference for physics.

## Language

**Component**: A node or a link.

**Material**: A named substance with temperature-dependent specific heat and thermal conductivity. The DR library supplies instances from data files.

**Node**: A thermal mass with a temperature, e.g. a plate or a linked mass. It is made of materials, each with a mass.

**Boundary node**: An ordinary node that the model holds at its initial temperature, e.g. the 295 K surroundings (top flange and vacuum can).

**Link**: Anything that transfers heat. It takes one of three shapes: between two nodes (struts, joints, gas-gap switches, radiation); among many nodes sharing one medium (exchange gas); or a source or sink at nodes (heaters, the cooler). Given the temperatures of its nodes, it returns the heat flow into each of them. Sources and sinks add or remove heat; the rest move it between nodes.

**Graph**: The topology. It holds the nodes and links in order, and for each link the indices of the nodes it touches. A component's index is its identity in the engine; names are for display.

**Model**: The graph plus the conditions of a run: every node's initial temperature, which nodes are boundary nodes, the stop condition (an end node reaching its end temperature, or a maximum time) and solver options.

**State**: The quantities the solver integrates: one flat array of node temperatures, ordered by the graph. Link heat flows are computed from the state.

**Result**: The output of a run: its status, the model, and every node's temperature and every link's heat flows over time.

**Engine**: The generic numerical layer. It knows components, states and equations.

**DR library**: The physics layer. It holds the dilution-refrigerator component types (plates, pulse tube, gas-gap switches, …) and the templates that wire them into a fridge.

## Positions

- Primary question: how long does a system with given specs take to cool from 300 K until the MXP turns on (typically at 40–50 K)?
- The model reports a (time, temperature) trajectory for every node. Link heat flows are also recorded over time.
- It is a design tool for R&D first, moving toward a per-unit twin of individual fridges.
- It models dilution refrigerators, including future DR architectures that may be very different and are not yet known.
- The DR library depends on the engine, never the reverse, so the solver is testable with stand-in components and new link types need no engine changes. The data model (materials, nodes, links, graph, model) is shared by both.
- The architecture supports parameter sensitivities, e.g. of cooldown time to a component parameter.
- Use is one run at a time. Runs should be efficient enough to group a handful together.
- Linked masses are nodes.
- Temperature is the only integrated field in v1. A part with its own thermal dynamics, such as a sorption pump, is a node.
- Changes in a link's behavior during a run are state-dependent and live inside the link. Exchange-gas pressure is an empirical function of CFP temperature, zero below the pumpout cutoff (default 70 K).
- The v1 DR link types are the brownfield set: cooler, gas-gap switch, exchange gas, radiation shields, strut conduction, MLI, bolted-joint contact resistance.
- The design spec covers v1: 300 K to MXP turn-on, the 6-plate reference system, the CLI. It feeds a task list detailed enough for agents to work autonomously, without aligning with the user at each step.
- The design spec defines the cooler link's interface only. The cooler model is a later decision.
- Input sanity checks live in the config layer. Checks on the numerics live in the tests.
- Python is the language.
- Brownfield physics models carry over largely unchanged. The rebuild changes the architecture.
- The 6-plate fridge (PRP, RGP, CFP, STP, ICP, MXP, plus linked masses) and its measured cooldown traces are the reference system for validation.
- Available hardware data: the manufacturer's 2D cooler load map (`data/SN9_PTR_load_map.png`) and abundant measured cooldown data.
- Every system has thermometers (up to ~6). Their positions vary between systems.
- Calibration is reserved for cases where it is necessary. Each calibrated value absorbs uncertainty that could belong to whatever is under test.
- The brownfield's external gas-gap switch model becomes a link type in the DR library.
- The engine supports arbitrary initial conditions, so a warm-up simulation can be added later.
- Configs are written fresh in TOML. A one-time converter from brownfield saved systems is a nice-to-have.
- Material and cooler libraries are data files in the repo, each with a recorded source.
- A saved run is one folder holding the config used, the results, and metadata (code version, git commit).
- The v1 interface is the CLI. The data structures are designed to support a GUI later.
- Docstrings are short and say what code does. Derivations and equations live in `docs/`.
- Tests are unit tests per component, analytic test cases for the solver, and a regression test against a reference run.
- The CI/CD decision is deferred until after initial development.
- The implicit solver is motivated by stiffness. An off-the-shelf stiff integrator is acceptable.
- Each cooler stage flange is a node with its own mass and heat capacity.
- Using scipy to construct Jacobians is preferred over having components calculate their own partial derivatives.
- Sensitivity will be handled with external wrappers if needed later.
