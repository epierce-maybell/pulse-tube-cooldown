# Pulse-Tube Cooldown Simulator

A ground-up rebuild of the pulse-tube cooldown model for dilution refrigerators. It predicts how long a fridge with given specs takes to cool down. The previous implementation (`../pulse_tube_cooldown_simulation`) is the reference for physics.

## Language

**Component**: Any part of the model that carries state. Every component is either a node or a link.

**Node**: A thermal mass with a temperature, e.g. a plate or a linked mass.

**Boundary node**: A node held at a fixed temperature, e.g. the 295 K surroundings (top flange and vacuum can).

**Link**: Anything that transfers heat. It takes one of three shapes: between two nodes (struts, joints, gas-gap switches, radiation); among many nodes sharing one medium (exchange gas); or a source or sink at nodes (heaters, the cooler).

**Internal state**: The time-varying quantities a component tracks for itself. For a node this includes its temperature. For a link it can include things like sorption-pump temperature or gas pressure.

**Event**: Something that happens during a run and changes the system, e.g. pumping out the exchange gas or turning a heater on. Events are first-class, alongside nodes and links. An event fires on a state condition or at a set time, and changes link parameters or modes; the topology stays fixed.

**Engine**: The generic numerical layer. It knows components, states, equations and events.

**DR library**: The physics layer. It holds the dilution-refrigerator component types (plates, pulse tube, gas-gap switches, …) and the templates that wire them into a fridge.

## Positions

- Primary question: how long does a system with given specs take to cool from 300 K until the MXP turns on (typically at 40–50 K)?
- The model reports a (time, temperature) trajectory for every node. Link heat flows are also recorded over time.
- It is a design tool for R&D first, moving toward a per-unit twin of individual fridges.
- It models dilution refrigerators, including future DR architectures that may be very different and are not yet known.
- The engine is generic, and the dilution-fridge knowledge lives in the DR library. Physics is kept separate from numerical methods.
- The architecture supports parameter sensitivities, e.g. of cooldown time to a component parameter.
- Use is one run at a time. Runs should be efficient enough to group a handful together.
- Linked masses are nodes.
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
- The engine supports arbitrary initial conditions and events, so a warm-up simulation can be added later.
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
