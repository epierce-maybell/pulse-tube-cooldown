# Pulse-Tube Cooldown Simulator

An orientation to the problem, the physical system, and the repository. Specifications and decisions live in `docs/`. Where this file and those documents disagree, those documents are current.

## The problem

A dilution refrigerator (DR) starts a cooldown at room temperature. A pulse-tube cooler carries it from about 300 K to the point where condensing begins and the fridge is more or less operational. How long that takes depends on the fridge: its masses, its materials, and the thermal links between its parts.

This project predicts that time for a fridge with given specs, along with every part's temperature and heat flows along the way. It is a design tool for R&D, with a longer-term aim of a per-unit twin of individual fridges. It should also accommodate future DR architectures, which may differ substantially from today's.

The project is a ground-up rebuild of an earlier model at `../pulse_tube_cooldown_simulation`, referred to as the brownfield. The brownfield is the reference for physics.

## The physical system

The reference system is a 6-plate fridge: PRP, RGP, CFP, STP, ICP and MXP, plus masses linked to those plates. The room-temperature surroundings are the top flange and the vacuum can. Measured cooldowns of this system are the validation target.

The pulse-tube cold head removes heat at three points: its first stage, its second stage, and its regenerator clamp. Each capacity depends on the temperatures at those points. The manufacturer load map characterizes the cooler near base temperature only; most of a cooldown happens far from base.

Heat moves between parts of the fridge through:

- gas-gap heat switches,
- exchange gas, which fills the can early in the cooldown and is pumped out partway through,
- thermal radiation between plates and shields, reduced by multi-layer insulation (MLI),
- conduction through support struts,
- contact resistance at bolted joints, and
- heaters.

A few features of the physics shape the model:

- The heat capacity of copper is concentrated above about 100 K, so cooler capacity at the warm end sets most of the cooldown time.
- Thermal timescales vary widely across parts and temperatures, which makes the governing equations stiff.
- Each system carries up to about six thermometers, in positions that vary between systems. They are the measurements a simulation is compared against.

## Language

**Component**: A node or a link.

**Material**: A named substance with temperature-dependent specific heat and thermal conductivity.

**Node**: A thermal mass with a temperature, such as a plate or a linked mass. It is made of materials, each with a mass.

**Boundary node**: A node held at its initial temperature, such as the room-temperature surroundings.

**Link**: Anything that transfers heat. A link connects two nodes (a strut, a joint, a switch, radiation), several nodes sharing one medium (exchange gas), or acts as a source or sink at nodes (a heater, the cooler). Given the temperatures of its nodes, it returns the heat flow into each of them.

**Graph**: The topology: the nodes, the links, and which nodes each link touches.

**Model**: A graph plus the conditions of a run: initial temperatures, boundary nodes, a stop condition and solver options.

**State**: The quantities the solver integrates, chiefly node temperatures.

**Result**: The output of a run: node temperatures and link heat flows over time.

**Engine**: The generic numerical layer. It knows components, states and equations.

**DR library**: The physics layer. It holds the dilution-refrigerator component types, such as plates, the cooler and gas-gap switches, and the data that parameterizes them. It builds on the engine.

## The repository

- `AGENTS.md`: rules for contributing and interacting.
- `docs/DESIGN.md`: the design specification.
- `docs/proposals/`: proposals under consideration.
- `docs/slides/`: an overview presentation of the project.
- `data/`: measured hardware data (cooldown traces, cooler load maps) and material properties.
- `src/pulse_tube_cooldown/`: the package.
- `simulations/`: validation experiments and other runs.

## Related resources

- `../pulse_tube_cooldown_simulation`: legacy code with similar physics. The reference version is the `Advanced_Feature_SIM` branch.
- `../../dev/cold_control`: routines and code for the refrigerators.
- `../../dev/cold_control/python/src/cold_control/routine_framework/diagnostics.py`: diagnostics saved from each run.
