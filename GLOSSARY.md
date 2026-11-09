# Glossary
Shared terminology for the pulse-tube cooldown project.

## Rules for LLM agents
1. Use the terms defined in this glossary where appropriate.
2. When a project concept needs a term that is missing from the glossary, stop to discuss it with the human, then add the agreed term and definition.
3. Use common sense about what needs an entry. Focus on terms where shared project language matters.

## Language
**Model**: A general term for a representation of a system or process.

**Material**: A named substance with temperature-dependent specific heat and thermal conductivity.

**Node**: A thermal mass with a single temperature, made of materials, each with a mass.

**Boundary node**: A node held at its initial temperature.

**Link**: Anything that transfers heat between nodes, couples several nodes through a shared medium, or acts as a source or sink at nodes.

**Graph**: The nodes, the links, and which nodes each link touches.

**Phase**: A part of a procedure with a particular set of active links. A phase may have a condition that ends it and advances to the next phase.

**Procedure**: The ordered phases and the condition that ends the simulation.

**Run**: One execution of a configured simulation.

**State**: The quantities the solver integrates, chiefly node temperatures.

**Result**: Node temperatures and link heat flows over time.

**Engine**: The generic numerical layer that calculates energy balances and integrates the state.

**DR library**: The dilution-refrigerator physics layer: component types and the data that parameterizes them. It builds on the engine.

**Head**: One pulse-tube cold head.

**Archetypal head**: A representative head model used when measurements for the particular head are unavailable.

**Manifold**: The compressors and heads sharing one gas supply.

**Stage 1, stage 2, clamp**: The head's two cold stages and its regenerator clamp.

**PRP**: Primary plate.

**RGP**: Regenerator plate.

**CFP**: Cold-foot plate.

**STP**: Still plate.

**ICP**: Intermediate cold plate.

**MXP**: Mixing chamber plate.

**Load map**: A head's measured stage temperatures against stage heater powers.

**Heat flow**: The rate of heat transfer, in W. Heat flow into a node is positive when it heats the node and negative when it cools the node.

**Cooling power**: The rate at which heat is removed, in W, positive for heat removal. Removing 1 W from a node gives a heat flow into that node of -1 W.
