# Tasks

## Current objective

1. Release a v1.0.0 cooldown model, including physics-based acceptance tests for individual components and transient simulations plotted against real fridge data. The objective is to get a reasonable answer (read: better fit against time-to-condensing and per-stage temperature profiles) for the right reasons (physics that we understand).

2. De-slop the repository. Clean up the directory structure so that the codebase records what exists currently, rather than the mess of decisions and tests that we took to get here. Build a test suite that enforces behavior, not implementation. Document the code, including how to use/modify/extend it. Record major decisions, upgrades, and overall findings in a slide deck. 

## Instructions for LLM agents

There are three kinds of entry. If a task does not specify, default to **Discuss**. 
1. **Build**: first, ensure alignment between the agent and the user. Discuss any outstanding points. Then build and hand back the work for review. Iterate as needed.
2. **Discuss**: the output of the task is a decision, more detailed follow-up tasks, an ADR, or other documentation. Ensure that the result is produced in collaboration with the user; a decision that the user cannot defend and take ownership over is a failure state.
3. **Human**: work that agents cannot do with what they have access to on this machine. Help the user identify where to track down the information. Remember that doing so takes time and effort; focus on finding the minimum set of necessary outside input, rather than comprehensive details that may not be needed.

Claim tasks with [o], mark finished tasks with [x], and update each section to have a running count of open vs. completed items, like (5/10). If a task was started and then later put on hold, mark it as [-] with a note for the blocker and do not update progress. If a task is unclear or has unsettled decisions, interview the user relentlessly until you are aligned on all decisions, before doing the work. Do not use this file as a dumping-ground for notes or reporting. Agents may not add tasks to the list.

## Progress `[##########################----------]` 72.2% (26/36)

## 1. Rebuild test suite (13/13)
#### [x] 1.1 Public/user interface tests
#### [x] 1.2 End-to-end acceptance tests
#### [x] 1.3 Pulse tube head acceptance tests
#### [x] 1.5 Conduction acceptance tests
#### [x] 1.6 Radiation acceptance tests
#### [x] 1.7 Joint-contact resistance acceptance tests
#### [x] 1.8 MLI acceptance tests
#### [x] 1.9 Gas-gap heat switch acceptance tests
#### [x] 1.10 Exchange gas acceptance tests
#### [x] 1.11 Gas pushing acceptance tests
#### [x] 1.12 Material properties acceptance tests
#### [x] 1.13 Test policy doc, guardrails to protect tests from agents
#### [x] 1.14 Prune existing test suite

## 2. Reconcile public documentation (4/5)
#### [x] 2.1 Add GLOSSARY.md
#### [x] 2.2 Revise DESIGN.md
#### [x] 2.3 Clean up directory structure
#### [x] 2.4 Prune data/, docs/, and simulations/
#### [o] 2.5 Document how to use/modify/extend

## 3. Pulse tube model (5/5)
#### [x] 3.1 Retire ADR-03
#### [x] 3.2 Plot pulse tube data
#### [x] 3.3 Derive reduced pulse tube model
#### [x] 3.4 Implement dynamical systems approach
#### [x] 3.5 Build archetypal pulse tube head

## 4. Transient simulations (4/6)
#### [x] 4.1 Build canonical configuration file for SN09
#### [x] 4.2 Build canonical configuration file for FP01
#### [x] 4.3 Write runner for transient simulations
#### [x] 4.4 Run transient simulations, record results
#### [o] 4.5 Build cal/val dashboard
#### [ ] 4.6 Calibrate **honestly**, revise as needed

## 5. Reporting (0/4)
#### [ ] 5.1 Comparison with legacy model: audit for good result *and* right reasons
#### [ ] 5.2 Current state of the PTC performance
#### [ ] 5.3 Model architecture and lessons learned
#### [ ] 5.4 Schedule deep-dive

## 6. Release mechanics (0/3)
#### [ ] 6.1 Audit public interface
#### [ ] 6.2 Rename functions to avoid calc_ and get_ filler
#### [ ] 6.3 Prepare release, including tag, changelog, etc.
