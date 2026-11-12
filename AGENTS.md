# Working in this repository

*Your job is to care deeply about the code. Treat design and architecture as fundamental to correctness, not optional polish. Prioritize clarity, simplicity, and maintainability. Make thoughtful decisions, build coherent abstractions, and leave the codebase better than you found it.*

## Start here
1. Before writing code or proposing structure, read your task's entry in `TASKS.md`, the sections of `docs/DESIGN.md` it touches, and any applicable ADRs in `docs/adrs/`. Read `GLOSSARY.md` for shared terminology and follow its interaction rules. `CONTEXT.md` maps each work area to its specification, code, tests, and data. Read further when the work crosses into another area or raises a specific question. `CONTEXT.md` is navigation only; where it disagrees with those files, they are current. 

## Design authority
2. In general, only the human edits markdown files. Linting is allowed. Editing is allowed when given explicit permission.
3. Recorded decisions can be challenged, not bypassed. If work on a task raises new information or if the logical, elegant approach would contradict DESIGN.md or an ADR, raise that as a discussion with the human. Only decisions that are consequential, non-obvious, and difficult to reverse are recorded in an ADR.

## Code
5. An `__init__.py` holds only a docstring and re-exports. Definitions live in named modules, not `__init__.py` files. Names have no leading underscore; `__all__` lists the agreed public surface.
6. Keep the directories organized. Caring about the codebase extends to caring about the repository structure. That is the entry-point and logical map for both humans and agents, so it needs to be clear and simple to navigate.
   Exploratory models, fits, diagnostics, and generated research outputs belong in `../pulse-tube-research/`. This repository holds the reusable library in `src/`, its tests and fixtures in `tests/`, source and curated inputs in `data/`, and maintained validation runs in `simulations/`. Historical extraction scripts live in the research archive. Promote research into the library against the agreed design and acceptance criteria.
7. Values that could change between systems, materials, or runs live in configuration or library data files, not in magic numbers or module-level CONSTANTS. Do not treat a chosen default as an invariant. Physical values are always in SI units. A source entry in a table lets a reader look up *consequential* or *controversial* values; cite the publication, file, or drawing the value is read from, and nothing else.
8. Handle an error only when it could reasonably happen, would corrupt something downstream, and needs handling beyond what Python already raises. Parse at the boundary rather than validate throughout, and only where parsing solves a real problem.
9. Comments and documentation say what is, not what isn't, and not why a choice was made. Defensive comments, guards against hypothetical misuse, and explanations aimed at future agents are code smells. Assume later agents are as capable as you.
10. Follow the Google style guides for [Python](https://google.github.io/styleguide/pyguide.html) and [documentation](https://google.github.io/styleguide/docguide/best_practices.html).
11. Use modern Python idioms and syntax. The code you were trained on is the median quality of what is publicly available; we strive for far higher quality than the median. 

## Git
12. Commit only the files you changed in this session. Reviews and commits should exclude other sessions' uncommitted work.
13. Commit messages are short, imperative, and not co-signed.

## Working with the human
14. Define every symbol in a way a human can read. There are no standard conventions for physics notation.
15. Accept ambiguity where a human interacts with the repository; reserve rigor for agent-only surfaces.
16. State assumptions up front, so that the human has an opportunity to discuss and revise the mental model with you.
17. Be economical. Use cheap sub-agents when sub-agents are called for. Trust results from other sessions, rather than wasting tokens verifying statements from the human or other agents. Never start parallel processes on the machine without asking the human how many cores you may use.
18. When you edit library code, open a Plannotator review session for the human to review changes and provide feedback. Include only your changes in the review.
