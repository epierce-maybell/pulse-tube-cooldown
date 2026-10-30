# Working in this repository

## Start here
1. Before writing code or proposing structure, read `docs/DESIGN.md` (the specification), your task's entry in `TASKS.md`, and the relevant ADR in `docs/adrs/`, if one exists for your task. `CONTEXT.md` is orientation only; where it disagrees with those files, they are current.
2. When sources disagree: specific ADRs, then DESIGN.md, then the task entry, then the code and its tests. Where all are silent, the code and its tests are the record.
3. Other sessions are working this checkout at the same time. Claim your task in `TASKS.md` with `[o]` before starting; mark it `[x]` only when it lands. 
4. One task per session. When it lands, hand off: the next session starts from the repository, not from this chat.

## Design authority
5. Only the human edits `docs/DESIGN.md`. Linting it is allowed.
6. Recorded decisions can be challenged, not bypassed. If logic, elegance, or pragmatism contradict DESIGN.md or an ADR, say so and propose the replacement wording; until the human records the change, build to the decision as written.
7. A decision made in chat is not recorded until it is in DESIGN.md, an ADR, or the task entry. Before a handoff, say what is still unrecorded and supply the proposed text. ADRs are for significant or controversial decisions only.

## Code
8. An `__init__.py` holds only a docstring and re-exports. Definitions live in named modules, such as `library/materials/load.py`, even where a task or DESIGN names an `__init__.py`.
9. Values that could change between systems, materials, or runs live in configuration or library data files, with their sources. Values that define a formula live in the code at their point of use. Physical values are SI.
10. Handle an error only when it could reasonably happen, would corrupt something downstream, and needs handling beyond what Python already raises. Parse at the boundary rather than validate throughout, and only where parsing solves a real problem.
11. Comments and documentation say what is, not what isn't and not why a choice was made. Defensive comments, guards against hypothetical misuse, and explanations aimed at future agents are code smells; assume later agents are as capable as you.
12. Define every symbol in a way a human can read. There are no standard conventions for physics notation.
13. Follow the Google style guides for [Python](https://google.github.io/styleguide/pyguide.html) and [documentation](https://google.github.io/styleguide/docguide/best_practices.html).

## Git
14. Commit only the files you changed in this session. Reviews and commits exclude other sessions' uncommitted work.
15. Commit messages are short, imperative, and not co-signed.

## Working with the human
16. Accept ambiguity where a human interacts with the repository; reserve rigor for agent-only surfaces.
17. Work with what you have. If a problem can be worked without a missing fact, do so and say what you assumed. If a missing fact blocks a foundational assumption, name the blocker.
