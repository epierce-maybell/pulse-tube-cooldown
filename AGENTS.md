# Rules specific to LLMs
1. DESIGN.md is the current source of truth. LLM agents do not have permission to edit it. Linting is allowed.
2. Commit messages must be short, imperative, and not co-signed.
3. Assume CONTEXT is deprecated where it contradicts other files or a decision under consideration. It only exists to orient new agents.
4. Do not use comments and documentation to justify choices. Defensive comments are a code smell.
5. An `__init__.py` holds only a docstring and re-exports. Definitions live in named modules, such as `library/materials/load.py`. This holds even where a task or design names an `__init__.py` as the location.
6. Assume that multiple agents are operating simultaneously. Scope Plannotator reviews and git commits to only the changes made in the current session.

# Rules for contributing
1. Documentation says "what is," not "what isn't."
2. Architectural Decision Records (ADRs) record significant or controversial decisions only.
3. View previous decisions as mutable recommendations. If logic, elegance, or pragmatism contradict a decision encoded in the repository, the decision should be revisited.
4. Keep values as close to their point of use as is reasonably possible.
5. Error handling is a slippery slope into insanity. Focus on errors that could (a) reasonably happen, (b) mess something else up down the line, and (c) need special handling beyond what native Python would throw. Lean towards parse-not-validate, but make sure parsing is solving a real problem.

# Rules for interacting
1. If you are an AI, assume that all future agents are at least as smart as you are. They have the same ability to interpolate missing information, make reasonable decisions, catch and handle discrepancies, and avoid pitfalls that you do. 
2. Accept ambiguity in places where a human will interact with the repository. Reserve rigor and strict guidelines for agent-only surfaces.
3. Symbols are meaningless unless they are defined in a way that the human can understand. Do not assume that there are standard conventions for physics notation the way that we might for mathematics.
4. You probably don't need that piece of information you're fixated on. If a problem can be worked without having all the details, do so. Get as far as you can, or if a specific fact is blocking a foundational assumption, specify the actual blocker.

# External resources
We follow the Google Style Guide for [Python](https://google.github.io/styleguide/pyguide.html) and [documentation](https://google.github.io/styleguide/docguide/best_practices.html?utm_source=chatgpt.com).
