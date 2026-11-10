# Maintained validation

Run from the repository root:

```text
uv run python -m simulations.validation.head_comparison.run
uv run python -m simulations.validation.legacy_baseline.run
```

- `validation/head_comparison/`: selected measured cooldowns, load maps, and all
  four stage-response derivatives. `records.toml` records selection and run
  conditions; the comparison figure is the head acceptance target. The runner
  accepts a prediction function with `--candidate` and an output stem with `--output`.
- `validation/legacy_baseline/`: the SN09 and FP1 baseline traces, complete run
  configurations, inputs, scores, plots, and their derivation script `run.py`.
  Its README describes the retained evidence and the required legacy repository.

Refrigerator reading and log-temperature scoring live in the package's
`measurements/` modules. System comparison conventions live in
`data/<system>/cooldown.toml`. Their tests live in `tests/measurements/`; legacy
baseline tests live in `tests/validation/`. Run all tests with `uv run pytest`.

Exploratory models, fitting campaigns, diagnostic scripts, and historical
comparisons live in [`../pulse-tube-research/`](../../pulse-tube-research/README.md).
The source snapshots retain the original paths and inputs of earlier investigations.
Future refrigerator validation runs belong here once their configurations and
comparison conventions are established.
