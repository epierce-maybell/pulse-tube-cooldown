# Maintained validation

Run from the repository root:

```text
uv run python -m simulations.validation.head_comparison.run
uv run python -m simulations.validation.legacy_baseline
```

- `validation/head_comparison/`: selected measured cooldowns, load maps, and all
  four stage-response derivatives. `records.toml` records selection and run
  conditions; `comparison.png` and `comparison.svg` are the comparison target.
  The runner accepts a candidate prediction function with `--candidate` and
  an output stem with `--output`.
- `validation/legacy_baseline/`: retained SN09 and FP1 reference runs, inputs,
  scores, plots, and interpretation. The runner `validation/legacy_baseline.py`
  uses the sibling legacy repository identified in its inputs.
- `validation/measurements.py`: refrigerator measurement windows, channel
  mapping, and scoring used by the legacy comparison.

Validation tests live in `tests/validation/` and run with `uv run pytest`.

Exploratory models, fitting campaigns, diagnostic scripts, and their generated
outputs belong in [`../pulse-tube-research/`](../../pulse-tube-research/README.md).
The research snapshot retains the original paths of earlier head investigations.
Future refrigerator validation runs belong here once their configurations and
comparison conventions are established.
