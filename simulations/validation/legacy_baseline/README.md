# Legacy baseline

The legacy model run on SN09 and FP1 and scored with `measurements.residual`. `../legacy_baseline.py` makes every
file here except this one and `inputs.toml`; run it from `simulations/validation/`.

## Runs

| Run | Configuration | Pooled RMS |
|---|---|---|
| `sn09_tuned` | `Fridge` as saved: its masses and `PT1.1 (1X F-70)`, the configuration legacy was tuned to SN09 with | 0.459 |
| `sn09_corrected` | `Fridge` with SN09's BOM masses and `PT1.1`, two F-70s | 0.806 |
| `fp1` | `Fridge` with FP1's BOM masses and two `PT1.5 (1X F-70)` | 0.495 |

Per-node scores are in `scores.toml`; `sn09.png` and `fp1.png` plot each system's runs over its measured cooldown. `sn09_tuned` matches SN09 through two offsetting differences from the system:
`Fridge`'s CFP and STP masses are about half the BOM's, and its cooler is the one-compressor estimate though SN09 runs
two. With `Fridge`'s masses and two compressors, SN09 scores 0.660; with the BOM masses and two, 0.806. FP1 with
`Fridge`'s masses scores 0.918.

## Legacy version

`Advanced_Feature_SIM` at the commit in `scores.toml`, the only legacy line with the six-plate model. The legacy
`main` branch, later in date, has the four-stage model and corrects the Lockheed MLI solid-conduction term
(`9f563c5`), which `Advanced_Feature_SIM` under-predicts by about the mean temperature.

## Run conditions

- Every node starts at the mean of the scored nodes' first valid readings, which the legacy engine also uses as the
  ambient temperature.
- Every legacy target is 0 K, so each run ends at `max_time`, one hour past the window, and DHX is not held at its
  20 K target.
- The configurations in the `.json` files load in the legacy GUI; they keep `Fridge`'s targets.

## Comparison with the slide deck

The last slide of `docs/slides/pulse_cooldown_simulator_slides.pdf` plots legacy curves over measured curves that
follow SN09's record through its first 40 h. `Fridge` as saved, from 295 K, reproduces the legacy curves within a few
kelvin for PRP, RGP, CFP, and STP. ICP and MXP end
3–5 K colder at 40 h. DHX cools faster than on the slide after about 15 h: 173 K at 20 h against about 210 K, and
57 K at 40 h against about 110 K. Exchange gas to linked masses, added in legacy commit `fa99b40`, accounts for part
of that difference; the rest is not recoverable from the committed configuration.
