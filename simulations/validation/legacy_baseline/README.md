# Legacy baseline

The legacy model run on SN09 and FP1 and scored with `measurements.residual`. `../legacy_baseline.py` makes every
file here except this one and `inputs.toml`. Run
`uv run python -m simulations.validation.legacy_baseline` from the repository root.

## Baselines

Each system's baseline is its best legacy score over committed legacy inputs: any committed legacy version, cooler
data, and mix of committed saved-system pieces, with no parameter fit to the measured cooldown and no change to legacy
code. `inputs.toml` names each system's baseline run.

| System | Baseline run | Pooled RMS |
|---|---|---|
| SN09 | `sn09_tuned` | 0.459 |
| FP1 | `fp1_generous` | 0.404 |

FP1's candidates were the 16 combinations of four changes in legacy commit `ddc79e8`: its `PT1.5 (1X F-70)`
regenerator curve, and from its saved system "FP1 - 2x PT1.5 (imported)" the RGP's added 10 kg of stainless steel and
6 kg of lead, the DHX He-4 gas link, and gas gaps without tabulated conductance. Every combination with the
regenerator curve scores between 0.574 and 0.904, and every one without it between 0.404 and 0.539; the saved system
as committed scores 0.574. SN09's candidates were the DHX gas link and untabulated gas gaps, each alone and together, on
`sn09_tuned` and `sn09_corrected`. With either change, each scores between 0.506 and 0.785, except `sn09_tuned` with
both, which the legacy engine ends at 57.0 h, before the window ends.

## Runs

| Run | Configuration | Pooled RMS |
|---|---|---|
| `sn09_tuned` | `Fridge` as saved: its masses and `PT1.1 (1X F-70)`, the configuration legacy was tuned to SN09 with | 0.459 |
| `sn09_corrected` | `Fridge` with SN09's BOM masses and `PT1.1`, two F-70s | 0.806 |
| `fp1` | `Fridge` with FP1's BOM masses and two `PT1.5 (1X F-70)` | 0.495 |
| `fp1_generous` | `fp1` with the DHX He-4 gas link of `ddc79e8`'s FP1 saved system and gas gaps without tabulated conductance | 0.404 |

Per-node scores are in `scores.toml`; `sn09.png` and `fp1.png` plot each system's runs over its measured cooldown. `sn09_tuned` matches SN09 through two offsetting differences from the system:
`Fridge`'s CFP and STP masses are about half the BOM's, and its cooler is the one-compressor estimate though SN09 runs
two. With `Fridge`'s masses and two compressors, SN09 scores 0.660; with the BOM masses and two, 0.806. FP1 with
`Fridge`'s masses scores 0.918. In `fp1_generous`, the RGP holds near 75 K from about 18 h while the measured RGP
falls to 8 K.

## Legacy version

`Advanced_Feature_SIM` at `inputs.toml`'s `legacy_commit`, the only legacy line with the six-plate model.
`legacy_baseline.py` runs the legacy tree at that commit, whatever is checked out. The legacy
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
