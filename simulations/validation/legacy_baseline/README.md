# Legacy baseline evidence

Run `uv run python -m simulations.validation.legacy_baseline.run` from the
repository root. The script uses the legacy tree at the commit in `inputs.toml`,
from the sibling `../pulse_tube_cooldown_simulation` repository, regardless of
its checkout. Each run uses its measured cooldown's comparison window and
channel mapping from `data/<system>/cooldown.toml`.
Use `--legacy-repository <path>` to select the source repository and
`--output <directory>` to derive a separate copy of the evidence.

| System | Retained baseline | Pooled RMS log-temperature residual |
|---|---|---|
| SN09 | `sn09_tuned` | 0.459 |
| FP1 | `fp1_generous` | 0.404 |

`inputs.toml` identifies the source legacy commit, saved system, and run changes.
Each baseline has a `.nc` temperature history and a `.json` record of its complete
legacy configuration. `scores.toml` records pooled and per-node residuals,
initial temperatures, and run durations. `sn09.png` and `fp1.png` plot only the
selected baseline over its measured record. `run.py` derives these artifacts.

Every legacy target is 0 K, so each run ends at its time limit and the DHX is
free to cool below its legacy target. All nodes start at the mean of the scored
nodes' first valid readings, which also supplies the legacy ambient temperature.
Each time limit extends one hour beyond the measured comparison window.

SN09's retained run is `Fridge` as saved. Its CFP and STP masses are about half
the BoM masses, and its cooler is the one-compressor estimate although SN09 runs
two compressors. These differences offset one another in its cooldown fit.
FP1's retained run uses its BoM masses and two PT1.5 heads, with the DHX He-4 gas
link and untabulated gas-gap options recorded in `inputs.toml`. Its RGP holds
near 75 K from about 18 h while the measured RGP falls to 8 K.

The earlier candidate comparisons, additional plots, and original derivation
script are preserved with their inputs in the local research snapshot
[`archive/e40e98a-directory-audit/`](../../../../pulse-tube-research/archive/e40e98a-directory-audit/).
