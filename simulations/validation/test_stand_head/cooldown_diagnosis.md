# Stand cooldown diagnosis, 2026-11-07

Why the fitted test-stand template's cooldown misses the stand cooldowns, for task 3.10. Every run used
`data/test_stand/fridge.toml` at commit 32a7992, with the PT1.5 head of 3.8. Variants changed one input in a scratch
copy: the head's `stage_coupling` or `node_count` in `pt1_5.toml`, the drive's `relief_swing`, or a1 and R1 in the
parameter file. The scratch scripts were not kept. Heat into a node is positive throughout.

## Summary

The cooldown runs to the end without a solver failure and reaches the measured base. T2 follows the measured spread.
T1 cools about twice as fast as measured, undershoots to 16 K, and recovers to its base.

At base the fit balances stage 1's acoustic power, E1 = 62 W, against 69 W of a1 regeneration loss. The load maps
constrain their difference and its slope. The regeneration loss scales as `(m c_p)^2 / (f C'_m) dT/dx` and nearly
vanishes while the head is warm, where the swing sits at relief and E1 = 117 W. Stage 1 then has about 100 W of net
capacity from 295 K down to 100 K.

Lowering E1 by about 1.5x, with a1 lowered to hold the base, brings T1 into the measured band on one and two
compressors. It raises stage 1's load-map residual from 0.7 to 3.0 W rms to 3.4 to 5.6 W, with too shallow a slope.
Under ADR 03 as written, the load maps and the cooldowns disagree on stage 1.

## Measured cooldowns

T1 / T2 in K after the compressor starts, from `data/test_stand/cooldown_set.csv`. T1 reads NaN above 219 K.

| Serial | Compressors | 30 min | 45 min | 60 min | 75 min | 90 min | 120 min | T2 below 4 K |
|---|---|---|---|---|---|---|---|---|
| S0007 | 1 | 164 / 133 | 120 / NaN | 85 / 22 | 64 / 7.9 | 51 / 4.0 | 44 / 2.3 | 90 min |
| S0009 | 1 | 205 / 112 | 188 / 52 | 171 / 28 | 156 / 21 | 140 / 17 | 105 / 9.7 | 149 min |
| S0019 | 1 | 190 / 103 | 161 / 41 | 128 / 19 | 98 / 9.3 | 76 / 4.2 | 52 / 2.6 | 92 min |
| S0028 | 2 | 165 / 126 | 124 / 60 | 90 / 18 | 64 / 5.0 | 50 / 2.8 | 38 / 2.3 | 78 min |
| S0035 | 1 | 169 / 140 | 138 / 85 | 108 / 38 | 82 / 16 | 62 / 7.2 | 48 / 2.9 | 105 min |
| S0038 | 1 | 176 / 152 | 148 / 101 | 122 / 56 | 96 / 25 | 69 / 12 | 49 / 3.3 | 114 min |

## The fitted template's cooldown

On one compressor the run takes 199 s and ends at T1 = 38.1 K and T2 = 2.42 K. The regenerator matrix lags its
stages. At 40 min stage 1 is at 17 K, its neighbouring Pb node at 99 K, and the stainless node at 272 K; stage 2 is at
98 K while the Er3Ni and HoCu2 nodes are at 256 and 283 K.

Stage terms along the run in W. Breakthrough loss is omitted; it is inactive at the fitted parameters (`fit.txt`).

| min | T1 | T2 | dP (bar) | E1 | E2 | F1- | -F2+ | coupling | Q1 | F2- | Q2 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 295 | 295 | 19.30 | 116.7 | 18.4 | 2.0 | -0.3 | 0.8 | -114.1 | 0.23 | -18.1 |
| 15 | 203 | 225 | 19.30 | 116.7 | 18.4 | 12.5 | 2.4 | 1.1 | -100.2 | 2.80 | -15.5 |
| 30 | 106 | 153 | 19.30 | 116.7 | 18.4 | 24.7 | 3.9 | 1.6 | -85.6 | 3.99 | -14.2 |
| 39 | 18 | 104 | 19.07 | 113.9 | 17.9 | 102.0 | 7.8 | 2.1 | -0.8 | 4.96 | -12.8 |
| 51 | 20 | 26 | 17.88 | 100.1 | 15.8 | 79.5 | 5.2 | 14.2 | -0.2 | 15.29 | -0.2 |
| 75 | 31 | 13 | 16.34 | 83.6 | 13.2 | 72.8 | 1.6 | 8.9 | +0.8 | 12.87 | 0.0 |
| 90 | 38 | 2.4 | 14.11 | 62.4 | 9.8 | 69.3 | -9.5 | 0.9 | -0.6 | 9.55 | 0.0 |

From about 39 to 75 min stage 1 sits at 16 to 31 K, where F1- balances E1. While stage 2's section is warm, its gas
swings little mass, so the stage-1 boundary's mass flow, and with it the regeneration loss, stays small until T1 is
far below its base. The coupling term adds 9 to 17 W to stage 1 at T2 of 13 to 26 K, outside the load maps' range.

## Capacity at the measured temperatures

S0007's stages held at its measured T1 and T2, each a straight-line fit over 5 min, with the matrix solved to steady
state and continued from the previous time. The plate's net cooling is its heat capacity times the fitted slope.

| min | T1 | T2 | Stage-1 capacity | Plate 1 | Stage-2 capacity | Plate 2 |
|---|---|---|---|---|---|---|
| 15 | 218.0 | 200.8 | -109.6 | -64.4 | -17.8 | -15.7 |
| 30 | 164.8 | 132.1 | -105.1 | -48.3 | -17.5 | -11.0 |
| 45 | 119.9 | 72.5 | -97.0 | -34.0 | -16.9 | -5.6 |
| 60 | 85.3 | 22.3 | -61.3 | -17.2 | -11.9 | -0.2 |

The matrix's own cooling along the steady profiles is about 10 to 20 W. It does not close stage 1's gap of 40 to
60 W.

## Variants

T1 / T2 in K on one compressor unless noted.

| Case | 30 min | 45 min | 60 min | 75 min | Base |
|---|---|---|---|---|---|
| As fitted | 106 / 153 | 16 / 65 | 24 / 18 | 31 / 13 | 38.1 / 2.42 |
| c = 0, not refitted | 103 / 147 | 16 / 47 | 22 / 15 | 28 / 12 | 37.4 / 2.46 |
| 4 nodes per layer | 112 / 152 | 18 / 65 | 23 / 22 | 32 / 15 | 53.9 / 9.56 |
| 8 nodes per layer | 112 / 153 | 18 / 66 | 23 / 22 | 31 / 15 | 54.0 / 9.57 at 120 min |
| Relief swing 15 bar | 186 / 213 | 130 / 171 | 55 / 127 | 16 / 73 | 38.1 / 2.43 |
| R1 x1.5, a1 = 0.065 | 176 / 153 | 115 / 60 | 73 / 17 | 46 / 12 | 36.1 / 2.43 |
| Two compressors, as fitted | 106 / 153 | 16 / 60 | 25 / 16 | 39 / 2.19 | 37.4 / 2.20 |
| Two compressors, R1 x1.5, a1 = 0.065 | 176 / 153 | 115 / 60 | 73 / 17 | 49 / 10 | 35.0 / 2.20 |

- The coupling term and the grid leave the early trajectory unchanged. With 4 or 8 nodes per layer the matrix settles
  on its warm branch, with stage 2 near 9.5 K.
- A 15 bar relief swing slows both stages, and T2 falls far behind the measurements. FP1's gauges put the swing at
  relief while the head is warm (`data/compressors.md`).
- At R1 x1.5 and a1 = 0.065, E1 is 78 W warm and 41 W at base, against 48 W of F1-. Stage 1's net capacity is 72 W at
  274 K and 64 W at 195 K. T1 ends 4 to 8 K below the one-compressor heads' base, and on two compressors T2 lags late.
- At R1 x2 and a1 = 0.05, T1 is 179 K at 40 min and 63 K at 120 min, and stage 2 holds near 9 K until 135 min.

## Load maps at R1 x1.5

The single-compressor points of `data/test_stand/fit_set.csv`, solved as 3.8 solves them, at a1 = 0.065,
R1 = 1.2e10 Pa s/m^3, and the fitted a2 and R2, with c re-solved in closed form: c = 0.29. S0007's 0 W, 0 W point has
no matrix steady state, as in the fit. RMS residual in W, against `fit.txt`:

| Serial | Stage 1, fitted | Stage 1, R1 x1.5 | Stage 2, fitted | Stage 2, R1 x1.5 |
|---|---|---|---|---|
| S0007 | 1.31 | 3.95 | 0.041 | 0.041 |
| S0008 | 1.05 | 3.78 | 0.116 | 0.114 |
| S0009 | 0.93 | 3.43 | 0.151 | 0.148 |
| S0019 | 0.69 | 3.44 | 0.087 | 0.085 |
| S0023 | 2.08 | 4.21 | 0.069 | 0.072 |
| S0034 | 1.67 | 4.32 | 0.154 | 0.159 |
| S0035 | 2.96 | 5.58 | 0.082 | 0.085 |

Stage 1's residual is +1.0 to +7.5 W at the 20 and 25 W heater points and -2.5 to -7.2 W at 0 W. On S0009 at 0.5 W
on stage 2, the model's dQ1/dT1 is about -0.9 W/K against -1.3 W/K measured. This is one point, not a refit. The
couplings against requirement 4 were not checked there.

## Other observations

S0023's Münster trace (`docs/pt_performance_overview_2026-10-30.pdf`, p. 33) logs only the return pressure: about
7 bar warm, a peak of 8.7 bar near 65 min, and 6.5 to 7 bar at base. In the model the return pressure is 7.35 bar
warm, about 7.3 bar mid-cooldown, and 7.67 bar at base, and the swing is smallest at base.

## Open

- In ADR 03 the only stage-1 term steep in T1 at base is the regeneration loss. A stage-1 loss that is steep at base
  and vanishes while warm would let E1 fall. The breakthrough loss has that form, since it scales with
  h(T0) - h(T1), but it is inactive: the stroke from the reservoir's mass swing is about 0.05 against a buffer of
  0.5. Waking it on stage 2 in 3.1 gave dL2/dT2 the wrong sign; stage 1 is untested.
- Whether stand cooldowns may constrain the canonical fit. ADR 03 requirement 2 prefers them for validation only.
