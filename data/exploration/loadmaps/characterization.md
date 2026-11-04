# Load-map characterization

What the existing load maps say about how stage heater powers and stage temperatures relate, before any functional form is chosen for the cooler. It takes these positions: Denver and Münster load maps are the cooler's fit inputs; Denver post-tuning cooldowns are its test; Münster cooldowns, taken before tuning, are unused; every run is an independent dataset; and the clamp is out of scope, since no load map has a clamp heater or thermometer.

## Data

There are 36 datasets from 21 heads, catalogued with provenance in `runs.csv`:

- 16 Denver maps, from the stand archive in `../../data/test_stand/loadmap_data/` relative to the repository.
- 19 Münster maps: S0007's numeric file, plus 18 maps digitized from the figures in `docs/pt_performance_overview_2026-10-30.pdf`.
- S0028's 11-point Münster stage-2 sweep at 20 W.

Excluded files are duplicates or retyped Münster data: the second SN009 May 14 and SN144 Apr 23 files, the Mar 27 four-point subset, and the SN008 Jan 8 and SN023 "Germany" files, which carry Münster key points in the Denver format.

Digitizing reproduces S0007's numeric Münster map within 0.02 K on T1 and 0.01 K on T2, and every printed key point within 0.08 K and 0.02 K, except two printing errors in the overview: S0048 at 25 W, 0.5 W prints T2 = 2.59 K where its figure shows 3.59 K, and S0023 at 20 W, 0 W prints 2.4 K where its figure shows 2.48 K (figure 10).

The maps span T1 from 33 to 70 K, T2 from 2.1 to 5.2 K, Q1 from 0 to 35 W, and Q2 from 0 to 2.5 W. Q is heater power, so each stand's parasitic heat into a stage is folded into it.

## Findings

**1. Maps of one head taken weeks apart differ about as much as maps of different heads, but the records cannot say why.** Median absolute differences in K between pairs of datasets (figure 03, `summary.txt`):

| Pair of datasets | T1 at 0 W, 0 W | T1 at 20 W, 0 W | T2 at 0 W, 0 W | T2 at 0 W, 0.5 W |
|---|---|---|---|---|
| Repeated point within a run | 0.06 | 0.06 | 0.03 | 0.03 |
| Same head and stand, 1 day apart | 0.03 | 0.31 | 0.01 | 0.01 |
| Same head and stand, 8 to 103 days apart | 1.6 | 2.9 | 0.02 | 0.23 |
| Same head, other stand | 4.2 | 7.7 | 0.10 | 0.17 |
| Different heads, Denver | 1.5 | 2.7 | 0.14 | 0.22 |
| Different heads, Münster | 2.5 | 3.5 | 0.14 | 0.13 |

The within-run row pools the four repeated points in the archive, whatever their powers. The weeks-apart pairs are S0007 from March to May across a hose change, S0009 from May to August with August's compressor count unrecorded, and SN144 eight days apart. Different heads are regular 1.5 heads on one compressor.

Same-day repeats bound measurement noise. The weeks-apart differences are not evidence of drift in the head itself. Two of those pairs span known drive changes. Only S0007's March maps record valve settings, and later Denver maps log none, so retuning between runs cannot be ruled out.

**2. Stage-2 power depends on T2 alone.** At fixed T2, the change of Q2 with T1 has a median of 0.000 W/K and quartiles within 0.002 W/K on one compressor, and about −0.005 W/K on two. Over the mapped 40 to 60 K of T1 that is under 0.05 W (figure 07). Seen as measured, 20 W on stage 1 raises T2 by up to 0.15 K with stage 2 unloaded, and lowers it by about 0.1 K with stage 2 loaded (figure 06).

**3. Every head shares one stage-2 curve for a given drive; heads differ by a constant.** The slope of Q2 against T2 rises with T2 and has the same value on every head at a given T2 (figure 07):

| T2 | 2.5–3 K | 3–3.5 K | 3.5–4 K | 4–4.5 K | 4.5–5 K |
|---|---|---|---|---|---|
| dQ2/dT2, one compressor (W/K) | 0.41 | 0.46 | 0.66 | 0.71 | 0.96 |
| dQ2/dT2, two compressors (W/K) | 0.86 | 1.04 | 1.20 | 1.44 | |

A curve shared by all 22 regular 1.5 and HX1 single-compressor runs, plus one constant per run, fits every stage-2 point with a residual of 0.059 W rms and 0.18 W at most. The constants span −0.03 to 0.25 W, and they set each run's unloaded T2. A slope per run lowers the residual to 0.038 W. The two-compressor runs share a steeper curve, with a residual of 0.096 W rms (figure 11). As measured, T2 rises 1.8 to 3.4 K/W for the first few tenths of a watt and about 1 K/W above 1 W; two-compressor maps run at roughly half those slopes (figure 05).

**4. Stage-1 power depends on both temperatures, and its shape varies between heads.** On one compressor, dQ1/dT1 at fixed T2 falls from 1.55 W/K at 40 to 45 K to 1.24 W/K at 55 to 60 K. dQ1/dT2 at fixed T1 runs from −1.0 W/K at 2.5 to 3 K to −1.45 W/K at 4.5 to 5 K. On two compressors it runs −2.1 to −3.2 W/K. Each watt on stage 2 costs 1.7 to 2.6 W of stage-1 capacity at fixed T1 on one compressor, and 2.3 to 3.1 W on two (`summary.txt`).

Unlike stage 2, no shared shape collapses stage 1. A shared curve with a constant per run leaves 1.6 W rms. A slope per run plus a shared dependence on T2 brings that to 0.71 W rms, with a 2.8 W maximum, over 0 to 25 W. T1 rises 13 to 19.5 K for 20 W across single-compressor regular 1.5 runs, and its slope steepens with load on most heads (figure 04).

**5. On some heads the first fraction of a watt on stage 2 lifts an unloaded stage 1 by 3 to 4 K.** With stage 1 unloaded, 0.5 W on stage 2 raises T1 by 3.1 to 4.1 K on S0010, S0022, S0023, and S0034, against 0.5 to 1.8 K on every other head. SN121 shows 3.0 K. The step belongs to the head, not the stand: S0023 and S0034 show it in Münster and again in Denver. With 20 or 25 W on stage 1, the same 0.5 W raises T1 by only 0.5 to 1.1 K on those heads. Their stage-2 response is ordinary. No single shipped valve setting orders this step or any other map feature (figure 09).

**6. A second compressor changes both stages.** S0006 and S0028 on both stands rise 10.1 to 11.7 K in T1 for 20 W on stage 1, below every single-compressor map. Their stage-2 slopes are about twice as large. Their unloaded T2, 2.3 to 2.6 K, is no lower than on one compressor. In S0009's August map, compressor count unrecorded, T2 rises 1.00 K for 0.5 W. That sits inside the single-compressor range of 0.84 to 1.52 K and above every two-compressor map's 0.56 to 0.88 K.

**7. T3 has two regimes.** T3, on the pulse tube below stage 1, sits at 33 to 35 K in most Denver runs whatever the load. In five single-compressor runs, and in S0009's August run, it climbs once stage 1 carries 20 or 25 W, reaching 58 K, close to T1. The five are SN121, both SN144 maps, S0008 in May, and the unnamed Sep 21 head. On S0028 with two compressors it sits higher, at 37 to 41 K, and rises steadily with T1. SN144's April 23 map is the only one where stage-2 load lowered T1 with stage 1 loaded: 0.5 W at 20 W dropped T1 by 1.9 K and T3 by 4.0 K at the same time (figure 08).

**8. Two Denver maps carry irregular points.** In S0007's March 27 map, stage-2 points at 0.35 and 0.5 W fall below their neighbors' trend, and about 0.25 K below the April maps at 0.5 W. At 10 W, its 0.35 W point reads 2.94 K. S0007's April 13 map logs its last point at 0.5 W on stage 2, which is read as 0 W from its temperature and the map's order.

## Functional forms

Figures 12 to 15 read every grid as Q1(T1, T2) and Q2(T1, T2). Slopes are local partials from grid cells. Cuts at a fixed temperature interpolate linearly between measured points. The fits are in `functional_forms.txt`. Everything here holds over T1 from 35 to 70 K and T2 from 2 to 5 K only.

**Q2 depends on T2 alone, and its slope grows as a power of T2.** After each run's constant is removed, points at one T2 do not sort by T1 (figure 12, first two panels). Cuts at fixed T2 stay within about 0.1 W across 20 K of T1. Above 4 K they stray to 0.3 W, where the measured stage-2 powers are 0.5 W apart and interpolation error dominates. Pooled over runs, dQ2/dT2 goes as T2 to the 1.38 power on one compressor and 1.34 on two, and two compressors double the prefactor:

```
dQ2/dT2 = 0.100 W/K * (T2 / 1 K)^1.38    one compressor, 115 cells
dQ2/dT2 = 0.196 W/K * (T2 / 1 K)^1.34    two compressors, 35 cells
```

Single runs give exponents from 0.9 to 2.1, most between 1.2 and 1.85, with nearly parallel lines in log-log (figure 12, last panel). Integrated, Q2 rises as about T2^2.4 above a constant for each run.

**At fixed T2, Q1 is close to linear in T1, with a slope that falls as T1 rises.** Lines of fixed Q2 keep their shape and shift to higher T1 as T2 rises (figure 13). Pooled, dQ1/dT1 goes as 1/T1, from about 1.6 W/K at 42 K to 1.1 W/K at 62 K on one compressor, and about 10% higher on two. Single runs give exponents from −0.3 to −2.0, so heads differ more in level than the pooled trend suggests (figure 15, first panel).

**At fixed T1, Q1 falls nearly linearly with T2, and the coupling weakens as T1 rises.** On most heads, a cut at fixed T1 is close to a straight line in T2 (figure 14). The step heads of finding 5 bend instead: most of their drop comes between the unloaded point and the first stage-2 power, at the coldest T1. Pooled over cells:

```
dQ1/dT2 = -1.40 W/K + 0.050 W/K^2 (T1 - 50 K) - 0.31 W/K^2 (T2 - 3.5 K)    one compressor, rms 0.65 W/K
dQ1/dT2 = -2.73 W/K + 0.020 W/K^2 (T1 - 50 K) - 0.49 W/K^2 (T2 - 3.5 K)    two compressors, rms 0.71 W/K
```

All 11 runs that span enough T1 show the coupling weakening as T1 rises, by 0.01 to 0.11 W/K per kelvin (figure 15, last panel). The scatter between heads is larger than either trend. The lines of fixed Q1 in the T1-T2 plane tilt by dT1/dT2 = 0.92 at the median on one compressor, with quartiles 0.65 to 1.22, and by 1.79 on two.

## What the maps cannot say

- **Nothing far from base.** Nothing is mapped above 70 K on T1 or 5.2 K on T2, where most of a cooldown happens.
- **Nothing about the clamp.** No map loads or measures the regenerator clamp.
- **Constants mix head and stand.** Heater power is not net cooling, so each run's constant contains its stand's parasitic heat.
- **Steady state is unverified.** Each point dwells 90 to 120 min on the stand. The continuous logs that cover each map could check it.
- **Most compressor counts are assumed.** Denver maps are taken as one compressor, from `data/compressors.md`, except S0028 and S0009's August map.

## Files

- `figures/16_capacity_ratio.png`: `(dQ1/dT2) / (dQ2/dT2)` with both partials at fixed T1, plotted against each cell's mean T1 and T2. This is `dQ1/dQ2 at fixed T1`, in W/W, from the same grid-cell derivatives as figure 07. The numerical summaries are labeled `dq1_dq2` in `summary.txt`.
- `runs.csv`: one row per dataset, with provenance and caveats. A source prefixed `archive/` is a file in the stand archive.
- `muenster_key_points.csv`: the overview's printed key points, with page numbers.
- `valve_settings.csv`: the overview's shipped valve settings, with page numbers.
- `extract_figures.py`: extracts the overview's figures into `muenster_figures/`, which is not tracked.
- `digitize.py` and `digitize.toml`: digitize those figures into `muenster_maps.csv`, with overlays in `digitized/`.
- `runs.py`: loads every dataset.
- `characterize.py`: draws figures 01 to 11 and 16 and writes `summary.txt`.
- `functional_forms.py`: draws figures 12 to 15 and writes `functional_forms.txt`.

Rebuild everything with `uv run --with pypdf python extract_figures.py`, then `uv run python digitize.py`, `uv run python characterize.py`, and `uv run python functional_forms.py`, all run from this folder.
