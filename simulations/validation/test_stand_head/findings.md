# Stage coupling in the ADR 03 head: findings, 2026-11-06

Evidence for `TASKS.md` section 3 and the proposed ADR 03 wording. Heat flow into a node is positive throughout. The
measured values are the one-compressor partials in `simulations/exploration/loadmaps/summary.txt`, converted from
heater power to heat into the node by a change of sign. Model values are for the test-stand template,
`data/test_stand/fridge.toml`, on one compressor at 1 Hz, with unfitted head parameters. The placeholder parameters
are a = (0.15, 0.05) and R = (1e10, 3e10) Pa s/m^3.

## Energy budget

With no clamp load at steady state, the energy flow is uniform through each regenerator section. The flow leaving
stage 1 into the stage-2 section therefore equals the flow arriving at stage 2, and

    Q1 + Q2 = F1(warm flange) + L1 + L2 + K1 (T0 - T1) + K2 (T0 - T2) - E1 - E2
    dQ1/dT2 = -dQ2/dT2 + d(Q1 + Q2)/dT2

at fixed T1. Let r = (dQ1/dT2) / (-dQ2/dT2) be the stage-1 capacity lost per watt of stage-2 capacity.

- Measured, one compressor: r has a median of 1.7 to 2.6 across the T2 bins from 2.5 to 5 K, with interquartile
  ranges between 1.3 and 3.0. Two compressors: medians of 2.3 to 3.1.
- Physics inside the stage-2 section moves heat between the stages one for one: its loss law, the gas-to-matrix
  capacity ratio, the matrix heat capacity, and the discretization. Alone, they give at most r = 1.
- r above 1 needs the total net refrigeration to fall as T2 rises. ADR 03 has three channels for that: the pressure
  swing setting E1 and E2, the stage-1 section's warm-end loss, and the tube losses. All three moved the other way
  at every state checked.

## ADR 03 as written

Model couplings at T1 = 58 K and T2 = 3.4 K, from the 2.10 diagnostics:

| Parameters | dQ1/dT1 | dQ1/dT2 | dQ2/dT1 | dQ2/dT2 | r |
|---|---|---|---|---|---|
| Placeholder | -1.351 | -0.425 | -0.0035 | -0.492 | -0.86 |
| Both R x4, both a / 4 | -0.344 | -0.128 | -0.0017 | -0.132 | -0.97 |
| a1 / 4 | -0.788 | -0.279 | -0.0031 | -0.492 | -0.57 |
| R1 x4, a1 = 0.02 | -0.310 | +0.067 | -0.0023 | -0.491 | +0.14 |
| As above, then R x4 and a / 4 | -0.080 | +0.015 | -0.0016 | -0.132 | +0.11 |

Measured at 3 to 3.5 K: dQ1/dT2 = +1.18 (interquartile +0.76 to +1.54), dQ2/dT2 = -0.46 (-0.42 to -0.53),
dQ2/dT1 = 0.000, and dQ1/dT1 = -1.32 at 50 to 55 K, all in W/K.

At the placeholder state:

- The drive is below relief: the swing is 16.95 bar against a 19.3 bar relief, so a lighter cold end raises the
  swing and every acoustic power.
- Breakthrough loss is dormant. The tube strokes, reservoir mass swing over tube gas mass, are 0.053 and 0.010
  against a buffer fraction of 0.5. A buffer fraction low enough to wake it gives 5 to 36 W of loss, with dL2/dT2 of
  the wrong sign.

Local fits of a1, a2, R1, and R2, at T1 = 50 K and T2 = 3.4 K without voids, never exceeded dQ1/dT2 = +0.11 W/K over
764 evaluations.

## Cold gas inventory

Helium's mass swing over the cycle, at 15.65 bar mean and 16.95 bar swing, in kg/m^3:

| T (K) | 3.0 | 4.0 | 5.0 | 6.0 | 8.0 | 10.0 | 15.0 | 30.0 | 58.0 |
|---|---|---|---|---|---|---|---|---|---|
| Isothermal | 17.9 | 21.4 | 27.1 | 37.1 | 73.5 | 78.4 | 54.1 | 25.0 | 13.1 |
| Adiabatic | 16.2 | 17.1 | 18.4 | 20.4 | 27.5 | 33.5 | 28.0 | 14.9 | 8.1 |

A warming cold end below about 8 K swallows more gas per cycle, which lowers the swing below relief and raises the
stage-1 section's flow at relief. The model places little gas there. It has no heat-exchanger voids, and tube 2's
gas lies on a linear profile from the warm flange.

CAD 910-000040 (`data/heads/pt1_5/geometry.csv`) gives cold heat-exchanger voids of 39.3 cm^3 at stage 1 and
10.8 cm^3 at stage 2. Pulse tube 2 clears the stage-1 flange by 2.0 mm, so it is anchored only at stage 2. T3, on
tube 2 about 51 mm below stage 1, reads 33 to 35 K at base. A linear profile from T2 to T1 predicts about 32 K
there, and ADR 03's profile from the warm flange about 130 K; that difference is unexplained.

Couplings with voids at stage 2 or both stages, at placeholder parameters and T1 = 58 K, T2 = 3.4 K:

| Voids (cm^3) | dQ1/dT1 | dQ1/dT2 | dQ2/dT1 | dQ2/dT2 | Stable |
|---|---|---|---|---|---|
| none | -1.35 | -0.43 | -0.004 | -0.49 | yes |
| 10 at stage 2 | -1.26 | +0.39 | +0.045 | +0.04 | no |
| 30 at stage 2 | -1.17 | +0.78 | +0.139 | +0.58 | no |
| 60 at stage 2 | -1.10 | +1.45 | +0.228 | +1.38 | no |
| CAD: 39.3 and 10.8 | -1.39 | +0.48 | +0.045 | +0.07 | no |

The voids turn the total-refrigeration channel the right way, since d(Q1 + Q2)/dT2 goes from -0.92 to +0.55 W/K at
the CAD volumes. But their flow passes through the stage-2 section and flattens stage 2's own slope.

Continuation with CAD voids at T1 = 50 K, T2 = 3.4 K, a = (0.15, 0.05), R1 = 1e10 Pa s/m^3, lowering R2:

| R2 (Pa s/m^3) | Stage-2 heater (W) | dQ2/dT2 | dQ1/dT2 | Growth rate (1/s) |
|---|---|---|---|---|
| 3.0e10 | -1.14 | +0.015 | +0.55 | +0.004 |
| 2.7e10 | -0.38 | -0.18 | +0.37 | -0.003 |
| 2.45e10 | +0.21 | -0.43 | -0.24 | +0.77 |
| 2.2e10 | +0.74 | -0.44 | -1.86 | +6.7 |
| 2.16e10 | | | | no steady state |

Lowering a2 from 0.05 behaves the same way, losing the steady state at a2 = 0.037. Every combination of a1 in
{0.02, 0.15} and R1 in {1e10, 3e10} reached those same limits. On the stable side, dQ1/dT2 falls toward 0 W/K as stage
2's slope approaches the data. Stage 2's slope reaches the data's range only on the unstable side. That instability
resembles the one 2.10's diagnosis traced to the point-sampled Gd2O2S heat capacity in the regeneration term; 3.3
and the 3.6 gate test whether it moves.

## Stage-coupling term

The term was added to stage 1's heat flow in a scratch copy, without voids, at T1 = 50 K and T2 = 3.4 K. Here
x = beta(T2) T2 / (1 + epsilon), with epsilon the coldest stage-2 node's gas-to-matrix capacity ratio.

Untapered, `c E2 x`: dQ1/dT2 rises in a straight line with c, from -0.42 to +1.02 at c = 1.5 for the placeholder.
The other three couplings and the growth rate are unchanged. It adds about 46 W to stage 1 at room temperature and
36 W at T1 = 80 K, T2 = 30 K.

Tapered, `c E2 x (1 - x)`:

| Parameters | c = 0 | c = 1.5 | c = 2.0 | c = 2.5 |
|---|---|---|---|---|
| a = (0.15, 0.05), R = (1e10, 3e10): dQ1/dT2 | -0.42 | +0.86 | +1.29 | +1.71 |
| a = (0.02, 0.05), R = (6.52e10, 2.98e10): dQ1/dT2 | +0.11 | +1.43 | +1.87 | +2.32 |

dQ1/dT1, dQ2/dT1 = +0.01, dQ2/dT2 = -0.49 and -0.50, and the growth rate are unchanged. Its added stage-1 load at
c = 2.0, in W:

| T1, T2 (K) | 295, 295 | 150, 120 | 80, 30 | 60, 12 | 55, 7 | 50, 3.4 |
|---|---|---|---|---|---|---|
| Added load | 0.8 | 1.9 | 11.3 | 12.8 | 5.3 | 2.4 |

The check covered one operating point on one compressor, with unfitted parameters, and stage-2 heater powers of 0.03
to 0.04 W against S0009's 0.61 W.

## Proposed ADR 03 wording

Requirement 4, replacing its first sentence:

> Coupling between stages should follow from model physics where it can. Where the physics as modeled cannot
> reproduce it, one design-level coefficient may carry it, tied to a modeled quantity rather than fitted as a free
> cross-term.

Governing equations, after the paragraph defining the Q1 and Q2 terms:

> Stage 1 also carries $c E_2 x (1 - x)$, where $x = \beta T_2 / (1 + \epsilon)$ is stage 2's gross-cooling fraction,
> taken at the stage-2 boundary with the coldest node's $\epsilon$. The term is the stage-1 cost of producing stage
> 2's cooling in the non-ideal helium regime, which the regenerator and drive model do not capture. $c$ takes one
> value per design, fitted to the population's load maps.

Mass swing, after "Pore gas exchanges heat with the matrix, so its density changes at its node's temperature.":

> So does the gas in each stage's cold heat exchanger, at its stage temperature.

Regenerator network, replacing the sentence that lists the volumes at the stage-1 end:

> At the stage-1 end of the stage-1 section, those volumes are the stage-1 heat exchanger, tube 1, reservoir 1, the
> whole stage-2 section, the stage-2 heat exchanger, tube 2, and reservoir 2.

Charge, in the list of gas volumes: "supply and return lines, buffers, head dead volumes, heat-exchanger voids, pore
gas, pulse tubes, and reservoirs."

Parameters, after the four fitted parameters:

> Each design also has a coupling coefficient $c$, fitted once to the population's load maps. Each design specifies
> the void of each stage's cold heat exchanger among its gas volumes.
