# Regenerator-2 loss in the stage-2 warm response

## Change

The stage-2 warm correction gains a loss term:

```text
s_2 = A_2 * (u2 / T0)^b_2 - G_R2 * (u1 - u2)
```

`u1` and `u2` are the clamp-shifted stage temperatures in K. `G_R2` is the regenerator-2 loss conductance in W/K,
standing for `mdot * c_p * (1 - epsilon)` of regenerator 2. It is a shared parameter, set to 0.2 W/K and not fitted.
The term enters through the existing warm blend `w`, so the cold map and its derivatives are unchanged.

Stand fits hold the stage-1 warm amplitude at 100 W, because the stand does not determine it.

## Evidence

Research is in `../pulse-tube-research/simulations/exploration/fridge_head_budget/`. While the CFP is no warmer than
the PRP and RGP, both refrigerators' stage-2 heads under ADR-06 are 2-3x short of the stage-2 side's measured heat
loss. The S0028 stand rejects a larger stage-2 amplitude. With the loss term, the stand fits hold and the transient
RMS improves:

| Run | RMS before | RMS with loss term |
|---|---|---|
| FP01 baseline | 0.944 | 0.665 |
| SN09 baseline | 1.661 | 0.752 |

`structural.py` holds a reference implementation (`RegeneratorLossColdMap`). `structural.json`
(`kappa_0.2_stage_1_100W`) and `sn09_check.json` hold the refit coefficients and the stand errors.

## Work

- Add `regenerator_2_loss_conductance` to `shared_parameters.toml`, `SharedParameters`, and `PulseTubeHead`, and
  apply it in `calc_energy_balance_terms`.
- Add held coefficients to `fit_parameters.toml` and `FitParameters`, honored by `fit_pulse_tube_head`.
- Regenerate `archetype.toml` and `data/test_stand/s0028/calibration.toml`.
- Rerun the four transient runs.
- The archetype refit in research missed the stage-2 3 K crossing by 12%. Check it against the acceptance test.
- ADR-06 needs a matching amendment.
