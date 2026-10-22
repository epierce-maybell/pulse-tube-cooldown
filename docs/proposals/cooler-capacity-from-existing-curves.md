# Proposal: cooler capacity from existing curves

Status: draft, 2026-10-22

## Summary

The cooler link removes heat at three points on the cold head: the first stage, the second stage, and the regenerator clamp. In general each capacity depends on all three temperatures: $Q_1(T_1, T_2, T_R)$, $Q_2(T_1, T_2, T_R)$, and $Q_R(T_1, T_2, T_R)$. This proposal builds the first cooler model from data already in hand. Each capacity depends only on its own temperature, with curves derived from a measured test-stand cooldown and checked against the base-temperature load map. A second, deliberately coupled run bounds how much the missing coupling can change predicted cooldown time. Together, the two runs decide whether coupled-capacity measurements are worth their cost, and where in temperature space to take them.

## Problem

The simulated window runs from 300 K until the MXP turns on at 40–50 K. At the end of the window the first stage is at 40–50 K and the CFP is at 15–25 K. Copper's heat capacity is concentrated above about 100 K, so capacity at the warm end sets most of the predicted cooldown time.

Measured capacity data covers two regions:

- **The base corner.** The manufacturer load map covers $T_1 \approx 41$–$67$ K and $T_2 \approx 2.2$–$4.7$ K, with $Q_1 = 0$–$30$ W and $Q_2 = 0$–$2$ W. Test-stand load maps at base fall inside the same region.
- **One cooldown path.** A test-stand cooldown of a bare cold head records both stage temperatures from room temperature to base, sampled every 2 s.

Within the base corner the coupling is one-sided:

- $\partial Q_2 / \partial T_1 \approx 0$. Lines of constant $Q_2$ are flat to within 0.04 K across the map.
- $\partial Q_1 / \partial T_2 \approx -1$ to $-4$ W/K. Raising $T_2$ from 2.2 K to 4.7 K costs 3–9 W of first-stage capacity.

The bare-head path and the fridge path diverge. On the bare head, $T_2 \approx 23$ K when $T_1 \approx 140$ K, and $T_2 \approx 2.5$ K when $T_1 \approx 50$ K. The fridge ends its window with $T_2$ at 15–25 K and $T_1$ at 40–50 K. That end state lies off the bare path and outside the load map. Extending the base-corner coupling linearly to it gives a first-stage correction larger than the first-stage capacity itself. The coupling away from the base corner is therefore the main open quantity.

The brownfield single-variable curves come from other coolers. Near room temperature they span 90–330 W on the first stage and 50–150 W on the second.

## Proposal

### Capacity along the measured path

At each cooled point, the energy balance on the test stand gives

$$
Q_i(t) = -\,C_i\big(T_i(t)\big)\,\frac{dT_i}{dt} + Q_{\mathrm{par},i}(t),
\qquad
C_i(T) = \sum_k m_{ik}\, c_k(T),
$$

where $m_{ik}$ is the mass of material $k$ at point $i$ and $c_k$ comes from the materials library. Pairing $Q_i(t)$ with $T_i(t)$ gives a path curve $Q_{i,\mathrm{path}}(T_i)$.

Processing steps:

1. Drop invalid samples (`nan` and zero readings).
2. Correct thermometer offsets using the readings before cooldown starts, when every sensor sits at room temperature.
3. Differentiate with a smoothing window of about 2 minutes.
4. Estimate parasitic loads from radiation and supports. They are small next to warm-end capacity and significant below about 80 K.

Near the end of the cooldown the path enters the load-map region. There, the path-derived capacity and the load map should agree to within the uncertainty in heat capacity. This is the main consistency check.

### Uncoupled run

$$
Q_1 \approx Q_{1,\mathrm{path}}(T_1), \qquad
Q_2 \approx Q_{2,\mathrm{path}}(T_2), \qquad
Q_R \approx Q_{R,\mathrm{placeholder}}(T_R).
$$

The load map supports dropping $T_1$ from $Q_2$. This model has the same structure as the brownfield curves, but measured on the production cooler rather than taken from literature on other coolers. The curves go into the cooler library as table-form data files with a recorded source. The cooler link still evaluates all three temperatures; the uncoupled model simply ignores the ones it does not use.

### Coupling bracket

A second run lowers first-stage capacity wherever $T_2$ sits above the bare path, matching the sign seen at base:

$$
Q_1 = Q_{1,\mathrm{path}}(T_1)\,\big(1 - f\,s\big),
\qquad
s = \min\!\left(1,\ \frac{\max\big(0,\ T_2 - T_{2,\mathrm{path}}(T_1)\big)}{\Delta T}\right),
$$

with $f$ = 0.2–0.3 and $\Delta T \approx 10$ K as starting values. Both runs are scored against a measured fridge cooldown using the design spec's log-temperature residual, and are compared on time to end temperature.

- If the two runs differ by less than the error budget, the uncoupled model stands for v1.
- If they differ by more, coupled measurements become the priority, aimed at the region of $(T_1, T_2)$ where the runs diverge.

## Data needed

1. **Masses and materials of stage 1 and stage 2 on the test stand, and of whatever is bolted to them.** Each stage includes its flange and its internal cold-end heat exchanger. These set $C_i$. Without them the cooldown yields only $dT/dt$, not capacity.
2. **A measured fridge cooldown in the repository.** Both runs are scored against it, and the decision between them rests on that comparison.

## Extension

Each step below adds data to the same cooler link without changing its interface. Each step is justified by what the coupling bracket shows.

1. **Copper-block cooldowns.** Bare-head cooldowns with a known copper mass on one stage at a time. Each run follows a different path, moving toward the fridge path. The block dominates $C_i$, which reduces sensitivity to the cold head's own mass inventory, and the slower cooldown is closer to how the fridge behaves. The analysis is the same as for the uncoupled run. This method suits the warm end, where heater holds need hundreds of watts on the first stage.
2. **Heater-held steady states below about 150 K.** Heaters under PID control on all three cooled points hold a chosen $(T_1, T_2, T_R)$. At steady state each capacity equals heater power plus parasitic load, with no heat-capacity model involved. Clamp-load steps at each point test a separable regenerator form,
   $$
   Q_R \approx G_R(T_1, T_2)\,\big(T_R - T_R^0(T_1, T_2)\big),
   $$
   where $T_R^0$ is the clamp temperature at zero load. That form reduces the grid from three dimensions to a $(T_1, T_2)$ grid with a few clamp levels per point. Points concentrate where the coupling bracket diverged. Each point settles on a drift threshold for temperature or heater power; a fixed dwell time is insufficient.
3. **Coupled map representation.** The cooler file grows from one-dimensional tables into scattered $(T_1, T_2, T_R) \to (Q_1, Q_2, Q_R)$ data with multidimensional interpolation. The uncoupled curves remain as the zero-coupling case and as a regression baseline.
4. **Quasi-static check.** A map built from steady states predicts a bare-head cooldown. Agreement confirms that the cold head's internal state (chiefly the regenerator temperature profile) follows the boundary temperatures closely enough for a static map to apply during cooldown.

This sequence spends measurement effort only where the simulation shows sensitivity, and every intermediate model is usable and testable against fridge data.
