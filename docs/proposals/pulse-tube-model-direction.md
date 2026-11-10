# Pulse-tube model direction

This document records the settled understanding for the pulse-tube model direction. It supports the direction decision for the fridge cooldown simulator.

The direction, thermal structure, and first-candidate approximations below are agreed. The complete constitutive equations and parameter-estimation strategy remain unfinished. This proposal is the current discussion record. [DESIGN](../DESIGN.md#open-decisions) summarizes the architecture and the in-progress head decision.

Session handoff, 2026-11-09: the next deliverable is one consolidated, executable candidate equation set, following the four bounded steps below. No numerical test of the new candidate has been performed. The agreed structure, provisional derivations, and remaining choices are distinguished below so that a fresh session can continue from this point.

## Language

**Stage temperatures (`T1`, `T2`)**: The temperatures at the head's PRP and CFP connections, respectively, in K.

**Clamp temperature (`TR`)**: The temperature at the head's RGP connection, in K. The regenerator clamp holds RGP.

**Head heat flows (`Q1`, `QR`, `Q2`)**: Heat flows from the head into PRP, RGP, and CFP, respectively, in W. Positive values heat the receiving node; negative values cool it. Every heat flow uses the positive-into-the-receiving-node convention.

**Differential-equation approach**: A model understood through dynamical systems and differential equations. Spatial temperature fields are optional.

**Bulk coefficients**: Parameters that describe aggregate behavior without necessarily corresponding to individual physical components or processes.

**Internal thermal state**: Evolving quantities representing thermal storage inside the pulse-tube head, in addition to its stage and clamp temperatures.

**Internal control volume (`a`, `b`)**: A material piece of the regenerator assembly represented by an entry `E_a` in the stored-energy vector, in J. Its energy is calculated from the material heat capacity and the temperature profile over its physical extent. Recovering the profile temperatures from the energy entries can require a coupled algebraic calculation.

**Regenerator matrix**: The solid material that stores and exchanges heat with the gas. The chosen representation has a state vector with seven material-piece energy entries. Use that wording rather than calling the representation "seven states" or "a seven-state representation." Seven energy entries do not imply seven independently fitted parameters.

**External connection (`p`)**: One of PRP, RGP, and CFP. The fridge model supplies its temperature `T_p`; the head returns heat flow `Q_p` into that node. The fridge model integrates the external node's stored energy.

**Exchange flow (`q(p --> a)`)**: Heat flow from source `p` into destination `a`, in W, positive into the destination. For a passive connection, `q(p --> a) = K_ap * (T_p - theta_a)` and `q(a --> p) = -q(p --> a)`, where `K_ap` is the connection's conductance in W/K.

The passive conductance expression illustrates the sign convention; a passive-only regenerator model is not an adopted closure.

**Manifold**: The compressors and pulse-tube heads sharing a gas supply.

## Working formulation

```text
dE/dt = F(T, E, d)
Q = G(T, E, d)
```

`T = (T1, TR, T2)` contains the PRP, RGP, and CFP temperatures in K. `E` contains stored energies in the modeled regenerator sections, in J. `t` is time in seconds. `d` describes the configured drive and shared manifold. `F` gives rates of energy accumulation in W. `Q = (Q1, QR, Q2)` contains heat flows into PRP, RGP, and CFP, in W. `G` calculates those heat flows from the complete state. Boundary conditions, including warm-flange temperature, and fitted parameters are implicit in these functions.

`F` and `G` share the heat flows in the component energy balances. Load maps constrain the steady state, where `F = 0`; test-stand cooldown curves constrain the evolution of internal energy and its effect on cooling power.

Cycle-averaged energy balances are the basis for the transport formulation. Gas enthalpy retains its temperature and pressure dependence across the full cooldown; a linearization appropriate near room temperature does not establish a closure at low temperatures.

## Agreed thermal representation

Integrate stored energies of the regenerator material pieces. Represent the mean temperature profile as piecewise linear, using material-piece center temperatures, boundary temperatures, and the clamp temperature within the regenerator. Calculate stored energy by integrating over each material's actual extent. At each evaluation, determine the profile temperatures algebraically from those energies and the boundary conditions. This inversion does not impose thermal equilibrium between pieces.

The seven material pieces, ordered warm to cold, are:

1. Regenerator 1: stainless-steel screens.
2. Regenerator 1: lead.
3. Regenerator 2: lead.
4. Regenerator 2: Er3Ni.
5. Regenerator 2: HoCu2 above the clamp.
6. Regenerator 2: HoCu2 below the clamp.
7. Regenerator 2: Gd2O2S.

Use the existing head geometry and material data. Include the known regenerator material inventory from the outset. Additional effective thermal storage is an optional extension, not part of the agreed first candidate. The exact booking of gas internal energy within the storage relations remains to be specified.

Use finite thermal contact to the external clamp plate:

```text
QR = K_R * (theta_c - TR)
```

`theta_c` is the internal regenerator temperature at the clamp, in K. `TR` is the external plate temperature, in K. `K_R` is the contact conductance, in W/K. Determine `theta_c` from the adjacent energy flows and this contact relation, without adding a separate clamp energy entry. Unloaded cold load maps alone do not determine `K_R`.

### Energy balances

For each material piece, the rate of stored-energy change equals incoming minus outgoing total energy flow, including any side extraction at the clamp. Summing the pieces gives:

```text
dE_R1/dt = H1h - H1c
dE_R2/dt = H2h - H2c - QR

Q1 = H1c - H2h - P1
Q2 = H2c - P2

Q1 + Q2 + QR = H1h - P1 - P2 - d(E_R1 + E_R2)/dt
```

`E_R1` and `E_R2` are the sums of the material-piece energies in each regenerator, in J. `H1h`, `H1c`, `H2h`, and `H2c` are total energy flows at their warm and cold ends, in W, positive toward the cold end. `P1` and `P2` are net powers removed through the two pulse tubes, including their losses, in W. These balances establish the structure; the laws that calculate every `H` and `P` remain to be completed.

Cold load-map Jacobians must be evaluated after internal energies and the unloaded clamp have relaxed. Holding the regenerator energies fixed would calculate a different response.

## Agreed cycle approximations

- Use acoustic power as the drive measure.
- Average over the individual cycles; evaluate the resulting expressions algebraically as cooldown proceeds.
- Use fast gas-matrix exchange, retaining the small temperature differences needed for transport and the finite matrix heat-capacity effect.
- Hold thermophysical properties at the current local mean temperature and pressure during a cycle calculation; update them during cooldown. This approximation was accepted for computational simplicity without a separate numerical check.
- Use the full represented matrix heat capacity in the local cycle response and in slow energy storage.
- Use a local temperature-gradient approximation rather than tracking gas parcels through the regenerator.
- Use favorable fixed phasing as an effective approximation. A 90-degree pressure-displacement offset corresponds to approximately in-phase pressure and flow; it does not mean pressure and flow are 90 degrees apart everywhere.

Cycle averaging and small-oscillation linearization are distinct approximations. The former is accepted; a global small-oscillation closure is not established. The complete drive calculation must stay within the intended reduced-model scope.

## Profile-dependent power delivery

The evolving regenerator temperature profiles are inputs to power delivery. This additional dependence is accepted as a necessary component of the first candidate.

For regenerator 2, use an effective linear viscous resistance with one aggregate packing scale:

```text
s2 = length_average[(mu(T(x), p) / rho(T(x), p)) / (mu_ref / rho_ref)]
R_flow,2 = R_ref,2 * s2
```

`x` is position along the regenerator, in m. `T(x)` is its current piecewise-linear mean temperature, in K; `p` is mean pressure, in Pa. `mu` is dynamic viscosity, in Pa s, and `rho` is gas density, in kg/m^3. Reference properties are evaluated at the warm reference condition. `s2` is dimensionless and equals one at that uniform reference condition. `R_ref,2` is the fitted reference pressure-drop/mass-flow resistance, in Pa s/kg. Its temperature dependence is calculated from properties, rather than fitting `s2` independently. Material-specific permeabilities and inertial pressure-drop terms are not part of this approximation.

### Two gas paths after regenerator 1

Regenerator 1 is common upstream. At its cold end, the stage-1 junction feeds:

1. Stage 1, pulse tube 1, and its warm flow-control/reservoir connection.
2. Regenerator 2, stage 2, pulse tube 2, and its warm flow-control/reservoir connection.

The clamp is a thermal side contact along regenerator 2, not another gas path.

For effective resistances referred to the stage-1 junction, the working split under the common-pressure, in-phase approximation is:

```text
R_branch,2 = R_other,2 + R_ref,2 * s2
f = R_branch,1 / (R_branch,1 + R_branch,2)
  = 1 / (1 + b + a*s2)

a = R_ref,2 / R_branch,1
b = R_other,2 / R_branch,1

W1 = (1-f) * A
W2h = f * A
```

`R_branch,1` is the effective resistance of the first path; `R_other,2` represents the second path's remaining effective resistance. Both have units Pa s/kg. `a`, `b`, and `f` are dimensionless. `A` is acoustic power reaching the stage-1 junction, in W. `W1` enters pulse tube 1; `W2h` enters the warm end of regenerator 2. The acoustic power reaching pulse tube 2 requires the still-unfinished regenerator transmission law.

Fitting relative resistance strengths is the agreed level of description. The remaining temperature dependence of the effective branch resistances is open. The algebra above is an effective approximation, not a completed acoustic circuit model.

## Shared compressor source: agreed first trial

Represent each identical compressor by fixed effective mean gas delivery `M`, in kg/s, and a fixed effective pressure-amplitude ceiling `p_cap`, in Pa, during cooldown. For `N` parallel compressors, available delivery is `N*M`; the pressure ceiling stays unchanged.

Let `p_a` be the effective peak oscillating pressure amplitude at the common head inlet, in Pa. Let `D(p_a; T, E)` be the total cycle-averaged supply demand of all connected heads, in kg/s, at their current boundary temperatures and stored energies.

```text
If D(p_cap; T, E) <= N*M:
    p_a = p_cap
Otherwise:
    solve D(p_a; T, E) = N*M for p_a below p_cap
```

In the first regime surplus compressor delivery bypasses the heads; in the second, delivery limits the attainable amplitude. This is the agreed approximation to try, not a claim that compressor delivery is physically constant under all conditions. `p_cap` is an effective head-inlet amplitude, not directly the measured compressor high-to-low pressure difference. Mean pressure and oscillating pressure amplitude are different quantities.

`D` describes gas transferred from the high-pressure supply to the return over complete cycles. It is different from oscillating regenerator mass flow, whose cycle average can be zero. Its closure remains open. The last discussion proposed including gas used to compress the internal gas inventory each cycle; that addition was not agreed before the session stopped. Do not silently adopt it or a linear demand law.

The joint source/head calculation must supply local acoustic powers at the regenerator and pulse-tube interfaces, mean-square oscillating mass flows in the regenerators, and mean pressure for helium properties. Acoustic power alone does not determine mean-square mass flow. One versus two compressors must be handled through this same source rule, rather than a universal multiplier on cooling power.

## Local transport derivation retained for completion

The following local reduction was worked through and understood in discussion. It is material for the consolidated transport closure, not a complete global model.

For a local gas/matrix section:

```text
Cg * dTg/dt = S - H*(Tg - Tm)
Cm * dTm/dt = H*(Tg - Tm)
alpha = Cm / (Cg + Cm)
Tg - Tm approximately equals alpha*S/H
```

`Cg` and `Cm` are gas and matrix heat capacities, in J/K. `Tg` and `Tm` are their cycle-dependent temperatures, in K. `H` here is the local gas-matrix exchange conductance, in W/K, distinct from the boundary energy flows above. `S` is local forcing, in W. `alpha` is dimensionless.

For the local gradient-only forcing, the reduction gives:

```text
S = -m_dot * cp * L * dTmean/dx
Jgrad = -Dgrad * dTmean/dx
Dgrad = (cp^2 * L / H) * average(m_dot^2) * [Cm/(Cg+Cm)]^2
```

`m_dot` is oscillating mass flow, in kg/s; `cp` is gas specific heat at constant pressure, in J/(kg K); `L` is section length, in m. `Tmean` is local mean temperature, in K. `Jgrad` is gradient-driven energy transport, in W, and `Dgrad` has units W m/K. The average is over one cycle. This derivation uses the local periodic approximation; completing the global closure requires consistent pressure-driven transport, acoustic transmission, and local flow amplitudes.

## Intended use

1. The user has a fridge with cooldown data and a pulse-tube serial number.
2. They obtain the test-stand cooldown and load map for that head.
3. They specify the compressor count and head configuration, or use a base configuration when these are unknown.
4. They fit parameters to that head's test-stand data.
5. They put the fitted parameters into the fridge configuration and run the simulation.

## Calibration targets

- The transient stage-temperature cooldown curves from room temperature (approximately 300 K) to base temperature (approximately 3 K) or the base-load operating condition.
- The measured load-map fixed points `(Q1, Q2, T1, T2)`, with head heat flows expressed using the positive-into-the-node convention.
- The four base-temperature Jacobian entries `dQ1/dT1`, `dQ1/dT2`, `dQ2/dT1`, and `dQ2/dT2`, expressed as slopes inferred from the head's cold load map.

The human will write acceptance tests for these targets before implementation. Equation and closure development should address all three targets.

The concrete comparison target is [head_comparison/comparison.png](../../simulations/validation/head_comparison/comparison.png): both stage cooldown traces, cold load maps, and all four Jacobian panels. The two-compressor record is a different head, not a controlled same-head compressor-doubling measurement.

Retain the observed sign and temperature dependence of the cross-coupling. In the heat-positive convention, `dQ1/dT2` is positive and weakens as `T1` rises; `dQ2/dT1` is small. See [load-map characterization](../../simulations/exploration/loadmaps/characterization.md). The contrary wording in the older `pulse-tube-model` context is not the target. Cooldown also includes the stage-2 slowdown and the moving balance as stage 1 cools.

## Positions

- Pursue a simpler pulse-tube model from a dynamical-systems and differential-equation perspective.
- Begin a new formulation from the component arrangement and energy balances.
- Represent the entire pulse-tube head as one link connecting PRP, RGP, and CFP, with internal thermal state belonging to the head model.
- The model supplies `(Q1, QR, Q2)` given the temperatures `(T1, TR, T2)`, its internal thermal state, and the configured head and drive parameters.
- The cooling powers may depend directly on other stage temperatures as well as internal state and drive; cross-coupling may be represented by effective laws.
- Include internal thermal state from the outset.
- Retain physical structure in the equations wherever possible, with coefficients allowed to represent bulk behavior.
- Retain the physical arrangement of components: the first regenerator connects the warm flange to stage 1, the second connects stage 1 to stage 2 with an intermediate clamp, and the two pulse tubes connect their respective stages to the warm end.
- Coefficients may describe bulk behavior without representing discrete physical things.
- Derive the differential equations and closures from physics as far as possible; introducing fitted coefficients requires a strategy for obtaining their values from available measurements.
- The reduced model describes slow thermal evolution with physically derived cycle-averaged transport closures. Explicit reconstruction of valve dynamics, pressure waveforms, and acoustic phase relationships is outside scope.
- Small algebraic systems within the model are acceptable. Prefer ODE evolution without DAE constraints on the integrated state.
- Low-temperature property behavior is a central numerical concern, including the sharp Gd2O2S heat-capacity peak near 5 K and helium-property behavior near 2 K.
- The working expectation is that the two stage-temperature cooldown traces constrain roughly two values and the load map roughly four. These are approximate estimates of useful fitting freedom, not established identifiability counts.
- Additional fitted values trade off the interpretability and determinacy of the others.
- Parameters may be fitted to base-temperature load maps and test-stand cooldown curves.
- Calibration targets the individual head's test-stand records. Research on the population of heads passing through the test stand is outside scope.
- Parameters must not be fitted to actual fridge cooldown data.
- Further development of high-fidelity pulse-tube process models is outside this project's scope.
- The model must support different thermal loads on the two stages.
- The model must support heat extraction at the regenerator clamp and its coupling to the two stages.
- The model must support changes in the number of compressors and heads on a manifold. This requirement may be revisited if the lack of test data makes it impractical.
- Preliminary inspection of cooldown curves suggests that stage 2 approaches a balance temperature that slowly declines as stage 1 cools. The model should capture this behavior.

## Open formulation

| Missing piece | Required completion |
| --- | --- |
| Regenerator transport | Assemble the local reductions into one complete energy-flow law, including pressure-dependent enthalpy transport, gradient-driven transport, and conduction. Specify consistent acoustic-power propagation. |
| Head flow and power distribution | Connect inlet excitation, the two gas paths, local acoustic powers, and oscillating mass-flow strengths. Finish the temperature dependence of the effective branch resistances. |
| Pulse-tube losses | Specify how acoustic powers become net removals `P1` and `P2`. |
| Shared source | Define head supply demand and mean-pressure treatment. Decide whether to include the proposed gas-compression contribution. |
| Storage and numerical specification | Finish gas-energy bookkeeping and the energy/profile inversion. Address material and helium-property behavior near the low-temperature limit. |
| Calibration | Enumerate fixed inputs and fitted parameters, and state which observations constrain each parameter. |

### Alternatives not adopted

- An additional effective thermal reservoir: the comparison data do not uniquely identify it. Earlier reservoir estimates depended on inferred heater histories and assumptions about power delivery. Begin with the represented material inventory.
- The fitted inverse-power temperature law for acoustic splitting from `pulse-tube-model`, including its fitted exponent, and the subsequently suggested warm resistance floor. That fit was joint with other closures; its internal split was not directly measured. Its stage-1 cross-curvature also disagreed with the target trend.
- A constant shared source resistance as the complete compressor model. The agreed first trial instead distinguishes pressure-limited and delivery-limited regimes.
- An assumed direct conversion from acoustic power to mean-square oscillating mass flow without a flow/pressure relation.
- Explicit gas-parcel tracking, valve-waveform reconstruction, or cycle-by-cycle integration.

## Four bounded steps to finish

1. **Produce one consolidated candidate equation set.** Put the agreed structure and all missing closures together in one document. For each open closure, give one recommended expression, its physical assumptions, its parameters, and the consequence for cooldown or cross-coupling. The deliverable must calculate every derivative and output from a supplied thermal state. Further discussion should resolve these specific gaps. This is the immediate next deliverable, rather than another isolated mechanism discussion.
2. **Review the complete candidate and parameter count together.** Decide whether its complexity is justified and whether the available measurements can constrain it. Resolve the shared-drive calculation as a whole, including whether gas-compression demand belongs in the first candidate. Finish the human's acceptance criteria before implementation.
3. **Implement a minimal numerical prototype and calibrate it jointly.** Use the same head parameters for its cold load map and cooldown. Calculate cold Jacobians after the internal regenerator energies have relaxed. Check conservation and the intended pressure-limited and delivery-limited behavior.
4. **Evaluate against `comparison.png` and make a model-level decision.** Examine both cooldown traces, all four Jacobian entries and their temperature dependence, and the one/two-compressor behavior. Revisions should address an identified discrepancy. Add storage or fitted laws only when the comparison supplies a concrete reason.

### How to resume the discussion

Continue using the hash-it-out format: explain assumptions and derivations in plain language, put one substantive technical decision to the human at a time, and preserve their ownership of the model. Frame questions as the actual physical or complexity tradeoff, with a recommendation and consequences. Keep the complete equation set and the bounded steps visible while resolving those decisions. Use the two repositories for evidence; distinguish their measured findings from assumptions and fitted laws. Define symbols and use terminal-readable plain-text math.

## Reference work

The new formulation starts afresh. Earlier models supply reference physics and failure findings, rather than an inherited model specification.

- [Archived ADR-03](../adrs/archive/03-coupled-capacity-maps.md): retired after failing to reproduce load maps, cross-coupling, and cooldowns with one parameter set.
- [Archived ADR-04](../adrs/archive/04-gas-wall-shuttling-losses.md): rejected gas-wall shuttling scaling with excessive warm losses; applicability of approximations matters across the full cooldown.
- [Exploratory model](../../../pulse-tube-model/MODEL.md) and [its context](../../../pulse-tube-model/CONTEXT.md): equations, assumptions, and subsequent findings. Its earlier restriction against fitting test-stand cooldowns is superseded for this proposal.
- [Warm-acoustic screens](../../../pulse-tube-model/results/phase_acoustic_check/interpretation.txt): unsuccessful warm continuations and stage-power allocation laws.
- [Storage and warm-power investigation](../../../pulse-tube-model/results/cooldown_model_screen/interpretation.txt) and [reservoir calibration findings](../../../pulse-tube-model/results/cooldown_calibration/B_shield_on_joint/interpretation.txt): storage/power tradeoffs and limitations of attributing an effective reservoir to a particular material mass.

- [Detailed regenerator-memory investigation](../../../pulse-tube-model/results/regenerator_memory_check/report.txt): assesses memory corrections under its fixed candidate laws; does not uniquely identify an additional physical reservoir.
- [Acoustic-split curvature investigation](../../../pulse-tube-model/results/curvature_mechanism_check/interpretation.txt): retains the evidence behind declining to transfer the older fitted split law.
- [Snodgrass et al., dynamic acoustic optimization (2024)](https://pmc.ncbi.nlm.nih.gov/articles/PMC11039843/): physical basis for pressure-limited versus delivery-limited operation and the proposed gas-compression contribution. Its complete circuit model is not adopted.
- [NASA/TM-2006-214391](https://ntrs.nasa.gov/api/citations/20060051853/downloads/20060051853.pdf), Eq. 8: the linear viscous pressure-drop term motivates the `mu/rho` factor when expressed in terms of mass flow. The aggregate resistance remains an effective approximation.
- [Compressor records](../../data/compressors.md): measured pressures, delivery estimate provenance, and manifold/reservoir volumes. The displacement-based delivery estimate is not a measured universal compressor capacity.
