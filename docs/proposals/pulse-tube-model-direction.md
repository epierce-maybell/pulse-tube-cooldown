# Pulse-tube model direction

This proposal records the candidate equation set agreed for handoff on 2026-11-10: five fitted coefficients, seven material-piece energy entries, and all three thermal connections. The candidate has since been implemented and tested; revisions remain under discussion. [DESIGN](../DESIGN.md#open-decisions) records the head-model replacement as in progress.

The objective is a compact dynamical model that can jointly fit the cooldown trajectories, load-map points, and cross-couplings. Physical derivations guide the functional forms and preserve useful constraints; reproducing pulse-tube engineering is not an objective. Added structure must address an identified limitation in representing those targets and explain why it improves on the failed formulations. Prioritize revisions and removals; additions need evidence of promise beforehand. This priority governs the closure choices below, including how far to pursue physical derivations.

The user approaches this as a dynamical-systems and differential-equations problem. The model is intended for operational use with limited development time. Explain decisions in those terms, make their complexity cost and target-relevant benefit explicit, and preserve the user's ownership of them. Additional acoustic structure requires a demonstrated need and an explanation of how it improves on the retired approaches. Uncertainty sweeps and further high-fidelity process modeling are outside the agreed work.

Eight fitted coefficients are too many; a replacement count has not been settled. The user has questioned whether unspecified temperature-dependent functions merely conceal complexity in a reduced formulation. That concern remains unresolved.

The user no longer supports the observation that originally motivated internal state and supports exploring a simpler formulation; the treatment of material storage remains unresolved. The user also expressed a preference for removing delivery from the transient and questioned the compressor inputs. The replacement drive description remains unresolved.

## Model interface

```text
dE/dt = F(T, E, d)
Q = G(T, E, d)
```

`T = (T1, TR, T2)` contains the PRP, RGP, and CFP temperatures in K. `E` contains stored energies in the modeled regenerator sections, in J. `t` is time in seconds. `d` describes the configured drive and shared manifold. `F` gives rates of energy accumulation in W. `Q = (Q1, QR, Q2)` contains heat flows into PRP, RGP, and CFP, in W. `G` calculates those heat flows from the complete state. Boundary conditions, including warm-flange temperature, and fitted parameters are implicit in these functions.

`F` and `G` share the heat flows in the component energy balances. Load maps constrain the steady state, where `F = 0`; test-stand cooldown curves constrain the evolution of internal energy and its effect on cooling power.

Stage temperatures `T1` and `T2` are the head's PRP and CFP connection temperatures; `TR` is the external RGP clamp-plate temperature. All are in K. The head heat flows `Q1`, `QR`, and `Q2` are in W, positive when heating the receiving plate. The fridge model evolves the external plate energies.

A **manifold** is the compressors and heads sharing one gas supply. A **bulk coefficient** describes aggregate behavior without necessarily identifying an individual physical component. The **regenerator matrix** stores and exchanges heat with the gas. Use the term **seven material-piece energy entries** for its representation; these entries are distinct from the five fitted coefficients.

Regenerator 1 connects the warm flange to stage 1. Its cold end feeds pulse tube 1 and regenerator 2. Regenerator 2 connects stage 1 to stage 2, with the clamp as an intermediate thermal side contact. Pulse tube 2 connects stage 2 to the warm end.

## Intended use

1. The user has a fridge with cooldown data and a pulse-tube serial number.
2. They obtain the test-stand cooldown and load map for that head.
3. They specify the compressor count and head configuration, or use a base configuration when these are unknown.
4. They fit parameters to that head's test-stand data.
5. They put the fitted parameters into the fridge configuration and run the simulation.

Calibration uses the individual head's test-stand measurements. Parameters are not fitted to fridge cooldown data; population-level head research is outside this task.

## Calibration targets

S0035 is the first case study. All four cold-map derivatives are critical and must not be sacrificed for transient fits. The model must retain both the clamp temperature `TR` and heat flow `QR`, including when the test-stand condition is `QR = 0`.

- The transient stage-temperature cooldown curves from room temperature (approximately 300 K) to base temperature (approximately 3 K) or the base-load operating condition.
- The measured load-map fixed points `(Q1, Q2, T1, T2)`, with head heat flows expressed using the positive-into-the-node convention.
- The four base-temperature Jacobian entries `dQ1/dT1`, `dQ1/dT2`, `dQ2/dT1`, and `dQ2/dT2`, expressed as slopes inferred from the head's cold load map.

The human will write acceptance tests for these targets before implementation.

The concrete comparison target is [head_comparison/comparison.png](../../simulations/validation/head_comparison/comparison.png): both stage cooldown traces, cold load maps, and all four Jacobian panels. The two-compressor record is a different head, not a controlled same-head compressor-doubling measurement.

Retain the observed sign and temperature dependence of the cross-coupling. In the heat-positive convention, `dQ1/dT2` is positive and weakens as `T1` rises; `dQ2/dT1` is small. See [load-map characterization](../../../pulse-tube-research/simulations/exploration/loadmaps/characterization.md). The contrary wording in the older `pulse-tube-model` context is not the target. Cooldown also includes the stage-2 slowdown and the moving balance as stage 1 cools.

## Candidate equations

### State, temperature profile, and storage

The integrated head variables are `E_j`, in J, for pieces `j = 1,...,7` in this warm-to-cold material order:

1. Regenerator 1: stainless-steel screens.
2. Regenerator 1: lead.
3. Regenerator 2: lead.
4. Regenerator 2: Er3Ni.
5. Regenerator 2: HoCu2 above the clamp.
6. Regenerator 2: HoCu2 below the clamp.
7. Regenerator 2: Gd2O2S.

In each regenerator, position `x` increases from warm to cold. The lengths are `L_R1` and `L_R2`, in m. The warm-flange temperature is `T0`, in K.

`theta_c`, in K, is the internal regenerator temperature at the clamp. Reconstruct piecewise-linear temperature profiles `T_R1(x)` and `T_R2(x)` through the material-piece centers, the end temperatures `(T0, T1)` and `(T1, T2)`, and the internal clamp temperature `theta_c`. The seven center temperatures and `theta_c` satisfy the seven storage equations and the clamp balance below:

```text
E_j = integral_over_piece_j e_j(T(x)) dm
de_j/dT = c_j(T)
```

Here `dm` is matrix mass in kg, `c_j` is matrix specific heat in J/(kg K), and `e_j` is matrix specific energy in J/kg with a fixed arbitrary reference zero. Mass and pore volume are distributed uniformly within each material piece, using the existing inventories. Split HoCu2 mass and pore volume according to the lengths above and below the clamp.

Slow stored energy is the represented matrix energy. Gas heat capacity enters the cycle response below; slow gas and wall storage are neglected in this candidate. A uniform initial temperature `T_init` gives `E_j(0) = mass_j * e_j(T_init)`.

### Mean pressure and local properties

```text
p_mean = prescribed constant
```

`p_mean`, in Pa absolute, is the effective mean pressure used for helium properties throughout a run. Constant mean pressure is an accepted first-candidate approximation; its adequacy for the fitting targets is unestablished. The candidate has no total-charge pressure solve.

Evaluate helium density `rho` in kg/m^3, isobaric specific heat `cp` in J/(kg K), expansion coefficient `beta = -(1/rho)*(partial rho/partial T)_p` in 1/K, and dynamic viscosity `mu` in Pa s at the current local `T` and `p_mean`. Use the existing real-helium and material property data within their supported domains.

For each regenerator:

```text
C_m(x) = (dm/dx) * c_j(T(x))
C_g(x) = rho(T(x), p_mean) * cp(T(x), p_mean) * (dV_pore/dx)
alpha(x) = C_m(x) / (C_m(x) + C_g(x))
g_i = G_i / L_Ri                         for i = 1, 2
```

`C_m` and `C_g` are heat capacities per length in J/(m K). `V_pore` is pore volume in m^3. `alpha` is dimensionless. `G_i`, in W/K, is a fitted constant whole-regenerator gas-matrix exchange conductance; `g_i` is its distributed value in W/(m K).

The cycle calculation uses a linear small-oscillation response, fixed favorable phasing, fast gas-matrix exchange, and properties frozen at the current local mean condition. These approximations apply across the first-candidate cooldown.

### Power distribution and transmission

Reference viscosity `mu_ref` and density `rho_ref` are evaluated at `(T0, p_mean)`. Define the dimensionless profile measure:

```text
s2 = (1/L_R2) * integral_0^L_R2 [(mu/rho) / (mu_ref/rho_ref)] dx
```

The suffix `_cap` denotes power evaluated at the prescribed pressure-amplitude ceiling `p_cap`, in Pa. For each head:

```text
W1_cap  = W_star
W2h_cap = W_star / (b + a*s2)
A_cap   = W1_cap + W2h_cap

tau_i(x) = exp(integral_0^x alpha(s)*beta(s)*(dT_Ri/ds) ds)

W_R1_cap(x) = A_cap * tau_1(x) / tau_1(L_R1)
W_R2_cap(x) = W2h_cap * tau_2(x)

W_in_cap = A_cap / tau_1(L_R1)
W2c_cap  = W2h_cap * tau_2(L_R2)
```

`W_star` is the fitted power scale in W; `a` and `b` are dimensionless fitted branch coefficients. `A_cap` is power at the stage-1 junction. `W1_cap` enters pulse tube 1, `W2h_cap` enters regenerator 2, `W2c_cap` reaches pulse tube 2, and `W_in_cap` enters regenerator 1. `W_Ri_cap(x)` is local power within regenerator `i`. All powers are in W; `tau_i` is dimensionless. The integration variable `s` is position in m.

Equivalently, the fraction entering regenerator 2 is `1/(1+b+a*s2)`. The branch coefficients describe effective allocation. The transmission law is `dW/dx = alpha*beta*W*dT/dx`; the candidate neglects explicit viscous and thermal-relaxation dissipation in that law.

### Shared source and oscillating flow

Each identical compressor has prescribed effective mean delivery `M`, in kg/s. There are `N` compressors and a configured set of connected heads. Each head has its specified valve frequency `nu`, in Hz. With `rho0 = rho(T0, p_mean)`, calculate each head's demand at the pressure ceiling:

```text
C_gas = integral_over_active_head_gas (partial rho/partial p)_T dV

D_pass = 2*rho0*W_in_cap / (pi*p_cap)
D_fill = 2*nu*p_cap*C_gas
D_cap  = D_pass + D_fill

z   = min(1, N*M / sum_heads(D_cap))
p_a = z*p_cap
```

`C_gas`, in kg/Pa, is the isothermal mass response to pressure. It is distinct from the heat capacity `C_g`. `D_pass`, `D_fill`, and `D_cap` are average supply demands in kg/s. `z` is the common dimensionless pressure fraction; `p_a` is effective peak oscillating pressure amplitude in Pa. For a disabled drive, `z = 0`.

The active inventory comprises regenerator pores, heat-exchanger voids, pulse tubes, and valve-to-head dead volume. Pore temperatures follow the reconstructed profiles; exchanger temperatures equal their stage temperatures; pulse-tube temperatures vary linearly from their stage to `T0`; inlet dead volume is at `T0`. Terminal reservoir pressure oscillations are neglected. Supply/return lines and terminal reservoirs do not enter `C_gas`.

Through-flow and filling demands add as an effective approximation. Filling is evaluated at fixed local temperature. The through-flow factor `1/pi` converts peak sinusoidal mass flow to its positive-flow average over a complete cycle. At a fixed thermal state, both demands scale linearly with `p_a`.

All actual powers scale with `z^2`:

```text
W1       = z^2 * W1_cap
W2h      = z^2 * W2h_cap
W2c      = z^2 * W2c_cap
W_Ri(x)  = z^2 * W_Ri_cap(x)

mean(m_dot_i(x)^2) = 2*z^2 * (rho(T_Ri(x), p_mean)*W_Ri_cap(x)/p_cap)^2
```

`m_dot_i` is local oscillating mass flow in kg/s; `mean` averages over one cycle. This relation uses the common effective pressure amplitude and in-phase pressure/flow approximation. Compression contributes to supply demand only; its local oscillating-flow contribution is neglected in the transport closure.

### Regenerator energy-flow closure

For either regenerator, define:

```text
Phi(x)  = (1 - alpha(x)*beta(x)*T(x)) * W(x)
Dgrad(x) = alpha(x)^2 * cp(x)^2 * mean(m_dot(x)^2) / g_i
K(x)    = Dgrad(x) + A_Rwall_i*k_Rwall_i(T(x))

H = Phi - K*dT/dx
```

`Phi` is pressure-dependent energy transport in W. `Dgrad` and `K` have units W m/K. `A_Rwall_i` is regenerator wall area in m^2 and `k_Rwall_i` is its thermal conductivity in W/(m K). Axial conduction uses the regenerator wall alone in this candidate.

Use one shared discrete energy flow between neighboring profile knots `u` and `v`, ordered warm to cold:

```text
H_uv = [T_u - T_v + integral_u^v (Phi/K) dx] / integral_u^v (1/K) dx
```

Here `u` and `v` are positions in m and `T_u`, `T_v` are their temperatures in K. Evaluate the integrals on the reconstructed profile, split at material boundaries. This expression specifies the discrete flux approximation. Between material centers it supplies the shared material-interface flow; at each regenerator end it connects the boundary to the nearest center. The clamp is a separate knot with a flow on each side. `H_uv` is in W, positive coldward.

### Seven energy derivatives and three outputs

Let `H_12`, `H_34`, `H_45`, and `H_67` denote the shared flows between the corresponding material-piece centers. `H_5c` runs from piece 5 to the clamp and `H_c6` from the clamp to piece 6. `H1h` and `H1c` are the warm- and cold-end flows of regenerator 1; `H2h` and `H2c` are those of regenerator 2. Every flow is calculated by the same discrete closure above.

```text
dE_1/dt = H1h  - H_12
dE_2/dt = H_12 - H1c
dE_3/dt = H2h  - H_34
dE_4/dt = H_34 - H_45
dE_5/dt = H_45 - H_5c
dE_6/dt = H_c6 - H_67
dE_7/dt = H_67 - H2c

QR = K_R*(theta_c - TR)
H_5c - H_c6 = QR
```

`K_R`, in W/K, is the prescribed finite contact conductance. The clamp has no separate stored-energy entry. Solve its balance together with the seven storage/profile relations, using the same drive calculation as the rest of the head.

Pulse-tube wall heat leaks and net removals are:

```text
L_i = (A_PTwall_i / length_PTi) * integral_Ti^T0 k_PTwall_i(T) dT

P1 = W1  - L_1
P2 = W2c - L_2

Q1 = H1c - H2h - W1  + L_1
QR = K_R*(theta_c - TR)
Q2 = H2c       - W2c + L_2
```

For stage `i = 1,2`, `A_PTwall_i` is pulse-tube wall area in m^2, `length_PTi` is tube length in m, and `k_PTwall_i` is wall conductivity in W/(m K). `L_i`, `P1`, and `P2` are in W. The pulse-tube loss closure consists of these calculated wall heat leaks, with no additional fitted loss coefficient.

All three `Q` outputs are positive into the receiving plate. The supplied `TR` changes the internal profile and therefore transport, power delivery, and source demand. On the unloaded test stand, impose `QR = 0`, so `H_5c = H_c6` and `TR = theta_c` for positive finite `K_R`; the clamp temperature floats.

With `E_R1 = E_1 + E_2` and `E_R2 = E_3 + ... + E_7`, the equations give:

```text
dE_R1/dt = H1h - H1c
dE_R2/dt = H2h - H2c - QR

Q1 + QR + Q2 + sum_j(dE_j/dt) = H1h - W1 - W2c + L_1 + L_2
```

## Parameters and status

The five fitted coefficients per head are:

| Coefficient | Units | Constraint | Role |
| --- | --- | --- | --- |
| `W_star` | W | positive | Power scale |
| `a` | dimensionless | nonnegative | Profile-dependent branch allocation |
| `b` | dimensionless | positive | Remaining branch allocation |
| `G1` | W/K | positive | Regenerator-1 transport scale |
| `G2` | W/K | positive | Regenerator-2 transport scale |

Fixed inputs are the geometry, material inventories and properties, helium properties, warm temperature `T0`, valve frequency `nu`, effective mean pressure `p_mean`, pressure ceiling `p_cap`, delivery `M`, compressor count `N`, connected heads, and clamp conductance `K_R`. Numerical values for the effective pressure inputs and `K_R` were not finalized in this discussion. Unloaded stand data cannot identify `K_R`.

The same five fitted values apply jointly to a head's stand cooldowns and cold load map. Cold-map derivatives are taken after all internal energies and the unloaded clamp have relaxed. The acceptance targets remain both cooldown traces, load-map points, and all four cross-couplings in `comparison.png`. A satisfactory joint fit to those targets has not been established.

The fitting budget reflects the limited information in the two transient curves and the load map. No identifiability claim has been established. Known material inventories are included from the outset; extra effective storage remains an optional revision tied to a demonstrated discrepancy.

## Alternatives not adopted

- An additional effective thermal reservoir: the comparison data do not uniquely identify it. Earlier reservoir estimates depended on inferred heater histories and assumptions about power delivery. Begin with the represented material inventory.
- The fitted inverse-power temperature law for acoustic splitting from `pulse-tube-model`, including its fitted exponent, and the subsequently suggested warm resistance floor. That fit was joint with other closures; its internal split was not directly measured. Its stage-1 cross-curvature also disagreed with the target trend.
- A constant shared source resistance as the complete compressor model. The agreed first trial instead distinguishes pressure-limited and delivery-limited regimes.
- An assumed direct conversion from acoustic power to mean-square oscillating mass flow without a flow/pressure relation.
- Explicit gas-parcel tracking, valve-waveform reconstruction, or cycle-by-cycle integration.
- A separate fitted stage-2 pulse-tube loss: the first candidate uses calculated pulse-tube wall heat leaks.
- A total-charge solve for changing mean pressure: the first candidate uses constant effective mean pressure. A rough pressure-change estimate did not establish a corresponding error in the fitting targets.

## Reference work

The new formulation starts afresh. Earlier models supply reference physics and failure findings, rather than an inherited model specification.

- [Archived ADR-03](../adrs/archive/03-coupled-capacity-maps.md): retired after failing to reproduce load maps, cross-coupling, and cooldowns with one parameter set.
- [Archived ADR-04](../adrs/archive/04-gas-wall-shuttling-losses.md): rejected gas-wall shuttling scaling with excessive warm losses; applicability of approximations matters across the full cooldown.
- [Exploratory model](../../../pulse-tube-model/MODEL.md) and [its context](../../../pulse-tube-model/CONTEXT.md): equations, assumptions, and subsequent findings. Its earlier restriction against fitting test-stand cooldowns is superseded for this proposal.
- [Warm-acoustic screens](../../../pulse-tube-model/results/phase_acoustic_check/interpretation.txt): unsuccessful warm continuations and stage-power allocation laws.
- [Storage and warm-power investigation](../../../pulse-tube-model/results/cooldown_model_screen/interpretation.txt) and [reservoir calibration findings](../../../pulse-tube-model/results/cooldown_calibration/B_shield_on_joint/interpretation.txt): storage/power tradeoffs and limitations of attributing an effective reservoir to a particular material mass.

- [Detailed regenerator-memory investigation](../../../pulse-tube-model/results/regenerator_memory_check/report.txt): assesses memory corrections under its fixed candidate laws; does not uniquely identify an additional physical reservoir.
- [Acoustic-split curvature investigation](../../../pulse-tube-model/results/curvature_mechanism_check/interpretation.txt): retains the evidence behind declining to transfer the older fitted split law.
- [Snodgrass et al., dynamic acoustic optimization (2024)](https://pmc.ncbi.nlm.nih.gov/articles/PMC11039843/): physical basis for pressure-limited versus delivery-limited operation and gas-compression demand. Its complete circuit model is not adopted.
- [NASA/TM-2006-214391](https://ntrs.nasa.gov/api/citations/20060051853/downloads/20060051853.pdf), Eq. 8: the linear viscous pressure-drop term motivates the `mu/rho` factor when expressed in terms of mass flow. The aggregate resistance remains an effective approximation.
- [Compressor records](../../data/compressors.md): measured pressures, delivery estimate provenance, and manifold/reservoir volumes. The displacement-based delivery estimate is not a measured universal compressor capacity.
