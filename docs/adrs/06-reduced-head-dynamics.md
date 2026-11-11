# ADR-06: Reduced nonlinear dynamical system model of the pulse-tube head

Status: accepted, 2026-11-10.

## Decision

Represent each head as a reduced nonlinear dynamical system. Its internal state
consists of seven material-piece energies. Stage and clamp temperatures are
inputs; the corresponding cooling powers are outputs. A vector field governs
the internal state, and an algebraic output map gives cooling power from that
state and the temperature inputs. Coupling the head to the external node energy
balances gives the assembled system of differential equations.

The reduction retains the slow evolution of stored thermal energy. Material
inventories and specific heats determine the energy-temperature relations.
Measured performance data determine functions and parameters in the vector field
and output map, including the cycle-averaged effects of acoustic and gas dynamics.

The head's measured cold load map anchors equilibrium performance and cold
cross-coupling. A prescribed power law supplies each stage's warm correction.
Stand cooldowns calibrate seven scalar coefficients: a warm amplitude and
exponent for each stage, an exchange coefficient for each regenerator, and one
coefficient coupling stage-1 performance to departure from internal equilibrium.
The functional form is shared across heads and drive conditions; the seven
coefficients are fitted separately for each head and drive condition.

Represent clamp loading by a load coordinate, published temperature-response
tables, loaded regenerator profiles, and a configured external contact
resistance. Zero clamp load recovers the same unloaded stage surface and storage
equations. The stand's unmeasured clamp temperature follows from that condition.

Adopt this formulation for production implementation. Acceptance uses the
available load-map points, all four cell-derived stage responses, cooldown
traces, and physical energy balances at the fidelity of the other simulator
components. Calibration applies to the head and drive conditions represented by
its records. The human authors the acceptance tests.

## Dynamical system

### State, inputs, and outputs

```text
E = (E_1, ..., E_7)
T = (T1, TR, T2)
Q = (Q1, QR, Q2)

dE/dt = f(E,T)
Q = g(E,T)
```

`E` is the internal state, with each component in J. `T` is the temperature
input vector, in K, and `Q` is the cooling-power output vector, in W. The vector
field `f` returns seven energy rates in W; the output map `g` returns three
cooling powers in W. The head data, calibration parameters, and prescribed
warm-flange temperature `T0` are fixed arguments of these functions. Internal
equilibrium satisfies `f(E,T) = 0`.

`T1` and `T2` are the head's stage temperatures; `TR` is the external clamp-plate
temperature. Time `t` is in s. `Q1`, `QR`, and `Q2` are delivered cooling powers
in W, positive for heat removal. The link returns heat flows `-Q1`, `-QR`, and
`-Q2` under the engine's positive-into-the-node convention.

The refrigerator connections are stage 1 to PRP, stage 2 to CFP, and the clamp
to RGP. Clamp load `ell`, in W, is an algebraic variable with `QR = ell`.
It is determined by the contact equation below and eliminated when evaluating
`f` and `g`. The seven energies carry the head's thermal history.

## Vector field and output map

### Internally equilibrated cooling response

```text
q0_i(T1,T2) = q_base,i(T1,T2) + w(T1,T2) * s_i(T_i)
q_base,i(T1,T2) = H_i(B1(T1),B2(T2)) + L_stand,i(T1,T2)
s_i(T_i) = A_i * (T_i / T0)^b_i
```

For stage `i = 1,2`, `q0_i` is internally equilibrated cooling power at zero clamp load,
`H_i` interpolates measured heater power, and `L_stand,i` is the modeled stand
parasitic heating, all in W. Adding the stand parasitics converts heater power
into inferred head cooling power. The resulting head performance is used with
the target refrigerator's own external loads.

Interpolate the cold heater-power map exactly with a thin-plate radial-basis
function and a linear polynomial, normalizing each temperature coordinate by
its measured span. The continuation `B_i` equals `T_i` inside that stage's
measured range and approaches a finite limit outside it through exponential
tails tangent to the identity at each boundary. Tail lengths are configured
relative to the measured temperature spans.

`s_i` is the warm correction in W. Its positive amplitude `A_i`, in W, is the
correction at `T_i = T0`; its positive dimensionless exponent `b_i` sets the
temperature dependence. `T0` is the prescribed warm-flange temperature in K.
At fixed `T0`, a small fractional change in `T_i` produces `b_i` times that
fractional change in `s_i`. The power law is monotone for positive temperatures.

The dimensionless blending function `w` is:

```text
D_i = T_max,i - T_min,i
outside_i = max(T_min,i - T_i, 0, T_i - T_max,i) / D_i
w(T1,T2) = 1 - exp(-(outside_1^2 + outside_2^2) / beta^2)
```

`T_min,i` and `T_max,i` are the measured cold-map temperature bounds in K;
`D_i` is their span in K; `outside_i` is the dimensionless distance outside that
stage's interval; and `beta` is a positive dimensionless blend width.
Thus `w` is zero throughout the cold-map temperature bounding box and approaches
one outside it. This construction preserves the cold map and its first
derivatives while supplying a prescribed continuation of coupling into the
warm region. The demonstrated fits use `beta = 1` and continuation tail lengths
equal to the respective measured spans, held fixed during calibration.

Throughout that bounding box, including its boundary, `w = 0` and
`partial(w)/partial(T_j) = 0` for either stage temperature `T_j`, `j = 1,2`.
Consequently:

```text
q0_i = q_base,i
partial(q0_i)/partial(T_j) = partial(q_base,i)/partial(T_j)
```

The power-law calibration preserves the complete cold surface, all four entries
of its temperature Jacobian, and the cell-derived responses of the earlier
cold-map interpolant. This is the relaxed, zero-clamp-load response: internal
energies follow their equilibrium values and the free clamp temperature follows
the unloaded condition. At equilibrium the storage rates vanish and `p = 1`,
so the fitted exchange and internal-energy coupling coefficients also leave
this response unchanged.

### Clamp-load response

The input tables contain the changes in stage-1, stage-2, and clamp temperatures
in a published fixed-stage-load record, relative to its zero-clamp-load point.
`d_1(ell)` and `d_2(ell)`, in K, interpolate the two stage changes. The clamp
temperature increment `d_c(ell)`, also in K, accounts for the source stages'
temperature changes through the unloaded profile. All three functions equal
zero at `ell = 0`. The default input is Falter et al.'s Fig. 4 curve with stage
heater powers of 20 W and 0.5 W. This transfers the measured response shape to the
particular head's unloaded surface.

```text
u1 = T1 - d_1(ell)
u2 = T2 - d_2(ell)

F_c = F(x_R)
d_c_at_source_knot = Delta_TR_source
                     - F_c * Delta_T1_source - (1 - F_c) * Delta_T2_source
r_eq(T1,T2,ell) = F_c * T1 + (1 - F_c) * T2 + d_c(ell)
```

`u1` and `u2` are the unloaded-equivalent temperature coordinates in K.
`x_R` is the configured normalized clamp position, `F` is the unloaded profile
fraction defined below, and `r_eq` is the relaxed head-side clamp temperature
in K. `Delta_T1_source`, `Delta_T2_source`, and `Delta_TR_source` are the
published temperature changes at a source load knot, in K. The published
differences preserve the load dependence along the chosen source curve. Their
transfer to other temperatures and heads is a prescribed approximation.
Interpolate the stage differences with shape-preserving cubic Hermite curves.
For `d_c`, use shape-preserving interior slopes and the first and last interval
secants as endpoint slopes. Continue all three curves linearly along their
endpoint tangents outside the tabulated load range.

### State evolution and target profiles

The seven material-piece energy entries, ordered from warm to cold, are:

1. Regenerator 1: stainless steel.
2. Regenerator 1: lead.
3. Regenerator 2: lead.
4. Regenerator 2: Er3Ni.
5. Regenerator 2: HoCu2 above the clamp.
6. Regenerator 2: HoCu2 below the clamp.
7. Regenerator 2: Gd2O2S.

Masses, material properties, layer lengths, and clamp position come from the
configured head design. Mass is distributed uniformly within each layer; the
HoCu2 inventory is split by its lengths above and below the clamp. External
stage heat capacities and matrix storage use separate inventories.

For piece `j`, define:

```text
E_j = m_j * e_j(theta_j)
de_j/dtheta = c_j(theta)

E_eq,j(T1,T2,ell) = integral_over_piece_j e_j(T_profile(x)) dm
theta_eq,j = inverse_e_j(E_eq,j / m_j)

dE_j/dt = m_j * g_r * (theta_eq,j - theta_j)
```

`E_j` and `E_eq,j` are actual and target energies in J; `m_j` is mass in kg;
`e_j` is material specific energy in J/kg; and `c_j` is material specific heat
in J/(kg K). `theta_j` is the temperature with the piece's actual mean specific
energy, and `theta_eq,j` is its target counterpart. `dm` is an element of matrix
mass in kg. `g_r` is the positive fitted exchange coefficient for the piece's
regenerator `r = 1,2`, in W/(kg K).

Each `g_r` is effective thermal exchange per unit matrix mass. Near a fixed
target temperature, the local relaxation time is approximately
`tau_j = c_j(theta_eq,j) / g_r`, in s. The two coefficients allow the two
regenerators to exchange stored energy at different rates.

Positive `c_j` makes the relation between `E_j` and `theta_j` invertible.
The caloric temperatures `theta_j` are alternative coordinates for the same
seven-dimensional internal state.

Within each regenerator, dimensionless position `x` runs from zero at the warm
end to one at the cold end. The target temperature profile of regenerator 1 is:

```text
T_profile,R1(x) = (1 - x) * T0 + x * T1
```

`F(x)` is the dimensionless warm-temperature fraction obtained by normalizing
Falter et al.'s simulated unloaded regenerator profile in Fig. 7 to its endpoint
temperatures. The same normalized unloaded shape applies across the cooldown.

For regenerator 2, use two dimensionless profile functions `U(x,ell)` and
`V(x,ell)`. `U` goes from one at stage 1 to zero at the clamp; `V` goes from one
at the clamp to zero at stage 2:

```text
T_profile,R2(x) = r_eq + (T1 - r_eq) * U(x,ell)  for x <= x_R
T_profile,R2(x) = T2 + (r_eq - T2) * V(x,ell)    for x >= x_R

U(x,0) = (F(x) - F_c) / (1 - F_c)
V(x,0) = F(x) / F_c
```

The zero-load definitions recover `F(x)*T1 + (1-F(x))*T2` exactly. The loaded
shape inputs are the normalized 1, 2, and 3 W profiles in Falter's Fig. 7, with
the source position mapped separately on each side of the clamp to the
configured clamp position. Normalize each source profile between its own
stage and clamp temperatures. Interpolate the shapes in load with
shape-preserving cubic Hermite curves, setting endpoint slopes to zero and
holding endpoint shapes beyond the tabulated range. The resulting target
profile passes through both stage temperatures and `r_eq`.

A uniform initial temperature `T_init` sets
`E_j(0) = m_j * e_j(T_init)` for every piece.

### Algebraic clamp relation

Use the HoCu2 piece immediately below the clamp, piece `h = 6`, to carry local
thermal memory in the clamp temperature:

```text
theta_c(T1,T2,E,ell) = r_eq(T1,T2,ell) * theta_h / theta_eq,h(T1,T2,ell)
TR = theta_c(T1,T2,E,ell) + R_contact * ell
QR = ell
```

`E` denotes the seven material-piece energies. `theta_c` is the head-side clamp
temperature in K. Its prescribed caloric-temperature ratio represents local
matrix disequilibrium and becomes one at internal equilibrium. `R_contact` is the
positive configured thermal resistance between the external clamp plate and
the modeled head-side clamp, in K/W. This contact is included in the head's
clamp port. Its value comes from the mounted system; unloaded stand records
leave it unconstrained.

Given `E` and `T`, solve this scalar equation for `ell`:

```text
0 = theta_c(T1,T2,E,ell) + R_contact * ell - TR
```

On an operating branch where
`partial(theta_c)/partial(ell) + R_contact != 0`, the relation defines `ell`
locally as a function of `E` and `T`. Substitution gives the vector field `f`
and output map `g`. The load changes both stage performance and the
regenerator-2 target energies, while the internal state changes the clamp
temperature through the caloric-temperature ratio.

### Dependence of the output on internal state

```text
Delta_E_R1 = sum_over_regenerator_1(E_j - E_eq,j)
p = exp(-k * Delta_E_R1 / E_scale)

Q1 = p * q0_1(u1,u2) + sum_j A_1j * dE_j/dt
QR = ell
Q2 =     q0_2(u1,u2) + sum_j A_2j * dE_j/dt
```

`Delta_E_R1` is excess regenerator-1 energy in J. The nonnegative fitted
coefficient `k` and performance multiplier `p` are dimensionless. `E_scale` is
a fixed positive energy normalization in J, stored with the calibration. The
demonstrated fits use `E_scale = 70289.80061191239 J` for all six heads. This
reference comes from the S0035 regenerator-1 energy difference between its
uniform initial state and the target profile at its measured final temperatures.
It remains fixed when fitting or evaluating another record.

For `k > 0`, the equivalent energy scale `E_penalty = E_scale / k`, in J,
gives `p = exp(-Delta_E_R1 / E_penalty)`. An excess energy of `E_penalty`
reduces this stage-1 cooling contribution to `exp(-1)` of its internally
equilibrated value. At `k = 0`, `p = 1`.

`A_ij` is the dimensionless fraction of piece `j`'s energy rate assigned to
stage `i`; the single-index `A_i` above is a warm-correction amplitude in W.
For a piece centered at normalized position `x_bar,j`, use:

```text
Regenerator 1: A_1j = x_bar,j;     A_2j = 0
Regenerator 2: A_1j = 1 - x_bar,j; A_2j = x_bar,j
```

The remaining fraction `1 - A_1j - A_2j` is assigned to the warm flange.
The geometric fractions are a prescribed approximation for allocating transient
matrix-storage power, chosen independently of the target temperature profiles.
Falter's steady measurements and REGEN results constrain the equilibrated
clamp-load response; they leave these transient routing coefficients
undetermined. At internal equilibrium, `dE_j/dt = 0` for every piece, so the
routing terms vanish and steady clamp-load coupling follows from the shifted
cooling surface.

A cooling matrix has negative `dE_j/dt`, reducing the cooling delivered to the
stages. Positive excess regenerator-1 energy also reduces stage-1 performance
through `p`.

The energy-rate terms in the output equations are the components of `f(E,T)`.
Thus `g` is evaluated from the current internal state and temperature inputs.

## Coupling to the external system

The external stage and clamp-plate temperatures evolve according to:

```text
C_i(T_i) * dT_i/dt = P_external,i - Q_i
C_R(TR) * dTR/dt = P_external,R - QR
```

`C_i` and `C_R` are the external stage and clamp-plate heat capacities in J/K,
calculated from their material inventories. `P_external,i` and `P_external,R`
are their external heating powers in W, including applied loads. For one head
and its three attached nodes, the combined state is:

```text
z = (T1, TR, T2, E_1, ..., E_7)
dz/dt = F_system(z,P_external)
```

`P_external` collects the three external heating powers, and `F_system` is the
assembled vector field defined by the node balances and `f`. With fixed
boundary conditions and state-dependent loads, this is an autonomous system.
Explicitly time-dependent loads act as forcing. Additional refrigerator nodes
extend the assembled state through their own energy balances.

At a fixed point, internal energies and external temperatures are stationary:
`f(E,T) = 0` and each external heating power equals its corresponding cooling
power. Internal equilibrium gives `p = 1`.

### Unloaded stand condition

The stand imposes the algebraic constraint `QR = ell = 0`. Its dynamical state
is `(T1,T2,E_1,...,E_7)`, and the free clamp temperature is:

```text
TR_free(T1,T2,E) = theta_c(T1,T2,E,0)
```

At zero load, `u1=T1`, `u2=T2`, and both target profiles recover the unloaded
formulation. All stage cooling powers and storage rates therefore recover the
same functions for every internal state throughout the cooldown. At internal
equilibrium, the measured cold-map balances and their cell-derived slopes are
retained. The stand's pulse-tube thermometer is a separate measurement.

Cold-map derivatives restrict the input-output relation to internal equilibrium
`f=0` and zero clamp load `QR=0`. Eliminating the free clamp temperature gives
the Schur complement of the clamp entry in the relaxed input-output Jacobian:

```text
dQ_i_unloaded/dT_j = partialQ_i/partialT_j
                     - (partialQ_i/partialTR)
                       * (partialQR/partialT_j) / (partialQR/partialTR)
```

Here `i,j = 1,2`, and all partial derivatives are of the internally equilibrated
three-output map. This is the derivative along the unloaded condition; the
external clamp temperature follows the stage temperatures through that condition.

## Calibration and use

### Seven fitted coefficients

The following table gives the physical interpretation and the final S0035
calibration to temperature-crossing times. Values are rounded here; the linked
calibration artifact retains full precision.

| Coefficient | Units | Meaning | S0035 |
|---|---|---|---:|
| `A_1` | W | Stage-1 warm-correction scale at `T1 = T0`. | 144.789 |
| `b_1` | dimensionless | Stage-1 warm-correction temperature exponent. | 1.42393 |
| `A_2` | W | Stage-2 warm-correction scale at `T2 = T0`. | 26.5995 |
| `b_2` | dimensionless | Stage-2 warm-correction temperature exponent. | 0.667621 |
| `g_1` | W/(kg K) | Effective exchange per unit mass in regenerator 1. | 0.0376894 |
| `g_2` | W/(kg K) | Effective exchange per unit mass in regenerator 2. | 0.120623 |
| `k` | dimensionless | Sensitivity of stage-1 performance to excess regenerator-1 energy. | 1.97665 |

The four warm coefficients set the scale and shape of each stage's dependence
on `T_i / T0`. The amplitudes describe the correction before blending, internal
state effects, and storage loads are applied. The exchange and coupling
coefficients describe effective behavior of the reduced system. Their physical
interpretations do not establish that a cooldown uniquely identifies each one.

### Forward calibration

Fit the seven coefficients by integrating the unloaded stand equations with
`ell = 0`, starting from a uniform initial temperature. The measured cold map,
material inventories, target profiles, routing fractions, stand loads, blend
and continuation scales, and energy normalization remain fixed. The simulated
temperatures and internal energies determine cooling power throughout the
integration. Measured temperatures supply comparison targets.

For the demonstrated fits, set `T_init = T0` using the last available prestart
pulse-tube thermometer reading as the initial room-temperature proxy. Exclude
nonfinite readings, stage-1 readings at or above 219 K, and stage-2 readings
above `T0` from the comparison. The source records and reader supply timing and
channel selection. The pulse-tube thermometer and the modeled clamp remain
distinct quantities.

Calibrate to the first downward crossings of these temperatures:

- Stage 1: 200, 150, 100, 80, 65, 55, and 50 K.
- Stage 2: 200, 100, 50, 30, 10, 5, and 3 K.

Crossing times use linear interpolation between the retained samples. For each
crossing `n`, minimize the sum of squared relative timing residuals:

```text
r_n = (t_model,n - t_measured,n) / t_measured,n
objective = sum_over_14_crossings(r_n^2)
```

`t_model,n` and `t_measured,n` are predicted and measured elapsed times in s;
`r_n` and the objective are dimensionless. Each crossing receives equal relative
weight. Earlier fits to the temperature traces supply initial coefficient
guesses. Full temperature traces are also compared after calibration.

The cold map constrains the equilibrium input-output relation and its measured
slopes. The cooldown constrains state evolution and the dependence of cooling
power on departure from equilibrium. The fitted coefficients, prescribed
profiles, and energy normalization are stored head data. All heads use the
same equations with separately calibrated coefficients.

The clamp temperature-response tables and loaded profile shapes are prescribed
design data with their published source and digitization recorded. The external
contact resistance is a system input. The unloaded stand calibration determines
the seven scalar coefficients independently of the chosen contact resistance.

### Demonstrated calibration results

The six selected head records give the following forward results with their
own fitted coefficients. Temperature errors are root mean square differences
over the retained observations; the timing column is the largest absolute
relative error among that head's fourteen crossings.

| Head | Stage-1 temperature error (K) | Stage-2 temperature error (K) | Worst crossing-time error |
|---|---:|---:|---:|
| S0007 | 4.05 | 4.62 | 8.4% |
| S0009 | 1.65 | 5.62 | 12.5% |
| S0019 | 1.64 | 3.88 | 3.3% |
| S0028 | 2.86 | 6.73 | 7.8% |
| S0035 | 1.82 | 4.66 | 3.4% |
| S0038 | 0.98 | 6.11 | 2.7% |

Of the 84 crossings, 83 are within 10% of the measured time. The remaining
crossing is S0009 stage 2 at 10 K: 6197.8 s predicted against 7083.2 s measured,
12.5% early. Reconstructed cold-map fixed-point temperatures differ from the
measured knots by at most `7.2e-14 K`. These are calibration checks on the fitted
records. Generalization to other records, loaded-clamp operation, and assembled
refrigerators requires separate validation against the human-authored criteria.

The S0028 fit places `A_1` at about 1999.5 W, close to the optimizer's 2000 W
upper search bound. Its successful timing fit leaves that coefficient poorly
established. Search bounds and optimizer settings are recorded with the research
scripts and are numerical calibration choices.

### Transient shape and coefficient interpretation

The small stage-2 shoulders examined in the earlier S0035 temperature fit are
changes in cooling rate while temperature continues to fall. The power-law
corrections themselves are monotone. The shoulders arise in the coupled model
as matrix-storage loads change and the warm correction blends into the cold
response. An independent, tighter DOP853 integration agreed with LSODA to within
`5.5e-5 K` across both stage traces.

Changing the prescribed blend width moved the fine transient features. Refitting
the same seven coefficients while omitting 5100-6600 s from the temperature
objective changed the stage-2 prediction in that interval by less than 0.9 K,
while some individual coefficients changed substantially. These diagnostics
support the usefulness of the fixed functional form for cooldown prediction;
the detailed shoulders and individual coefficient values have more limited
physical identification. The diagnostic artifacts apply to the earlier
temperature fit; the table above reports the final crossing-time fits.

## Scope and implementation

The vector field and output map include material-energy relations, calibrated
coefficients, and interpolated response functions. The prescribed continuation
of those functions defines the model outside the measured cold region. Available
records constrain their combined behavior at the simulator's component fidelity.

The model's three-input, three-output formulation includes the clamp-load
response, local thermal memory, loaded target profiles, and external contact
closure. The current prescribed clamp data come from Falter's integrated heat
exchanger; an external clamp is represented through its configured contact
resistance. The transfer of those published response shapes and their
continuation into the warm region remain model assumptions. The algebraic
zero-load limit recovers the unloaded stage equations for every internal state.
The six-head results above evaluate that unloaded model with the power-law warm
response. The earlier three-port consistency report used the preceding
tabulated warm response and documents the clamp construction at that stage of
the investigation.

Task 3.3 is complete with this accepted formulation. Task 3.4 implements it
against the human-authored acceptance tests; task 3.5 supplies the archetypal
head. Head-dependent values and calibration choices belong in configuration or
library data with their sources.

## References

- [Head comparison records](../../simulations/validation/head_comparison/records.toml)
  and [stand records](../../data/test_stand.md).
- [Power-law model and forward integration](../../../pulse-tube-research/simulations/exploration/parametric_warm_response/run.py)
  and [crossing-time calibration](../../../pulse-tube-research/simulations/exploration/parametric_warm_response/fit_crossings.py).
- [Investigation guide](../../../pulse-tube-research/simulations/exploration/parametric_warm_response/README.md),
  including artifact locations, calibration settings, and reproduction commands.
- [Six-head results](../../../pulse-tube-research/simulations/exploration/parametric_warm_response/crossing_fits/summary.json),
  [comparison figure](../../../pulse-tube-research/simulations/exploration/parametric_warm_response/crossing_fits/comparison.png),
  and [S0035 fitted coefficients](../../../pulse-tube-research/simulations/exploration/parametric_warm_response/crossing_fits/S0035/fit.json).
- [S0035 transient diagnostics](../../../pulse-tube-research/simulations/exploration/parametric_warm_response/diagnostics/),
  including solver, structural-sensitivity, and omitted-interval checks.
- [Supporting inventory, profile, and cold-map implementation](../../../pulse-tube-research/simulations/exploration/empirical_head/run.py).
- [Earlier three-port realization](../../../pulse-tube-research/simulations/exploration/empirical_head/clamp_extension/verify_clamp.py)
  and [consistency report with tabulated warm response](../../../pulse-tube-research/simulations/exploration/empirical_head/clamp_extension/results/report.json).
- [Falter et al. (2024), Cryogenics 144, 103976](https://doi.org/10.1016/j.cryogenics.2024.103976),
  Fig. 4: measured loaded temperature response; Fig. 7: simulated regenerator-2
  profiles at zero and finite clamp load.
