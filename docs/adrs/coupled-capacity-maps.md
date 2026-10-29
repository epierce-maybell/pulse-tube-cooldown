# Decision: How to model pulse tube cooling power?

Target: Derive a model for pulse tube cooling power from first principles and what we can measure for each pulse tube unit, accounting for a regenerator, where each of Q1, Q2, and QR are appropriately coupled. 

## Requirements
1. We want a model that characterizes pulse tubes per unit - there's variability, and our model should capture that;
2. For post-tuning data, we have continuous cooldown data on a test stand and load maps at cold temperatures only (including TR, but only at QR = 0);
3. We may ask for targeted experiments iff we can show they are high-leverage, but we cannot ask for more data for every pulse tube unit;
4. As much as possible, we don't want to hide behind calibration: a model that fits one fridge well but doesn't generalize is not useful, and if a model fails to fit a fridge, we need to be able to attribute that failure to a particular system;
5. There is a regenerator clamp between stages 1 and 2, so the true map is of the form Q1(T1, T2, TR), Q2(T1, T2, TR), QR(T1, T2, TR). Our model should get as close to that fidelity as is reasonable, but it may not be possible to do so with the data we have;
6. Fridge-plus units have two pulse tubes, so we should make sure our model can handle interactions between them (this probably doesn't change most designs, but worth keeping an eye on).

## Governing equations
Stage 1 is the primary plate (PRP), stage 2 is the cold-foot plate (CFP), and there is a regenerator plate (RGP) clamped at some fraction $x_R \in [0, 1]$ of the distance between them; quantities on these are denoted $\Box_1$, $\Box_2$, and $\Box_R$, respectively. The temperature of each stage, $T_i$, is held fixed over one valve period, so that the energy flows (gas enthalpy plus axial conduction) are averaged across cycles. The energy balances for each stage are:

$$
Q_1 = H_1 - L_1  - K_1 (T_a - T_1) + L_R, \qquad
Q_R = L_2 - L_R, \qquad
Q_2 = H_2 - L_2 - K_2 (T_a - T_2),
$$

where $Q_i$ is heat flow into each stage, $H_i$ is the enthalpy flow away from each stage, $K_i$ is the thermal conductance of the tube ending at stage $i$, and $L_i$ is the loss down the regenerator ending at plate $i$, such that:

$$
L_i = (1 - \epsilon_i) \dot{m}_i c_p \Delta T_i,
$$

where $\dot{m}_i$ is the mean mass flow through the section, $c_p$ is the specific heat of helium, $\Delta T_i$ is the temperature difference across the section (e.g., $\Delta T_1 = T_a - T_1$), and $\epsilon_i = 1 - a_i u_i$ is the effectiveness in terms of a packing constant $a_i$ of order one and utilization $u_i$, defined by:

$$
u_i = \frac{\dot{m}_i c_p}{f M_i c_{m,i}(\bar{T_i})},
$$

where $f$ is the valve frequency, $M_i$ is the mass matrix of the section, $c_{m,i}$ is the specific heat of the mass matrix, and $\bar{T_i}$ is the mean temperature of the section. This is essentially the heat capacity of the gas moving through a section per cycle divided by the heat capacity of the matrix it passes over. At $T_2$, helium is far from an ideal gas, so we add a correction to of the form $L_2 += (1 - \beta(T_2, \bar{P}) T_2) H_2$, where $\beta$ is helium's isobaric expansion coefficient. 

One further correction is required for the regenerator clamp, which sits at some normalized position $x_R \in [0, 1]$ between stages 1 and 2. The two sections on either side of $x_R$ are one matrix of mass $M_{col}$ divided into $M_R = x_R M_{col}$ and $M_2 = (1 - x_R) M_{col}$. 

Then, mass flows are given by:

$$
\dot{m}_i = f \sum_{j \in D_i} V_j \left[\rho(T_j, P_h) - \rho(T_j, P_l) \right],
$$

where $D_i$ is the set of gas-filled components including $i$ and everything colder than it, $P_h$ and $P_l$ are the supply and return pressures, respectively, and $\rho$ is the density of helium. The mass flow of the entire system is set by the drive constraint:

$$
\dot{m}_1 = \dot{m}_c(P_h, P_l) = \dot{V}_c \rho(T_a, P_l),
$$

where $\dot{m}_c$ is the drive's flow characteristic and $\dot{V}_c$ is the compressor's volumetric displacement rate. The mean pressure is fixed by the gas charge, so $\bar{P} = (P_h + P_l) / 2$ always, and the swing in pressure $\Delta P = P_h - P_l$ is set by the drive constraint, subject to $\Delta P \leq \Delta P_{max}$, falling as the head cools. Where two pulse tube heads share one compressor, $D_1$ contains both heads' components and they interact through the pressure wave. The pressure amplitude is then $P_a = \frac{1}{2} \kappa \Delta P$, where $\kappa$ is a waveform factor of order one that accounts for non-sinusoidal behavior.

From this, we get the cycle-averaged enthalpy flow in each stage:

$$
H_i = \frac{1}{2} C_i P_a^2,
$$

where $C_i$ is an effective conductance for each valve.

## Data requirements
Each pulse tube unit needs $C_1$, $C_2$, $a_1$, and $a_2$ from post-tuning load map data.

Drive values $\dot{V}_c, \kappa, $\Delta P_{max}$, $\bar{P}$, and $f$ can be obtained once, as they should not vary much between compressors. 

The values $V_j, x_R, M_1 and M_{col}, c_{m,i}, and K_i depend on the geometry and materials of each fridge design.

## Assumptions and simplifications
1. Helium properties come from he34pack.
2. Enthalpy flow is evaluated at the cold end of each element.
3. Temperatures are fixed over one cycle.
4. There is no pressure drop along the regenerator column.
5. Each regenerator section has a linear temperature profile.
6. $C_i$ is constant because the valve sits at $T_a$.
7. $a_i$ is constant because the temperature dependence is carried entirely by $c_{m,i}$.
8. The response of $Q_R$ to $T_R$ is captured by $a_2$ and is only tested by the RGP trace.
