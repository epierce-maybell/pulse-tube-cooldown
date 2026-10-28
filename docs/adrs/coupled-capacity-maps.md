# Decision: How to model pulse tube cooling power?

## Requirements
1. We want a model that characterizes pulse tubes per unit - there's variability, and our model should capture that;
2. For post-tuning data, we have continuous cooldown data on a test stand and load maps at cold temperatures only;
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
L_1 = \Lambda_1 (T_a - T_1), \qquad
L_R = \Lambda_R (T_1 - T_R), \qquad
L_2 = \Lambda_2 (T_R - T_2) + (1 - \beta(T_2, P_0) T_2) H_2,
$$

for the effective conductances $\Lambda_i$, ambient temperature $T_a$, and helium's isobaric expansion coefficient $\beta$. The comparable terms $\beta_i$ evaluated at $T_1$ and $T_R$ are dropped, since helium at those stage's temperatures behaves within a few percent of an ideal gas. Enthalpy flow is given by:

$$
H_i = \frac{1}{2} \, Z_i \, R \, T_i \, m_i \, \frac{P_a}{P_0} \cos \phi_i,
$$

where $Z_i$ is helium's compressibility factor at $(T_i, P_0)$, $R$ is the helium gas constant, $m_i$ is the mass-flow amplitude, and $\phi_i$ is the phase lead. The mean pressure $P_0$ and the amplitude $P_a$ are set by the supply and return pressures, $P_h$ and $P_l$, such that:

$$
P_0 = \frac{1}{2}(P_h + P_l), \qquad
P_a = \kappa \frac{1}{2}(P_h - P_l),
$$

where $\kappa$ is a waveform factor of order one. The pressure swing falls as the gas cools and the downstream demand rises. This leads to a drive constraint, where:

$$
\Delta m = \sum_j V_j \left[ \rho(T_j, P_h) - \rho(T_j, P_l) \right], \qquad
f \, \Delta m = \dot m_c (P_h, P_l),
$$

where $\rho$ is the density of helium, $\Delta m$ is the mass exchange between the compressor and the head per cycle, related to the sum of gas volume $V$ and temperature $T$ of all gas-filled components $j$, $P_h$ and $P_l$ are the supply and return pressures, respectively, $f$ is $1 /$ the period, and $\dot m_c$ is the drive's flow characteristic.

## Assumptions and simplifications
1. Helium properties come from he34pack. We concern ourselves with the temperature-dependence, under the assumption that the mean pressure does not change much.
2. Enthalpy flow is evaluated at the cold end.
3. Temperatures are fixed over a cycle.
4. There is no pressure drop along the regenerator column.
