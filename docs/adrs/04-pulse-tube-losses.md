# ADR-04: Gas-wall thermal-shuttling loss

## Observed behavior
Our load maps exhibit asymmetric coupling between the two pulse tube stages. Let $Q_1(T_1, T_2), Q_2(T_1, T_2)$ represent a coupled capacity map. In our simulations, when $Q_1$ increases for fixed $Q_2$, $T_1$ increases and $T_2$ stays relatively constant. But, when $Q_2$ increases for fixed $Q_1$, both $T_1$ and $T_2$ increase. In terms of the partial responses:

$$
\frac{\partial Q_1}{\partial T_2} < 0 \qquad
\frac{\partial Q_2}{\partial T_1} \approx 0.
$$

Define the Jacobian $J$ as:

$$
\begin{bmatrix}
\frac{\partial Q_1}{\partial T_1} & \frac{\partial Q_1}{\partial T_2} \\
\frac{\partial Q_2}{\partial T_1} & \frac{\partial Q_2}{\partial T_2}
\end{bmatrix},
$$

or more simply:

$$
\begin{bmatrix}
a & b \\
c & d
\end{bmatrix}.
$$

If we invert this and look at the determinant, we see that:

$$
\frac{\partial T_1}{\partial Q_2} \approx -b / (ad).
$$

To satisfy our observations, we then have the following three conditions:
1. The partial derivative of $Q_1$ with respect to $T_2$ must be negative: $b < 0$. 
2. The second partial of $b$ with respect to $T_1$ must be positive: $\frac{\partial b}{\partial T_1} > 0$. 
3. The partial derivative of $Q_2$ with respect to $T_1$ must be approximately zero: $c \approx 0$.

Let $Q_1 = C_1 - L_1$, where $C_1$ is everything that contributes to cooling power, and $L_1$ is everything that contributes to losses. Then:

$$
b = \frac{\partial Q_1}{\partial T_2} = \frac{\partial C_1}{\partial T_2} - \frac{\partial L_1}{\partial T_2},
$$

where

$$
\frac{\partial L_1}{\partial T_2} > \frac{\partial C_1}{\partial T_2}
$$

ensures that requirement 1, above, is satisfied. In plain language: warming stage 2 needs to increase the loss rate in stage 1 more than the cooling rate, while still ensuring positive cooling above the stage 2 base temperature.

## Hypothesis
`Here, I hypothesize that stronger oscillatory flow increases the length scale of heat exchange and the exchange rate between gas and the pulse tube wall; in a particular regime, these losses can grow faster with temperature than useful cooling.`

In a simplified, one-dimensional model, let $x$ increase from the cold end of the tube to the warm end. The gradient $G = \partial T \partial x > 0$. For oscillatory flow, the gas velocity $u(t) = U \, cos(\omega t)$ for velocity amplitude $U$. Let $\theta$ denote some perturbation in temperature from gradient-driven effects only (neglecting the compression-driven temperature repsonse), and then we can write:

$$
\frac{\partial \theta}{\partial t} + a\theta = -Gu, 
$$

for some gas-wall relaxation rate $a$. The average quantities are therefore:

$$
\bar{(u\theta)} = -\frac{aGU^2}{2(a^2 + \omega^2)},
$$

and the loss $L$ is:

$$
L = \frac{A\rho c_p a G U^2}{2(a^2 + \omega^2)},
$$ 

for cross-sectional area $A$ and helium density $\rho$ and specific heat $c_p$. 

Consider then the log-change in loss:

$$
\frac{\partial ln(L_1)}{\partial ln(U)} = 2 + n\frac{\omega^2 - a^2}{\omega^2 + a^2}, 
$$

where we define

$$
n = \frac{\partial ln(a)}{\partial ln(U)}.
$$

This leads to an interesting regime where $a << \omega$. Then:

$$
\frac{\partial ln(L_1)}{\partial ln(U)} \approx 2 + n,
$$

and so $L_1 \propto U^{2 + n}$, while we from acoustic power theory, we approximate to first order that $C_1 \propto U^2$ locally.

For a small relative flow increase, $dU/U$, we can approximate $dQ_1 \approx \left[2 C_1 - (2 + n) L_1 \right] dU/U$. There is a velocity coupling coefficient:

$$
\epsilon = \frac{1}{U}\frac{\partial U}{\partial T_2},
$$

where $b = \epsilon \left[2 C_1 - (2 + n) L_1 \right]$. 

From here, we have all the necessary pieces for a first reasonable test of our hypothesis. To satisfy requirement 1, that $b < 0$, we need $L_1 / C_1 > 2 / (2 + n)$ while positive cooling demands that $L_1 / C_1 < 1$. To satisfy requirement 2, we can look at $G$ over the length $l$ of tube 1, so $G \approx (T_0 - T_1) / l$. When $T_1$ increases, if we hold all other local state variables constant, the gradient along the tube $G$ decreases, which in turn causes the loss $L_1$ to decrease. For $b < 0$, decreased loss pushes $b$ back towards zero, as $\partial b / \partial T_1 > 0$. To satisfy requirement 3, we need the coupling to act primarily through the stage 2 velocity amplitude, $U_2$, such that $c \approx (\partial Q_2 / \partial U_2) (\partial U_2 / \partial T_1)$. Then, a sufficiently small $\partial U_2 / \partial T_1$ gives $c \approx 0$, even while the converse response $\partial U_1 / \partial T_2$ is appreciable.
