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

Let $Q_1 = C_1 + L_1$, where $C_1$ is everything that contributes to cooling power, and $L_1$ is everything that contributes to losses. Then:

$$
b = \frac{\partial Q_1}{\partial T_2} = \frac{\partial C_1}{\partial T_2} - \frac{\partial L_1}{\partial T_2},
$$

where

$$
\frac{\partial L_1}{\partial T_2} > \frac{\partial C_1}{\partial T_2}
$$

ensures that requirement 1, above, is satisfied. In plain language: warming stage 2 needs to increase the loss rate in stage 1 more than the cooling rate, while still ensuring positive cooling above the stage 2 base temperature.
