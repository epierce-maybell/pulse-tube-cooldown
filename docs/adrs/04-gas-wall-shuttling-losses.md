# ADR-04: Gas-wall thermal-shuttling loss

## Observed behavior
Our load maps exhibit asymmetric coupling between the two pulse tube stages. Let $Q_1(T_1, T_2), Q_2(T_1, T_2)$ represent a coupled capacity map. (Here, $Q_i$ denotes positive cooling capacity; this is opposite from the convention in ADR-03.) In our measurements, when $Q_1$ increases for fixed $Q_2$, $T_1$ increases and $T_2$ stays relatively constant. But, when $Q_2$ increases for fixed $Q_1$, both $T_1$ and $T_2$ increase. In terms of the partial responses:

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

If we invert this, holding $c \approx 0$, we see that:

$$
\frac{\partial T_1}{\partial Q_2} = -b / (ad - bc) \approx -b / (ad). 
$$

To satisfy our observations, we then have the following three conditions:
1. The partial derivative of $Q_1$ with respect to $T_2$ must be negative: $b < 0$.
2. The partial derivative of $b$ with respect to $T_1$ must be positive: $\frac{\partial b}{\partial T_1} > 0$.   
3. The partial derivative of $Q_2$ with respect to $T_1$ must be approximately zero: $c \approx 0$.

Let $Q_1 = C_1 - L_1$, where $L_1$ is the total gas-transport loss, including shuttling and breakthrough, and $C_1$ is the cooling capacity not accounting for these losses. Then:

$$
b = \frac{\partial Q_1}{\partial T_2} = \frac{\partial C_1}{\partial T_2} - \frac{\partial L_1}{\partial T_2},
$$ 

where

$$
\frac{\partial L_1}{\partial T_2} > \frac{\partial C_1}{\partial T_2}
$$

ensures that requirement 1, above, is satisfied. Now, $C_1$ includes the regenerator flow that stage 1 sends down the stage-2 column. At steady state with no clamp load, that flow is uniform along the column and arrives at stage 2. Warming stage 2 lowers the flow from stage 1, which costs stage 1 about one watt of capacity per watt gained at stage 2. So we can expect that $\frac{\partial C_1}{\partial T_2} < 0$ holds. 

## Hypothesis
`Gas-wall shuttling transports heat towards the cold end as oscillating gas exchanges heat with the pulse tube wall. We hypothesize that warming stage 2 changes the shared pressures and stage 1 gas motion enough to increase this shuttling loss more than the stage 1 cooling power increases.`

In a simplified, one-dimensional model, let $x$ increase from the cold end of the tube to the warm end. The gradient $G = \partial T / \partial x > 0$. For oscillatory flow, the gas velocity $u(t) = U cos(\omega t)$ for velocity amplitude $U$. Let $\theta$ denote some perturbation in temperature from gradient-driven effects only (neglecting the compression-driven temperature response), and then we can write:

$$
\frac{\partial \theta}{\partial t} + \alpha\theta = -Gu,
$$ 

for some gas-wall relaxation rate $\alpha$. The average quantities are therefore:

$$
\bar{(u\theta)} = -\frac{\alpha G U^2}{2(\alpha^2 + \omega^2)},
$$

and the shuttling loss $L_s$ is:

$$
L_s = \frac{A\rho c_p \alpha G U^2}{2(\alpha^2 + \omega^2)},
$$ 

for cross-sectional area $A$ and helium density $\rho$ and specific heat $c_p$. 

At a fixed temperature profile and mean pressure $\bar{P}$, laminar shuttling grows as $U^2$. Warming stage 2 also changes mean pressure and pressure amplitude $P_a$, so the loss response is:

$$
\frac{\partial L_1}{\partial T_2} = 
(\frac{\partial L_1}{\partial P_a})(\frac{\partial P_a}{\partial T_2}) +
(\frac{\partial L_1}{\partial \bar{P}})(\frac{\partial \bar{P}}{\partial T_2}),
$$

where we note that the pressure derivatives include their effects on velocity and helium properties. The complete stage balance is then required to determine the cooling response. 

## Specific parameters
The key parameter in this model is $\alpha$, which in a cylindrical tube is modeled as:

$$
\alpha = \frac{2h_g}{\rho c_p r},
$$

where $h_g$ is the gas-wall heat-transfer coefficient and $r$ is the inner radius of the cylinder. If flow is laminar, the heat-transfer coefficient is:

$$
h_g = \frac{k}{\delta_k}, \qquad
\delta_k = \sqrt{2k / (\rho c_p \omega)},
$$ 

for thermal conductivity $k$ and other parameters as defined above. The angular frequency is $\omega = 2 \pi f$, where $f$ is the set valve frequency.

## Decision
Given the temperature readings from pulse tube walls and the physics of the system, we should not ignore gas-wall heat exchange in the pulse tube model. Implementation of this ADR is straightforward: each pulse tube $i$ carries a shuttling loss $L_{s,i}$ into its stage, evaluated at the cold end of the tube with $G_i = (T_0 - T_i) / l_i$, where $l_i$ is the tube length. It depends on physical properties of helium, evaluated at $T_i$ and mean pressure, the tube's geometry, and a cold-end velocity amplitude $U_i$.
