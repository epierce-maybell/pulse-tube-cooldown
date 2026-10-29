# Decision: Gas-pushing model

Target: Derive a gas-pushing model from first principles, where the inputs are related to hardware and routines, rather than emergent properties of the gas-pushing cycle.

## Governing equations
Consider an ordered set of gas volumes, $j = 0, 1, ..., N$, from the warm manifold inward. Each element has a volume $V_j$, an average temperature $T_j$ of the nodes it touches (log-mean, to be precise), and we divide the warm side into an open volume $V_w$ and a buffer $V_b$. Then, the mole capacity of the ordered set is:

$$
A = \frac{V_w}{R T_a} + \sum_{j \geq 1} \frac{V_j}{R T_j},
$$

where $R$ is the molar gas constant and $T_a$ is the ambient air temperature. The still-side volume holds $n_s$ moles at pressure $P_s = n_s / A$ and the buffer holds $n_b$ moles at pressure $P_b = n_b R T_a / V_b$. We can also define the buffer's mole capacity as $B = V_b / (R T_a)$. The total inventory of gas is then $n = n_s + n_b$. 

During pressurization, we have:

$$
\frac{dn_b}{dt} = -\frac{dn_s}{dt} = \frac{S P_s}{R T_a},
$$

or written more simply:

$$
\frac{dP_s}{dt} = -\frac{P_s}{\tau},
$$

where $S$ is the volumetric pump speed and $\tau$ is the time scale of pressurization. Pressurization runs until one of two conditions is met: (1) $P_b - P_s = \Delta P_{wave}$ or (2) $P_b = P_{cap}$, whichever comes first. The wave strength $\Delta P_{wave}$ and buffer cap $P_{cap}$ are set by the gas-pushing routine. The durations of pressurization and depressurization are $t_p$ and $t_d$, respectively; $t_d$ is set by the routine, while $t_p$ is emergent from the pressurization equation and stop conditions. 

The maximum pressure for each cycle is $P_{max} = P_{min} + B \Delta P_{wave} / (A + B)$ and refill begins when the still-side pressure reaches a $P_{min} \leq P_{fill} - \Delta P_{refill}$, where $P_{fill}$ and $\Delta P_{refill}$ are set by the routine. Then, $t_p = \tau ln(P_{max} / P_{min})$, where the characteristic time scale $\tau = A R T_a / S$. If we average over one fill-refill cycle, we can pin maximum pressure to $P_{max} = P_{fill} - \text{max}(\Delta P_{refill} - \Delta P, 0) / 2$.

Let:

$$
A_j = \sum_{k \geq j} \frac{V_k}{R T_k}
$$

be the mole capacity of $j$ plus everything further along the ordered set than $j$, such that $A = A_0$. Each cycle, $\Delta P A_j$ moles cross element $j$, such that the heat flow into $j$ is:

$$
Q_j = \frac{5}{2} R f \Delta P \left[A_j(T_{j-1} - T_j) + A_{j+1}(T_{j+1} - T_j) \right],
$$

where $\Delta P = P_{max} - P_{min}$, the frequency $f = 1 / (t_p + t_d)$ and $A_{N+1} = 0$. The initial temperature $T_0 = T_a$. Each element distributes $Q$ evenly between adjacent nodes.

## Assumptions and simplifications
1. Gas equilibrates completely in each volume that it enters.
2. Helium is approximated as an ideal gas at the temperatures we care about.
3. During a cycle, pressure is uniform across the set, and the two sides equilibrate while the valve is open.
4. Each element is a link that distributes heat evenly to its adjacent nodes.
5. The pump operates at a constant speed.
