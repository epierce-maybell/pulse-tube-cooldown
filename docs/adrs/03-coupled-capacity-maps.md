# ADR-03: Pulse tube cooling power

Target: Derive a model for pulse tube cooling power from first principles and what we can measure for each pulse tube unit, accounting for a regenerator, where each of Q1, Q2, and QR are appropriately coupled. 

## Requirements
1. The model should be able to predict cooling power from the same unit, on the same drive, for a different trajectory in (T1, T2, TR) space. It should also be able to predict cooling power under a different drive configuration, if the number of compressors or heads on a manifold changes.
2. Per-unit parameters should be fit from a single load map; in future coupled models, we can require a load map as input for the pulse tube head. In an ideal case, cooldowns on the test stand would be used for validation only, but if certain population-level parameters are required, we can consider fitting those once per design or group.
3. We are not concerned with run-to-run variability among units. If we decide that this can be *useful* to a model later on (e.g., testing behavior that should be invariant), we can revisit this. But we have no details about what changed between runs, and the list of possibilities is too long to constrain anything.
4. Coupling between stages should emerge from model physics, rather than being fit directly. With heat flow into a node positive, our data imply that dQ1/dT2 should be positive at fixed T1, about +1.4 W/K on one compressor and +2.7 W/K on two, weakening as T1 rises; dQ2/dT1 at fixed T2 is nearly zero.
5. QR is coupled to Q1 and Q2 in both directions, and its behavior follows from model physics and design parameters. Test stand load maps have nothing mounted there, so they are a point at QR = 0.
6. A prediction will be judged on the cooldown trajectory for T1 and T2 and cooling powers Q1 and Q2 at base.
7. This routine will run inside an implicit solver of some kind. The heat flows need to be functions of temperature alone with continuous derivatives. Internal solves are acceptable, but need to be weighted against the cost of nested iteration.
8. We can use our data, plus published data as a test of functional forms or other invariants (e.g., Snodgrass 2023, 2024). Do not propose targeted experiments or additional per-unit data collection.
9. Fitted parameters should be related to physical properties of the system, rather than acting as lumped or black-box coefficients.

## Governing equations
We model each cold head as one regenerator column, split into two sections (one from the warm flange to stage 1, and the other from stage 1 to stage 2), and two pulse tubes, connected to a shared gas supply. The regenerator clamp, holding the regenerator plate (RGP), is connected partway up the second (stage 1 to 2) regenerator tube. Cooling at each stage follows from the energy balance at its heat exchanger. The regenerator flows connect those balances, and the gas demand across the entire manifold sets the shared pressure swing.

We hold temperatures fixed over each valve cycle and calculate cycle-averaged powers. Acoustic power, $E$, describes the energy carried by the oscillating pressure and gas flow. Total energy flow, $F$, includes the energy transported through the regenerator by the gas and axial heat conduction. By convention, here, heat flow into a node is positive and similarly, $F$ is positive towards the cold end of the regenerator.

We write $F^+$ and $F^-$ for $F$ at the warm and cold ends of a regenerator section, respectively. At stage 2, the pulse tube removes acoustic power $E_2$, while the regenerator delivers energy flow $F^-_2$ and the tube wall conducts heat from the warm flange. At stage 1, the pulse tube removes $E_1$, the warm regenerator section delivers $F^-_1$, and the cold section carries $F^+_2$ away towards stage 2. The heat flows are then:

$$
Q_1 = F^-_1 + L_1 + K_1 (T_0 - T_1) - E_1 - F^+_2, \qquad
Q_2 = F^-_2 + L_2 + K_2 (T_0 - T_2) - E_2,
$$

Here, $Q_1$ and $Q_2$ are the cooler's heat flows into the stage nodes, in W. $K_1$ and $K_2$ are the pulse-tube wall conductances, in W/K, $T_1$ and $T_2$ are the stage temperatures, and $T_0$ is the warm-flange temperature. $L_1$ and $L_2$ are the total loss in each pulse tube. As of this ADR, we only consider breakthrough loss. 

Gas that crosses the whole pulse tube in one cycle, "breaking through" the buffer, carries its heat directly to the stage ("large-stroke loss" in Snodgrass et al., 2024). Each cycle, the mass that flows through the warm end of pulse tube $i$, into and out of its reservoir, is $\Delta m_i = V_{r,i} \Delta\rho_{r,i}$. Gas entering the tube's warm end must push the tube's gas out through its cold end ahead of it. The tube's own compression moves gas in and out of the cold end without carrying warm gas along, so warm gas reaches the stage only when $\Delta m_i$ exceeds the tube's buffer of gas, $\lambda m_i$. For a tube at one temperature, $\Delta m_i / m_i$ is the stroke $2|\xi|/L$ of Snodgrass et al. (2024) eq. (4), taken from the flow at the warm end, so $\lambda$ and $\sigma$ read from their Fig. 8 apply directly. The excess breaks through each cycle, so the loss is:

$$
L_i = f \left[h(T_0, \bar{P}) - h(T_i, \bar{P})\right] s(\Delta m_i - \lambda m_i),
$$

where $h$ is helium's specific enthalpy, $m_i$ is the tube's gas mass at $\bar{P}$, $s$ is a smooth positive part whose width, $\sigma m_i$, is the spread in how far gas parcels travel along the tube each cycle, and $\lambda$ is the fraction of the tube's gas that acts as a buffer. $\sigma$ and $\lambda$ each take one value per design. $\sigma$ comes from published data or an estimate and is not fitted. Snodgrass et al. (2024) eq. (4) motivates the large-stroke criterion; here we approximate it using a mass buffer and resulting enthalpy loss. 

Regenerator energy flow then combines pressure-driven transport and a contribution from temperature gradients, including imperfect gas-matrix heat exchange and axial conduction. At steady state, $F$ is constant between external heat inputs, but during cooldown, differences in $F$ cool or heat the matrix. For one-dimensional $x$ increasing from the cold end towards the warm end of each section, we have:

$$
F = F_{pressure} + F_{regeneration} + F_{conduction},
$$

where:

$$
F_{pressure} = \left[1 - \frac{\beta T}{1 + \epsilon} \right] E
$$

for helium's isobaric expansion coefficient $\beta$ and $\epsilon = C'_g / C'_m$ for gas $C'_g$ and matrix $C'_m$ heat capacity per length, in J/(m K) (see Snodgrass, 2023); and:

$$
F_{regeneration} = a \frac{(\dot{m} c_p)^2}{f C'_m} \frac{dT}{dx},
$$

for a dimensionless packing coefficient $a$, helium's specific heat $c_p$, the mass-transfer rate $\dot{m}$, in kg/s, and the valve frequency $f$. This is an approximation that basically details how much energy is carried towards the cold end of the regenerator, given imperfect exchange of heat between the gas and the matrix. Lastly:

$$
F_{conduction} = k_{ax} A \frac{dT}{dx},
$$ 

where $k_{ax}$ is the thermal conductivity in the axial direction and $A$ is the cross-sectional area. The product $k_{ax} A$ is the regenerator tube wall's conduction; conduction through the matrix is neglected.

Acoustic power varies along the regenerator as:

$$
\frac{d\ln E}{dx} = \frac{\beta}{1 + \epsilon} \frac{dT}{dx},
$$

using the low-frequency approximation of Snodgrass et al. (2023, Appendix C), which neglects viscous and thermal relaxation dissipation in the acoustic-power balance. The stage 2 section has cold-end acoustic power $E_2$. The stage-1 section has cold-end acoustic power $E_1 + E^+_2$, where $E^+_2$ is the acoustic power at the warm end of the stage-2 section. 

Each cycle, the mass entering head $n$ is:

$$
\Delta M_n = \sum_j V_j \, \Delta\rho_j
$$

summed over its gas volumes $j$, each with volume $V_j$ and temperature $T_j$, where $\Delta \rho_j$ is the change in density between the high (entry) pressure $P_+$ and low (return) pressure $P_-$. Pore gas exchanges heat with the matrix, so its density changes at its node's temperature. The other gas is compressed adiabatically, so its density difference is taken at constant entropy: the gas in each pulse tube, along a linear temperature gradient from the warm flange to its stage; the volume between the rotary valve and the head, at $T_0$; and each reservoir, also at $T_0$, whose pressure swing is smaller by a factor of $\sqrt{1 + (\omega R_i B_i)^2}$ (see below). 

The compressor(s) must provide that mass at each head's valve frequency $f_n$, so:

$$
\sum_n f_n \Delta M_n = \rho(T_0, P_-) \sum_k \dot{V}_k,
$$ 

where $`\dot{V}_k`$ is the volume per second drawn into compressor $k$ at its low-pressure inlet. This supply equation sets the swing $\Delta P = P_+ - P_-$ up to $\Delta P_r$, the relief point set by the compressor (or lowest among multiple compressors). Above that, the relief valve returns the excess flow and holds the swing at $\Delta P_r$, with a smooth transition between the two.

The entire manifold's helium charge $M_{He}$, in kg, sets the mean pressure, $\bar{P} = (P_+ + P_-) / 2$, where:

$$
M_{He} = \sum_j V_j \rho(T_j, \bar{P}),
$$

summed over all gas volumes $j$ on the manifold, each once: supply and return lines, buffers, head dead volumes, pore gas, pulse tubes, and reservoirs. Lines and buffers enter the inventory but not $\Delta M_n$. All pressures are absolute. The pressure amplitude is given by:

$$
P_a = \frac{1}{2} \kappa \Delta P,
$$

where $\kappa$ accounts for the shape of the wave (1 for a sinusoidal wave); we follow Snodgrass et al. (2024) and use $\kappa = \sqrt{2}$, corresponding to a square wave, and any waveform error beyond that ends up getting captured in the resistance (below). 

Then, the acoustic power available to each pulse tube $i$ is:

$$
E_i = e_i P_a^2 
$$

for some coefficient $e_i$ that depends on valve resistance, reservoir volume, and valve frequency. Each pulse tube connects to a reservoir through a needle valve which restricts gas flow; we treat the valve resistance $R_i$ as constant and the reservoir as an ideal gas held at $T_0$. Then:

$$ 
e_i = \frac{\omega^2 R_i B_i^2}{2 \left[1 + (\omega R_i B_i)^2 \right]},
$$

where $\omega = 2\pi f$ and the pressure response $B_i = Vr_i / (\gamma \bar{P})$ depends on reservoir volume $Vr_i$ and helium's heat capacity ratio $\gamma$. Here, $B_i$ is the gas volume entering the reservoir per unit pressure increase, in m^3/Pa. We use the resulting power as the pulse tube's cold-end acoustic power; the tube's heat losses enter the stage balance as $L_i$. Note that when $\omega R_i B_i >> 1$, the reservoir pressure changes very little, and the expression for acoustic power reduces nicely to $E_i = P_a^2 / (2 R_i)$.

## Regenerator network
Each regenerator section is divided along its length into matrix nodes. Node $k$ has matrix mass $M_k$ with specific heat $c_{m,k}$, pore volume $V_k$, and temperature $T_k$ that the solver integrates. We assume that the pore gas shares that temperature. The capacity ratio of each node is:

$$
\epsilon_k = \frac{V_k \rho(T_k, \bar{P}) c_p(T_k, \bar{P})}{M_k c_{m,k}(T_k)},
$$ 

and its heat capacity is $(1 + \epsilon_k) M_k c_{m,k}$. At each boundary between nodes, the energy flow and acoustic power relations hold, including at the boundary where a node meets a plate. The mass flow at a boundary is $f$ times the mass swing $V_j \Delta\rho_j$ of every gas volume on its cold side, so it is $f_n \Delta M_n$ at the warm flange. At the stage-1 end of the stage-1 section, those volumes are tube 1, reservoir 1, the whole stage-2 section, tube 2, and reservoir 2. The heat flow into a node is $F^+ - F^-$, the energy flows at its warm and cold boundaries. $`F^-_1`$, $`F^+_2`$, and $`F^-_2`$ are $F$ at the stage boundaries, and $Q_R$ is the heat flow into the node containing the clamp, which holds the RGP. At the stage-2 boundary, $F_{pressure}$ is evaluated at the plate temperature $T_2$ with the coldest node's $\epsilon$. Where the temperature gradient there is small, stage 2's net cooling is $\beta T_2 E_2 / (1 + \epsilon) - K_2 (T_0 - T_2) - L_2$.

## Parameters
Four parameters are fitted to the load map, per head: packing coefficients $a_i$ and valve resistances $R_i$. 

For the pulse tube head, each design specifies the geometry of the regenerator sections and pulse tubes, the gas volumes, clamp position, and regenerator matrix materials and properties. Each unit's valve frequency is its shipped setting. For each drive, we need the compressor volumetric flow rates, relief swing, fill pressure at rest, warm-line volume, and the configuration of compressors and heads on the manifold. 

Model parameters that we don't fit per run include $\kappa$, $\lambda$, and $\sigma$, all taken from Snodgrass et al. (2024) for now.

To start a run, we need the warm-flange temperature, $T_0$. 

## Decision
Accept. This is our best approximation of a pulse tube cooling power model, with the information we have available at the time of writing.
