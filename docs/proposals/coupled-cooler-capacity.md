# Proposal: coupled pulse-tube capacity during warm cooldown

Status: exploratory v1, audited 2026-10-20

## Summary

The cooler link gives each stage a cooling power that depends on both stage temperatures, $Q_1(T_1, T_2)$ and $Q_2(T_1, T_2)$. This proposal records an exploratory six-constant functional form and the physical assumptions requiring further derivation. The target is the full cooldown from 300 K until the MXP reaches 40–50 K; the CFP can then be at 15–25 K. The candidate equations below serve as a phenomenological starting point for a coupled flow-and-energy model. Independent parameter estimates and verification of the temperature dependence are prerequisites for a priori prediction.

## Terms

- **Cold head**: the cold end of the pulse-tube cooler. Its two stages bolt to the PRP (first stage) and the CFP (second stage). $T_1$ and $T_2$ are the cold-head stage temperatures.
- **Regenerator**: a packed bed of fine metal that the helium flows back and forth through. It stores heat from the gas on one half of the cycle and returns it on the other. Regenerator 1 runs from room temperature down to stage 1. Regenerator 2 runs from stage 1 down to stage 2.
- **Pulse tube**: an open tube above each stage. The cooling happens at its cold end.
- **Phase shifter**: the valves and reservoir at the room-temperature end of each pulse tube. It sets the timing between pressure and flow.
- **Pressure swing**: the amplitude of the pressure oscillation the compressor and rotary valve impose. $p_0$ is the swing entering the cold head; $p_{c,i}$ is the swing at the cold end of stage $i$.
- **$T_\mathrm{amb}$**: room temperature at the warm end, about 295–300 K.

## Background

The simulated window runs from 300 K until the MXP turns on at 40–50 K. The manufacturer load map (`data/SN9_PTR_load_map.png`) covers $T_1 \approx 41$–$67$ K together with $T_2 \approx 2.2$–$4.7$ K. At the end of the window the first stage is at about 40–50 K and the second stage (CFP) is at about 15–25 K. The whole cooldown path therefore lies above the load map, and the capacity above the map sets the predicted cooldown time.

The brownfield uses one-variable curves $Q_i(T_i)$ from the literature (Wang 2016, Hao 2022, Green 2015). Those curves describe Cryomech PT415/420/425 coolers. The production cooler is a different design, the PTD 420Q.

Warm operation motivates an ideal-gas starting approximation. Its accuracy depends on helium pressure and temperature. Regenerator effectiveness also depends on matrix heat capacity, gas inventory, flow, and frequency. The target domain reaches at least 15 K on the CFP and potentially lower temperatures on its cold head across contact resistance; property and regenerator approximations require assessment over that domain.

## Physics in the window

### Gross cooling

For an ideal, lossless pulse tube, gross cooling equals the acoustic power, which is conserved along the tube. In a harmonic approximation,

$$
W_i = \tfrac12\,\mathrm{Re}(p_{c,i} U_{c,i}^{*})
    = \tfrac12\,\mathrm{Re}(p_{w,i} U_{w,i}^{*}).
$$

Here $p$ and $U$ are peak-amplitude pressure and volume-flow phasors, and $c$ and $w$ denote the cold and warm ends. With a purely resistive warm-end termination $R_{o,i}$ and negligible reservoir pressure oscillation, $W_i=|p_{w,i}|^2/(2R_{o,i})$. The flows and pressures follow from the oscillating-flow circuit. This is the ideal result in Radebaugh's enthalpy-flow treatment, §§5.1–5.2 [4].

The candidate form below retains the ansatz $W_i=\alpha_i T_i\pi_i^2$, where $\pi_i=|p_{c,i}|/|p_0|$, for discussion. Its explicit temperature dependence is an empirical assumption requiring a separate justification. Conservation of acoustic power supplies the physical starting equation above.

### Regenerator pressure drop

The gas is forced through each regenerator every cycle, and the viscous pressure drop reduces the swing that reaches the cold ends. In the viscous (Darcy) regime the drop is proportional to viscosity times volume flow. For helium between 20 K and 300 K, viscosity scales as roughly $T^{0.67}$. For a fixed mass flow, the volume flow scales as $T$. The drop therefore scales as $\bar T^{\,n}$ with $n \approx 1.7$, where $\bar T$ is the regenerator's mean temperature.

For the assumed series regenerator arrangement, the candidate attenuation law is

$$
\pi_1 = 1 - \kappa_1 \left(\frac{\bar T_1}{T_\mathrm{amb}}\right)^{n},
\qquad
\pi_2 = \pi_1 - \kappa_2 \left(\frac{\bar T_2}{T_\mathrm{amb}}\right)^{n},
$$

$$
\bar T_1 = \tfrac{1}{2}(T_\mathrm{amb} + T_1),
\qquad
\bar T_2 = \tfrac{1}{2}(T_1 + T_2).
$$

This term can bend the capacity curves over at high temperature. It also couples the stages: within this ansatz, a warm first stage raises the prescribed attenuation across both regenerators and lowers $Q_2$ when $\pi_2>0$. The arithmetic mean temperatures approximate the spatial temperature profiles.

The scalar subtraction assumes fixed flow-to-pressure ratios and effective phases. A self-consistent circuit determines flow and pressure together, using complex impedances and mass conservation. The fixed-mass-flow argument motivates a resistance scaling; applying that scaling directly to fractional pressure attenuation is an additional assumption. The candidate law requires $0\leq\pi_1\leq1$ and $0\leq\pi_2\leq1$ throughout its domain, because $\pi_i$ represents a magnitude ratio.

### Regenerator leakage

Passive axial conduction can be approximated by thermal conductances $G_i$ between the specified endpoint nodes:

- Regenerator 1 carries $G_1 (T_\mathrm{amb} - T_1)$ from room temperature into stage 1.
- Regenerator 2 carries $G_2 (T_1 - T_2)$ out of stage 1 and into stage 2.

These passive heat flows enter the endpoint balances with equal and opposite signs. Regenerator ineffectiveness contributes cycle-averaged enthalpy transport, while pressure drop dissipates acoustic power. Their endpoint contributions require a flow-and-energy balance. The candidate form currently represents gross cooling and passive conduction; its loss model remains to be completed. Conductance values approximate temperature-dependent material conductivity over the chosen domain.

## Candidate phenomenological form

$$
Q_1 = \alpha_1 \, T_1 \, \pi_1^2 \;-\; G_1\,(T_\mathrm{amb} - T_1) \;+\; G_2\,(T_1 - T_2)
$$

$$
Q_2 = \alpha_2 \, T_2 \, \pi_2^2 \;-\; G_2\,(T_1 - T_2)
$$

with $\pi_1$, $\pi_2$, $\bar T_1$, $\bar T_2$ as above and $n \approx 1.7$.

Properties of the form:

- At $T_1 = T_2 = T_\mathrm{amb}$ both leakage terms vanish, and each stage's capacity is its gross cooling reduced by the full room-temperature pressure drop.
- The stages are coupled in two ways. $Q_2$ depends on $T_1$ through regenerator 1's pressure drop and through $G_2$. $Q_1$ depends on $T_2$ through $G_2$.
- The $G_2$ term changes sign by itself when $T_2 > T_1$, which can happen early in a cooldown.

### Comparison with the ideal-regenerator form

An earlier linear candidate is

$$
Q_1 = a_1 T_1 - b_1 (T_\mathrm{amb} - T_1), \qquad Q_2 = a_2 T_2 - b_2 (T_1 - T_2).
$$

It corresponds algebraically to setting $\kappa_1 = \kappa_2 = 0$ and assigning the interstage conductance contribution only to $Q_2$. With the other temperature fixed, each capacity is linear in its own temperature. The brownfield curve `PT1.5 (1X F-70)` in `user_coolers.json`, described there as "Test Data Used", has 20 W at 50 K, 50 W at 150 K, and 59 W at 250 K. A straight line through the two ends falls 21% below the value at 150 K. The candidate attenuation term can supply curvature. Establishing its physical cause requires the curve's measurement provenance and accompanying first-stage temperatures or loads.

## Constants

| Constant | Units | Meaning | Set by |
|---|---|---|---|
| $\alpha_1$, $\alpha_2$ | W/K | gross-cooling scale in the candidate temperature ansatz | effective drive and circuit response |
| $\kappa_1$, $\kappa_2$ | — | fractional attenuation coefficients in the candidate law | effective geometry, packing, flow, and phase |
| $G_1$, $G_2$ | W/K | passive axial conductance between specified endpoints | material conductivity and geometry |
| $n$ | — | pressure-drop temperature exponent | helium viscosity, $n \approx 1.7$ |
| $T_\mathrm{amb}$ | K | warm-end temperature | the lab |

A possible parameter-sharing hypothesis assigns common $\kappa_i$ and $G_i$ to one design and per-unit $\alpha_i$. This assumes that unit variation enters mainly through the gross-cooling scales. The effective attenuation coefficients can also vary with flow and phase, so this hypothesis requires evidence. A priori use requires an independent source or physical estimate for each parameter.

### Optional characterization from measured cooldowns

A summed energy balance can characterize total cooling. Choose a boundary containing all modeled thermal masses and define $Q_\text{outside}$ as positive heat input and $Q_\text{cooler,total}$ as positive heat removal:

$$
Q_\text{cooler,total} = Q_\text{outside} - \frac{dE_\text{stored}}{dt},
\qquad
\frac{dE_\text{stored}}{dt}=\sum_j C_j(T_j)\,\frac{dT_j}{dt}.
$$

Internal heat transfers cancel when both endpoints and their associated storage lie inside the boundary. Cold-head storage, linked masses, and any additional stored energy belong in $E_\text{stored}$. If the cooler acts through the two cold ends and the RGP, its total removal is $Q_1+Q_2+Q_R$, using a consistent boundary for cooler-internal storage and heat flows.

For the two-output candidate above,

$$
Q_1+Q_2=\alpha_1T_1\pi_1^2+\alpha_2T_2\pi_2^2-G_1(T_\mathrm{amb}-T_1).
$$

Thus $G_2$ cancels exactly and the total balance evaluated on measured temperatures identifies combinations of the remaining parameters. Separate stage balances or other independent information are needed to constrain $G_2$. Separate balances include every residual heat transfer across their partition, including radiation and supports after gas links switch off.

Values inferred this way absorb errors in thermal masses, heat capacities, outside loads, and unobserved storage. Their independence for validation depends on the characterization runs and assumptions. Calibration remains an optional route when required.

## Place in the model

- The cooler is a link acting on the two cold-head stage nodes. It evaluates $Q_1$ and $Q_2$ from the current $T_1$ and $T_2$.
- The $G_2$ conductance is part of the cooler link, since it is cooler physics.
- Cold-head and regenerator storage are represented by thermal masses with explicit temperature assignments. A lumped representation approximates the regenerator's spatial temperature profile and its evolution.
- The joint between each cold-head stage and its plate is a separate link.

## Scope and extensions

The target domain includes cold-head temperature pairs visited from 300 K until MXP turn-on, including CFP temperatures around 15–25 K and potentially lower cold-head temperatures across joints. The manufacturer map supplies a separate low-temperature region. The operating domain and any connection between the regions require explicit treatment.

The candidate is a phenomenological ansatz motivated by pulse-tube physics. Its gross-cooling temperature dependence, attenuation law, and regenerator-loss accounting remain subjects for derivation and comparison with independent data.

Coupled effects for a revised derivation include:

1. **Temperature-dependent flow split.** Cold volumes hold more gas, so the division of flow between the two stages shifts as they cool.
2. **Compressor response.** The pressure swing the compressor delivers changes as the cold head's gas demand changes.

A lumped oscillating-flow network can include both. Regenerator and tube elements carry resistance, gas storage, and inertance where relevant, together with thermal energy accounting. Geometry, helium properties, drive conditions, and phase-shifter characteristics determine the elements. Solving the network at each temperature pair gives pressures and flows together. Its useful reduction and required parameters are subjects for v2.

## Open questions

1. **Regenerator plate (RGP).** The RGP is clamped to regenerator 2 between the stages, where the local temperature lies between $T_2$ and $T_1$. Heat it gives up enters regenerator 2. How does the RGP enter this form?
2. **Per-unit drive.** Is the drive one shared scale ($\alpha_1/\alpha_2$ fixed by design), or two independent per-unit values?
3. **Loss temperature dependence.** $G_i$ approximates passive conduction here. What separate model describes regenerator ineffectiveness, and what temperature dependence does each contribution require?

## References

- Radebaugh, R., "Development of the pulse tube refrigerator as an efficient and reliable cryocooler," *Proc. Institute of Refrigeration* 96 (2000).
- Swift, G. W., *Thermoacoustics: A Unifying Perspective for Some Engines and Refrigerators*, Acoustical Society of America (2002). Energy flux in regenerators and pulse tubes.
- Arp, V. D., McCarty, R. D., Friend, D. G., *Thermophysical Properties of Helium-4 from 0.8 to 1500 K with Pressures to 2000 MPa*, NIST Technical Note 1334 (1998). Helium viscosity.
- [4] Radebaugh, R., ["Pulse Tube Cryocoolers"](https://trc.nist.gov/cryogenics/Papers/Review/2003-Pulse_Tube_Cryocoolers.pdf), *Low Temperature and Cryogenic Refrigeration*, pp. 415–434 (2003), §§5.1–5.2. Ideal acoustic-power conservation and regenerator enthalpy flow.
- Radebaugh, R. et al., ["Inertance Tube Optimization for Pulse Tube Refrigerators"](https://trc.nist.gov/cryogenics/Papers/Pulse_Tube_Cryocoolers/2006-Inertance_Tube_Optimization.pdf), *Advances in Cryogenic Engineering* 51, pp. 59–67 (2006). Complex flow impedances and gas compliance.
