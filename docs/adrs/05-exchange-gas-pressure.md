# ADR-05: Exchange-gas pressure during cooldown

## Decision

Represent exchange gas as one link coupling the facing surface pairs in the
vacuum can at a shared pressure. For the reference refrigerator runs, use the
`monitor_ovc` pumping threshold as a proxy for that pressure. The threshold is a
cubic in CFP temperature, with its coefficients converted to SI units and supplied
in the link configuration. The pressure node is configurable.

Clamp the calculated pressure at zero. The reference threshold reaches zero near
68 K, ending the exchange-gas coupling continuously as the fridge cools. Use N2
to represent residual air and evaluate each pair's gas conductance at its mean
surface temperature, combining molecular and viscous resistances as for a gas-gap
heat switch.

## Basis and assumptions

The refrigerator routine pumps the vacuum can according to a temperature-dependent
pressure threshold. Using this threshold as the gas pressure assumes that
outgassing and residual air keep the measured P1 pressure close to it. The
reference cooldown records provide node temperatures; the corresponding P1
diagnostics have not been retrieved. The threshold therefore supplies an
operationally defined approximation for the thermal coupling.

The pressure is uniform throughout the can. Individual heat paths depend on the
configured facing areas and gaps and the temperatures of their surfaces. The
same pressure applies to all paths in the link.

## Consequences

If the actual pressure is below the pumping threshold, the link overstates the
gas coupling. The discrepancy scales approximately with pressure in the
molecular regime. Comparison with measured P1 diagnostics would test the pressure
assumption directly.

The configured coefficients and pressure node describe the particular routine
and system. Their values belong with the run inputs and their sources. The
pressure reaches zero through the threshold calculation, so a separate procedure
phase is unnecessary for that transition.

## References

- [System records](../../data/systems.md): cooldown records, operating routines,
  and available pressure diagnostics.
- [Exchange-gas implementation](../../src/pulse_tube_cooldown/library/links/exchange_gas.py):
  shared pressure and heat transfer across facing surface pairs.
- [Gas-gap conductance](../../src/pulse_tube_cooldown/library/links/gas_gap_switch.py):
  molecular and viscous gas conduction.
