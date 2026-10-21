# Design Specification: Pulse Tube Cooldown Simulation

## Scope
Simulate temperature and heat fluxes in a Maybell dilution refrigerator during cooldown, roughly from 300 K to when the MXP turns on (40-50 K). Version 1.0.0 targets a full cooldown simulation of a 6-plate system (SN5), with error comparable to or lower than the legacy model. For this task, we'll define "error" more precisely as the least-squares residual to log-temperature across all thermometers mounted on the system (except for MXP). 

## Public API
For now, we'll call the library Pulse Tube Cooldown, and we can have a canonical import statement:
```
import pulse_tube_cooldown as ptc
```

#### Configuration files
Each system gets one, self-contained configuration file. Each model run takes exactly one configuration file as input.
Configuration files are written in TOML, with the following sections:
```toml
[fridge]

[simulation]

```

#### Public functions
The library exposes the following public functions:
- configure_model(path_to_config) --> Model
- simulate_cooldown(model) --> Result
- save_output(result, path_to_output) --> None

So the canonical run is:
```python
model = ptc.configure_model(path_to_config)
result = ptc.simulate_cooldown(model)
ptc.save_output(result, path_to_output)
```

## Data Model

## Engine

## Dilution Refrigerator Library

## Validation and Testing
