# Proposal: Toy run configuration

A complete configuration for the first-slice toy run (task 7.4). It exercises conduction, a fixed boundary node, and a phase switch, using the activation and run-input decisions from tasks 5.1 and 5.2.

```toml
[meta]
experiment = "Toy conduction with a phase switch"

[[fridge.nodes]]
name = "Plate"
composition = [["OFHC Cu RRR100", 1.0]]

[[fridge.nodes]]
name = "Bath"
composition = [["OFHC Cu RRR100", 1.0]]

[[fridge.links]]
name = "Rod A"
type = "conduction"
nodes = ["Plate", "Bath"]
material = "304 SS"
area = 1.0e-4
length = 0.1

[[fridge.links]]
name = "Rod B"
type = "conduction"
nodes = ["Plate", "Bath"]
material = "304 SS"
area = 1.0e-4
length = 0.1

[[procedure.phases]]
name = "single rod"
disable = ["Rod B"]
until = { condition = "exchange_gas_cutoff", node = "Plate", cutoff = 100.0 }

[[procedure.phases]]
name = "both rods"
enable = ["Rod B"]

[simulation]
solver = "Radau"
max_time = 86400.0
initial_temperature = [300.0, 4.0]
fixed_boundary_nodes = ["Bath"]
```

The Plate cools through Rod A toward the Bath, held at 4 K. When the Plate reaches 100 K, Rod B activates and the Plate cools faster. The run has no stop condition, so it ends at `max_time` with the time-limit status.

The conduction parameter names (`material`, `area`, `length`) are placeholders until task 8.1 defines the conduction Link. `exchange_gas_cutoff` serves as a generic temperature threshold because it is the only condition registered in the first slice.
