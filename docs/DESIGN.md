# Design Specification: Pulse Tube Cooldown Simulation

## Scope
Simulate temperature and heat fluxes in a Maybell dilution refrigerator during cooldown, roughly from 300 K to when the MXP turns on (40-50 K). Version 1.0.0 targets a full cooldown simulation of a 6-plate system (SN5), with error comparable to or lower than the legacy model. For this task, we'll define "error" more precisely as the least-squares residual to log-temperature across all thermometers mounted on the system (except for MXP). 

## Public API
For now, we'll call the library Pulse Tube Cooldown, and we can have a canonical import statement:
```python
import pulse_tube_cooldown as ptc
```

#### Configuration files
Each system gets one, self-contained configuration file. Each model run takes exactly one configuration file as input.
Configuration files are written in TOML, with the following sections:
```toml
[meta]

[fridge]
[[fridge.nodes]]
name = "PRP"
composition = [["OFHC Cu RRR50", 11.8]]

[[fridge.nodes]]
name = "RGP"
composition = [["OFHC Cu RRR100", 6.35], ["304 SS", 5.0], ["Lead", 3.0]]

[[fridge.links]]
name = "PRP-RGP Switch 1"
type = "gas_gap_switch"
nodes = ["PRP", "RGP"] # order matters

[simulation]
end_node = "MXP"
end_temperature = 50.0 # K
max_time = 0.0 # s

```

#### Public functions
The library exposes the following public functions:
```python
def configure_model(path_to_config) --> Model: ...
    
def simulate_cooldown(model) --> Result: ...

def save_output(result, path_to_output) --> None: ...
```

So the canonical run is:
```python
model = ptc.configure_model(path_to_config)
result = ptc.simulate_cooldown(model)
ptc.save_output(result, path_to_output)
```

## Data Model
The important structures are listed here. The Model holds information about the simulation. The Graph describes the connectivity between components of the system, Nodes and Links, which may be made of Materials that are responsible for their own thermal properties. Simulation output is stored in a Result.

```python
@dataclass(frozen = True, slots = True, kw_only = True)
class Material:
    name: str
    specific_heat: Callable[[float], float]
    thermal_conductivity: Callable[[float], float]

@dataclass(frozen = True, slots = True)
class Node:
    name: str
    composition: tuple[tuple[Material, float], ...]

class Link(Protocol):
    @property
    def name(self) -> str: ...
    def calc_heat_flows(self, temperatures: NDArray[np.float64]) -> NDArray[np.float64]: ...

@dataclass(frozen = True, slots = True)
class Graph:
    nodes: tuple[Node, ...]
    links: tuple[Link, ...]
    connectivity: tuple[tuple[int, ...], ...]

@dataclass(frozen = True, slots = True)
class Model:
    graph: Graph
    initial_temperature: NDArray[np.float64]
    fixed_boundary_nodes: tuple[int, ...]
    end_node: int
    end_temperature: float
    max_time: float
    solver_options: Mapping[str, Any]

@dataclass(frozen = True, slots = True)
class Result:
    status: str
    model: Model
    time: NDArray[np.float64]
    temperature: NDArray[np.float64]
    heat_flow: NDArray[np.float64]
```

## Engine
The engine is responsible for constructing the thermal properties, heat flows, and linear system. 

```python
def calc_heat_capacity(nodes: Sequence[Node], temperature: NDArray[np.float64]) -> NDArray[np.float64]: ...

def calc_link_heat_flows(graph: Graph, temperature: NDArray[np.float64]) -> NDArray[np.float64]: ...

def calc_net_heat_flow(graph: Graph, temperature: NDArray[np.float64]) -> NDArray[np.float64]: ...

def calc_temperature_rate_of_change(model: Model, time: float, temperature: NDArray[np.float64]) -> NDArray[np.float64]: ...

def find_dist_until_stop(model: Model, time: float, temperature: NDArray[np.float64]) -> float: ...

def build_jacobian_sparsity(graph: Graph) -> scipy.sparse.csr_array: ...

def integrate_temperature(model: Model) -> scipy.integrate.OdeResult: ...

def reconstruct_heat_flow(graph: Graph, temperature: NDArray[np.float64]) -> NDArray[np.float64]: ...

def build_result(model: Model, solution: scipy.integrate.OdeResult, heat_flow: NDArray[np.float64]) -> Result: ...
```

## Dilution Refrigerator Library
Materials are stored in TOML files, in data/materials. Each one should take the format:
```toml
name = "OFHC Cu RRR100"
source = "NIST Cryogenic Material Properties"

[specific_heat]
form = "log_poly"
coefficients = [...]
range = [4.0, 300.0]

[thermal_conductivity]
form = "rational_poly"
coefficients = [...]
range = [4.0, 300.0]
```

This library has to handle some mapping from the config file to actual objects:
```python
LINK_TYPES = {
    "cooler": Cooler,
    "gas_gap_switch": GasGapSwitch,
    "exchange_gas": ExchangeGas,
    "radiation": Radiation,
    "multi_layer_insulation": MultiLayerInsulation,
    "conduction": Conduction,
    "joint_contact": JointContact
    ...
}

MATERIAL_FITS = {
    "log_poly": eval_log_poly,
    "rational_poly": eval_rational_poly,
    "table": eval_fit_from_table,
    ...
}
```



Known things we need (not organized yet):
- different types of links
- different types of materials
- functions to get specific heat and thermal conductivity for different materials

## Validation and Testing
