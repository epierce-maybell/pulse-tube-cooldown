# Design Specification: Pulse Tube Cooldown Simulation

## Scope
Simulate temperature and heat fluxes in a Maybell dilution refrigerator during cooldown, from room temperature to when condensing begins. Version 1.0.0 targets a full cooldown simulation of two 6-plate systems (SN09 for the standard model and FP1 for the fridge-plus model), with error comparable to or lower than the legacy model. For this task, we'll define "error" more precisely as the least-squares residual to log-temperature across all thermometers mounted on the system (except for MXP).

## Public API
For now, we'll call the library Pulse Tube Cooldown, and we can have a canonical import statement:
```python
import pulse_tube_cooldown as ptc
```

#### Configuration files
Each system takes one configuration file as an input that references the data subdirectory where necessary.

Configuration files are written in TOML, with the following sections:
```toml
[meta]
experiment = "Validation on SN5"
version = # insert commit hash or version tag

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

[[fridge.links]]
name = "Exchange gas"
type = "exchange_gas"
nodes = ["PRP", "RGP", "CFP", "STP", "ICP", "MXP"]

[procedure]
until = { condition = "condensing_start", node = "MXP", threshold = 10.0 }

[[procedure.phases]]
name = "precool"
until = { condition = "exchange_gas_cutoff", node = "CFP", cutoff = 70.0 }

[[procedure.phases]]
name = "cooldown"
disable = ["Exchange gas"]

[simulation]
solver = "Radau"
max_time = 604800.0 # still returns a Result
save_every = 100.0
initial_temperature = 300.0 # scalar or array of size [nodes]
fixed_boundary_nodes = ["Vacuum can"]
```

#### Public functions
The library exposes the following public functions:
```python
def configure_model(path_to_config) -> Model: ...


def simulate_cooldown(model) -> Result: ...


def save_output(result, path_to_output) -> None: ...
```

So the canonical run is:
```python
model = ptc.configure_model(path_to_config)
result = ptc.simulate_cooldown(model)
ptc.save_output(result, path_to_output)
```

#### Command-Line Interface
V1 includes a command-line interface (CLI) to run simulations for a given configuration file. It should be easy to invoke, with only a few options:
```bash
run-ptc --config "path_to_config" --output "path_to_output"
```

## Data model
The important structures are listed here. The Model holds information about the simulation, including the order of operations (see Phases and Procedure, below), and is responsible for constructing graphs. A Graph describes the connectivity between components of the system, Nodes and Links, which may be made of Materials that are responsible for their own thermal properties. Simulation output is stored in a Result.

```python
@dataclass(frozen=True, slots=True, kw_only=True)
class Material:
    name: str
    specific_heat: Callable[[NDArray[np.float64]], NDArray[np.float64]]
    thermal_conductivity: Callable[[NDArray[np.float64]], NDArray[np.float64]]
    conductivity_integral: Callable[[NDArray[np.float64]], NDArray[np.float64]]


@dataclass(frozen=True, slots=True, kw_only=True)
class Gas:
    name: str
    molar_mass: float
    heat_capacity_ratio: float
    accommodation_coeff: float
    thermal_conductivity: Callable[[NDArray[np.float64]], NDArray[np.float64]]
    saturation_pressure: Callable[[NDArray[np.float64]], NDArray[np.float64]]
    adsorption_energy: float
    adsorption_constant: float
    adsorption_capacity: float # per kg of charcoal


@dataclass(frozen=True, slots=True)
class Node:
    name: str
    composition: tuple[tuple[Material, float], ...]


class Link(Protocol):
    @property
    def name(self) -> str: ...

    @property
    def nodes(self) -> tuple[str, ...]: ...

    def calc_heat_flow(self, temperature: NDArray[np.float64]) -> NDArray[np.float64]: ...


@dataclass(frozen=True, slots=True)
class Graph:
    nodes: tuple[Node, ...]
    links: tuple[Link, ...]
    connectivity: tuple[tuple[int, ...], ...]

    def get_node_index(self, name: str) -> int: ...
    def get_link_index(self, name: str) -> int: ...


@dataclass(frozen=True, slots=True)
class Model:
    nodes: tuple[Node, ...]
    procedure: Procedure
    initial_temperature: NDArray[np.float64]
    fixed_boundary_nodes: tuple[int, ...]
    solver_options: Mapping[str, Any]
    graphs: tuple[Graph, ...]


@dataclass(frozen=True, slots=True, kw_only=True)
class PhaseSolution:
    phase: Phase
    end: Literal["until", "stop_condition", "max_time", "solver_failure"]
    time: NDArray[np.float64]
    temperature: NDArray[np.float64]
    message: str


@dataclass(frozen=True, slots=True)
class Result:
    status: str
    model: Model
    data: xr.Dataset
```

A Result's data is an xarray Dataset, saved as netCDF by `save_output`. The dimensions of `data` are `time`, `node`, and `link`, using node and link names as coordinates. The variables are `phase(time)` (Phase names), `temperature(time, node)` (K), `heat_flow(time, link, node)` (W into the node, NaN where link does not touch node or is inactive), and `connected(link, node)` (bool). At a phase switch, time holds a duplicate entry, one for each phase. Attributes record status, units, and run provenance.

## Phases and Procedure
During different phases of the cooldown, the simulation needs to activate different physics for certain components. A Procedure holds a sequence of Phases, each of which specifies the Links that are operating during that Phase (note that these are instantiated Links, not Link types, so a Phase can change detailed information during the run). It has a defined stopping point, which should be a function that takes (time, temperature) and returns a float that equals zero when the phase should end. The Procedure has a stop_condition that operates the same way.

The design follows from this:
```python
@dataclass(frozen=True, slots=True)
class Phase:
    name: str
    links: tuple[Link, ...]
    until: Callable[[float, NDArray[np.float64]], float] | None


@dataclass(frozen=True, slots=True)
class Procedure:
    phases: tuple[Phase, ...]
    stop_condition: Callable[[float, NDArray[np.float64]], float] | None
```

The [fridge] part of the configuration defines the full graph, with every link active. Each phase specifies which links should activate or deactivate in each phase, cumulatively, including an initial phase (which should only deactivate links) and an end phase (which should not have an `until` condition). Link names are unique. Redundant enables or disables should raise an error rather than silently succeeding.

## Engine
The engine is responsible for constructing the thermal properties, heat flows, and linear system.

```python
def calc_heat_capacity(nodes: Sequence[Node], temperature: NDArray[np.float64]) -> NDArray[np.float64]: ...


def calc_link_heat_flows(graph: Graph, temperature: NDArray[np.float64]) -> NDArray[np.float64]: ...


def calc_net_heat_flow(graph: Graph, temperature: NDArray[np.float64]) -> NDArray[np.float64]: ...


def calc_temperature_rate_of_change(
    graph: Graph, fixed_boundary_nodes: tuple[int, ...], time: float, temperature: NDArray[np.float64]
) -> NDArray[np.float64]: ...


def integrate_phase(
    model: Model, phase: Phase, time: float, temperature: NDArray[np.float64]
) -> PhaseSolution: ...


def run_procedure(model: Model) -> Sequence[PhaseSolution]: ...


def reconstruct_heat_flow(graph: Graph, temperature: NDArray[np.float64]) -> NDArray[np.float64]: ...


def build_result(model: Model, solutions: Sequence[PhaseSolution]) -> Result: ...
```

## Dilution refrigerator library
As a general rule, physical values should always be specified in SI units.

Materials are specified in individual TOML files in materials/data. Each entry should take the format:
```toml
name = "OFHC Cu RRR100"

[specific_heat]
form = "log_poly"
coefficients = []
range = [4.0, 300.0] # out-of-range will clamp to ends and raise a warning
source = "NIST Cryogenic Material Properties"


[thermal_conductivity]
form = "rational_poly"
coefficients = []
range = [4.0, 300.0]
source = "NIST Cryogenic Material Properties"
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

PROPERTY_FITS = {
    "log_poly": make_log_poly,
    "rational_poly": make_rational_poly,
    "table": make_fit_from_table,
    ...
}

CONDITIONS = {
    "exchange_gas_cutoff": exchange_gas_cutoff,
    "condensing_start": condensing_start,
    ...
}
```

Specific Links are frozen dataclasses that match the Link Protocol, where their fields are set by the TOML configuration file. There is a generic loader:
```python
def load_link(config: Mapping) -> Link: ...
```
which relies on a shared materials loader used by both nodes and links:
```python
def load_material(name: str) -> Material: ...
```
Links share a common schema, even though there is considerable variation in physics between specific types. Shared quantities should have shared names between types. This should look something like:
```python
@dataclass(frozen=True, slots=True)
class SpecificLinkTemplate:
    name: str
    nodes: tuple[str, ...]
    example_arg1: ...
    example_arg2: ...

    def calc_heat_flow(self, temperature_at_nodes: NDArray[np.float64]) -> NDArray[np.float64]:
        """Temperature only on this link's nodes. Positive values indicate heat flow into node."""
```

Conditions should generate functions that are positive before the condition is met (during cooldown), negative after, and zero where the switch happens. To stay agnostic to the specific Graph, the condition generator should map a node name rather than an index. For example:
```python
def exchange_gas_cutoff(graph: Graph, node: str, cutoff: float) -> Callable[[float, NDArray[np.float64]], float]:
    """Exchange gas is pumped out once a certain node hits a certain cutoff."""
    idx = graph.get_node_index(node)
    return lambda time, temperature: temperature[idx] - cutoff
```
Conditions and their arguments should be named for the physical event they mark, and for that reason, they are allowed to share implementations between two functions with different names.

## Specific physical processes
Some decisions are not obvious, but have leverage over the rest of the design.

#### Coupled capacity maps
We do not have extended load map data during cooldown for our pulse tubes. What we do have is one run per unit, post-tuning, on a test stand for exactly one cooldown curve. This gives us T1 and T2 with respect to time over one path. We also have a load map for each unit, but only for T2 $\in (2, 5)$ K and T1 $\in (40, 68)$ K.

This is an ongoing question, but for now, we can set up a first-order uncoupled model, where each $Q(T)$ path is separable. For each stage $i$ and heat flow $Q_i(T_i)$, we have

$$
Q_i(t) = -C_i(T_i(t)) \frac{dT_i}{dt} + Qp_i(t)
$$

for heat flow $Q$, thermal capacity $C$, temperature $T$, and parasitic heat load $Qp$. Conduction to the regenerator can be modeled simply as

$$
Q_R = G_R (T_{RGP} - T_R)
$$

where $Q_R$ is the heat flow, $G_R(T_R)$ is the regenerator conductance, $T_{RGP}$ is the temperature of the RGP node, and $T_R = T1 + x_R (T2 - T1)$ is the temperature of the regenerator tube at the normalized clamp position $x_R$. This assumes that $Q_R$ is an independent sink, so it does not influence $Q1$ or $Q2$ (which is almost certainly only an approximation of reality).

#### Gas pushing
Gas pushing is its own specific link type, although future scenarios may extend this to consider different methods. Here, we focus on `gas_pushing_pm5`, which is currently deployed. Rather than fit to data, we derive a simple model from first principles, where gas pushing is effectively a conductance between two nodes. Pressure varies on short time scales relative to heat flow, so we average over the pumping cycle. Then, the effective conductance $G$ between two nodes $k$ and $k + 1$ is given by

$$
G = \frac{5}{2} \frac{\Delta p}{P} \sum_{j > k} \frac{V_j}{T_j}
$$

where $\Delta p$ is the pressure swing for each cycle, $P$ is the period of the cycle, and $V_j$ and $T_j$ are the gas volume and temperature at each node $j$ that the gas fills. Heat flow is then simply: $Q = \pm G(\Delta T)$ where $\Delta T$ is the difference in temperature between each node.

#### Gas-gap heat switch
Gas-gap heat switches are a specific link type, always connecting only two nodes, handled here with a reduced-complexity model. Heat flow $Q$ through the switch is given by:

$$
Q = \frac{(T_h - T_c)}{\dfrac{1}{hA} + \sum_b \dfrac{L_b}{k_b A_b}},
$$
where $T_h$ and $T_c$ are the temperatures of the warmer and colder adjacent nodes, respectively, $A$ is the gap area, and each bracket $b$ has a length $L_b$, area $A_b$, and conductivity $k_b$. The gap's conductance per unit area $h$ is a combination of the high-pressure and low-pressure behavior:

$$
\frac{1}{h} = \frac{d}{k_g} + \frac{1}{h_m}, \qquad h_m = \alpha \frac{\gamma + 1}{\gamma - 1} \sqrt{\frac{R}{8 \pi M \bar{T}}} P,
$$

for gap width $d$, gas thermal conductivity $k_g$, accommodation coefficient $\alpha$, heat capacity ratio $\gamma$, molar mass $M$, and ideal gas constant $R$.

The pressure $P$ follows from conservation. Each switch seals after filling to pressure $P_{fill}$ at temperature $T_{fill}$, so the total amount of gas is $n = P_{fill} V / (R T_{fill})$. The gas is either free or held by the charcoal, so we have:

$$
n = \frac{P V}{R \bar{T}} + n_s \frac{K P}{1 + K P}, \qquad K = K_0 \, e^{E / R \bar{T}},
$$

where $n_s$ is the charcoal mass times its capacity and $K_0$ and $E$ are the gas adsorption constants on the charcoal. $P$ is the positive root of this quadratic, capped at some $P_{sat}(T_c)$. Gas properties and $k_b$ are evaluated at the mean temperature $\bar{T} = (T_h + T_c) / 2$.

The configuration for the switch should look something like:
```toml
[[fridge.links]]
name = "CFP-STP Switch"
type = "gas_gap_switch"
nodes = ["CFP", "STP"]
gas = "He-4"
fill_pressure = 1.0e5      # Pa
fill_temperature = 295.0   # K
gas_volume = 2.5e-6        # m^3
gap = 5.1e-4               # m
area = 4.6e-3              # m^2
charcoal_mass = 1.8e-4     # kg
brackets = [{ material = "OFHC Cu RRR100", area = 1.0e-4, length = 0.02 }]
```

## Validation and testing
Validation experiments should live in simulations/validation/identifier, where identifier is usually a subdirectory that stores the runner script, the configuration file(s), and output.

In the future, we will want regression tests. The first validation experiment(s) that pass, leading us to v1 release, will become the regression tests for later iterations of the model.

Our test suite should include unit tests for each discrete unit of code that is worth testing (e.g., could reasonably not do what we expect it to). Analytic tests are appropriate where we have access to analytic solutions

## Package layout
The repository structure should have the following shape, although this is not a comprehensive description (specific files may be added or removed later, during implementation):
```bash
src/pulse_tube_cooldown/
    __init__.py
    __main__.py
    config.py
    output.py
    core/
        __init__.py
        types.py
        errors.py
    engine/
        __init__.py
        equations.py
        solver.py
    library/
        __init__.py
        conditions.py
        property_fits.py
        materials/
            __init__.py
            load.py
            data/
                ofhc_cu_rrr100.toml
                ss_304.toml
                ...
        gases/
            __init__.py
            load.py
            data/...
        links/
            __init__.py
            conduction.py
            radiation.py
            cooler.py
            gas_gap_switch.py
            gas_pushing.py
            exchange_gas.py
            multi_layer_insulation.py
            joint_contact.py
    plot/
        __init__.py
        graph.py
        result.py

tests/
    core/...
    engine/...
    library/...
    regression/...
```
