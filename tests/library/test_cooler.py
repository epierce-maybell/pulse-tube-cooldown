"""Requirements for a manifold link."""

from dataclasses import replace

import numpy as np
import pytest

from pulse_tube_cooldown.core.types import Graph, Material, Node
from pulse_tube_cooldown.engine.equations import calc_net_heat_flow
from pulse_tube_cooldown.library.gases import helium
from pulse_tube_cooldown.library.links.cooler import (
    Compressor,
    PulseTubeCooler,
    PulseTubeHead,
    RegeneratorLayer,
    calc_charge,
)


@pytest.fixture
def ideal_head(monkeypatch):
    # Synthetic ideal gas: rho = P/T, c_p = 10, beta = 1/T, and a matrix with c_m = T.
    monkeypatch.setattr(helium, "density", lambda temperature, pressure: pressure / np.asarray(temperature))
    monkeypatch.setattr(helium, "isobaric_specific_heat", lambda temperature, pressure: np.full_like(temperature, 10.0))
    monkeypatch.setattr(helium, "isobaric_expansion_coefficient", lambda temperature, pressure: 1 / temperature)
    material = Material(
        name="Matrix",
        density=1.0,
        specific_heat=lambda temperature: temperature,
        thermal_conductivity=None,
        conductivity_integral=None,
    )
    return PulseTubeHead(
        name="A",
        stages=("S1", "S2"),
        columns=(
            (RegeneratorLayer(node=Node("A matrix 1", ((material, 2.0),)), pore_volume=1e-4, length_fraction=1.0),),
            (RegeneratorLayer(node=Node("A matrix 2", ((material, 3.0),)), pore_volume=2e-4, length_fraction=1.0),),
        ),
        pulse_tube_volumes=(6e-5, 2e-4),
        clamp_fraction=0.5,
        packing_constants=(0.8, 1.3),
        conductances=(1e-10, 2e-12),
        exchanger_effectiveness=0.3,
        frequency=1.2,
    )


@pytest.mark.parametrize("head_count, displacement, relief", [(1, 1e-3, 1.9e6), (2, 1e-3, 1.9e6), (2, 0.1, 1.1e6)])
def test_manifold_flows_match_analytic_ideal_gas(ideal_head, head_count, displacement, relief):
    heads = (ideal_head,)
    if head_count == 2:
        layers = tuple(
            tuple(replace(layer, node=replace(layer.node, name=layer.node.name.replace("A", "B"))) for layer in column)
            for column in ideal_head.columns
        )
        heads += (replace(ideal_head, name="B", columns=layers, frequency=0.9, conductances=(2e-10, 3e-12)),)
    cooler = PulseTubeCooler(
        name="Cooler",
        warm_flange="Warm",
        heads=heads,
        compressors=(
            Compressor(displacement_rate=displacement, relief_swing=1.9e6),
            Compressor(displacement_rate=displacement / 2, relief_swing=relief),
        ),
        charge=180.0,
        warm_volume=0.03,
        waveform_factor=0.8,
    )
    states = [np.array([300.0, 210.0, 120.0, 70.0, 20.0]), np.array([300.0, 180.0, 120.0, 50.0, 20.0])][:head_count]
    by_name = {}
    for head, state in zip(heads, states, strict=True):
        by_name.update(
            zip(("Warm", head.matrix_nodes[0].name, "S1", head.matrix_nodes[1].name, "S2"), state, strict=True)
        )
    temperature = np.array([by_name[name] for name in cooler.nodes])
    # Gas inventory is P times sum(V/T); pulse-tube averages integrate 1/T over their linear profiles.
    inventory = [
        1e-4 / t[1] + 2e-4 / t[3] + 6e-5 * np.log(300 / 120) / 180 + 2e-4 * np.log(300 / 20) / 280 for t in states
    ]
    mean = 180.0 / (0.03 / 300 + sum(inventory))
    swing = min(
        relief,
        1.9e6,
        mean
        * (1.5 * displacement / 300)
        / (sum(h.frequency * v for h, v in zip(heads, inventory, strict=True)) + 1.5 * displacement / 600),
    )
    expected = dict.fromkeys(cooler.nodes, 0.0)
    for head, t, volume in zip(heads, states, inventory, strict=True):
        flows = (
            head.frequency
            * swing
            * np.array([volume, volume - 1e-4 / t[1], 2e-4 / t[3] + 2e-4 * np.log(15) / 280, 2e-4 * np.log(15) / 280])
        )
        capacities = np.array([2 * t[1], 2 * t[1], 3 * t[3], 3 * t[3]])
        loss = 2 * np.array([0.8, 0.8, 1.3, 1.3]) * (flows * 10) ** 2 * (-np.diff(t)) / (head.frequency * capacities)
        h1, h2 = np.array(head.conductances) * (0.8 * swing / 2) ** 2 / 2
        # The stage-2 column ends at its node, which sinks h2; the exchanger carries eta m c_p (T_m - T_2) on to
        # stage 2. With beta T = 1 the stage-1 boundary carries only the finite-capacity flow eps / (1 + eps) E_2, for
        # the stage-2 node's eps = V rho c_p / C and E_2 = h2 (T_1 / T_m)^(1 / (1 + eps)).
        exchange = 0.3 * flows[3] * 10 * (t[3] - t[4])
        ratio = 2e-4 * (mean / t[3]) * 10 / (3 * t[3])
        pressure_driven = ratio / (1 + ratio) * h2 * (t[2] / t[3]) ** (1 / (1 + ratio))
        values = [
            h1 + h2 - loss[0],
            loss[0] - loss[1],
            loss[1] - loss[2] - h1 - pressure_driven,
            loss[2] + pressure_driven - exchange - h2,
            exchange,
        ]
        for name, value in zip(
            ("Warm", head.matrix_nodes[0].name, "S1", head.matrix_nodes[1].name, "S2"), values, strict=True
        ):
            expected[name] += value
    result = cooler.calc_heat_flow(temperature)
    assert len(cooler.nodes) == 3 + 2 * head_count
    np.testing.assert_allclose(result, [expected[name] for name in cooler.nodes], rtol=1e-10, atol=1e-10)
    assert result.sum() == pytest.approx(0.0, abs=1e-9)


def test_cooler_flows_follow_graph_node_order(ideal_head):
    cooler = PulseTubeCooler(
        name="Cooler",
        warm_flange="Warm",
        heads=(ideal_head,),
        compressors=(Compressor(displacement_rate=0.1, relief_swing=1e6),),
        charge=calc_charge((ideal_head,), 1.7e6, 300.0, 0.03),
        warm_volume=0.03,
    )
    by_name = {node.name: node for node in ideal_head.matrix_nodes}
    composition = ideal_head.matrix_nodes[0].composition
    by_name.update((name, Node(name, composition)) for name in ("Warm", "S1", "S2", "Separate"))
    graph = Graph(
        nodes=tuple(by_name[name] for name in ("S2", "Separate", "A matrix 2", "Warm", "S1", "A matrix 1")),
        links=(cooler,),
    )

    # At uniform temperature the ideal gas carries no column loss. At the 1 MPa relief swing,
    # H = C * (dP / 2)^2 / 2 gives 12.5 W and 0.25 W of cooling at stage 1 and the stage-2 column's cold end, rejected
    # to the warm flange. The stage-2 column's pore gas, with eps = V rho c_p / C at 1.7 MPa and 300 K, still moves
    # eps / (1 + eps) H_2 from stage 1 to the cold end, and the exchanger carries nothing on to stage 2.
    ratio = 2e-4 * (1.7e6 / 300) * 10 / (3 * 300)
    moved = ratio / (1 + ratio) * 0.25
    np.testing.assert_allclose(
        calc_net_heat_flow(graph, np.full(6, 300.0)),
        [0.0, 0.0, -0.25 + moved, 12.75, -12.5 - moved, 0.0],
        atol=1e-12,
    )
