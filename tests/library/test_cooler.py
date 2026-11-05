"""Requirements for the cooler link: one manifold's heads and drive, configured from TOML."""

from itertools import accumulate

import numpy as np
import pytest

from pulse_tube_cooldown.library.gases import helium
from pulse_tube_cooldown.library.heads import load_head_model
from pulse_tube_cooldown.library.links.cooler import (
    Compressor,
    calc_charge,
    calc_gas_mass,
    calc_heat_flow,
    calc_mean_pressure,
)
from pulse_tube_cooldown.library.links.cooler.link import load_cooler, load_matrix_nodes
from pulse_tube_cooldown.library.links.load import LINK_TYPES, load_link_nodes

# Drive inputs from data/compressors.md; placeholder head parameters pending task 2.3's fit.
_ROOM_TEMPERATURE = 295.0
_WARM_VOLUME = 19.6e-3
_FILL_PRESSURE = 1.70e6
_DISPLACEMENT_RATE = 5.72e-3
_RELIEF_SWING = 1.93e6
_PACKING = (0.15, 0.05)
_RESISTANCE = (1e10, 3e10)
_FREQUENCY = 1.4

_MATRIX_NAMES = [
    "PT stage 1 304 SS",
    "PT stage 1 Pb",
    "PT stage 2 Pb",
    "PT stage 2 Er3Ni",
    "PT stage 2 HoCu2",
    *(f"PT stage 2 Gd2O2S {number}" for number in range(1, 9)),
]
# Pore volume in m^3 of each matrix node, from pt1_5.toml, with Gd2O2S split into eight nodes.
_PORE_VOLUME = [2.8459e-4, 6.0081e-5, 1.1956e-5, 2.3911e-5, 2.3911e-5, *[1.4945e-5 / 8] * 8]


def _write_drive(directory, *, compressors=1, warm_volume=_WARM_VOLUME, displacement_rate=_DISPLACEMENT_RATE):
    path = directory / f"drive_{compressors}_{warm_volume}_{displacement_rate}.toml"
    path.write_text(
        f"warm_volume = {warm_volume}\n"
        f"fill_pressure = {_FILL_PRESSURE}\n"
        f"fill_temperature = {_ROOM_TEMPERATURE}\n"
        + f"[[compressors]]\ndisplacement_rate = {displacement_rate}\nrelief_swing = {_RELIEF_SWING}\n" * compressors
        + '[sources]\nwarm_volume = "test"\n',
        encoding="utf-8",
    )
    return path


def _write_parameters(directory):
    path = directory / "pt.toml"
    path.write_text(
        f"packing = {list(_PACKING)}\n"
        f"resistance = {list(_RESISTANCE)}\n"
        f"frequency = {_FREQUENCY}\n"
        '[sources]\npacking = "placeholder"\n',
        encoding="utf-8",
    )
    return path


def _config(drive, parameters, *heads):
    """A cooler table whose heads each take a name and stage node names."""
    return {
        "name": "Cooler",
        "type": "cooler",
        "warm_flange": "Warm flange",
        "drive": str(drive),
        "heads": [
            {"name": name, "stages": list(stages), "model": "PT1.5", "parameters": str(parameters)}
            for name, stages in heads
        ],
    }


def _cooler(tmp_path, *heads, **drive):
    heads = heads or (("PT", ("Stage 1", "Stage 2")),)
    return load_cooler(_config(_write_drive(tmp_path, **drive), _write_parameters(tmp_path), *heads))


def _chain_temperature(stage_temperature, warm_temperature=_ROOM_TEMPERATURE):
    """Chain temperatures with each section's nodes linear in position between its ends' temperatures."""
    model = load_head_model("PT1.5")
    ends = (warm_temperature, *stage_temperature)
    sections = []
    for section, warm, cold in zip(model.regenerators, ends[:-1], ends[1:], strict=True):
        fractions = [
            layer.length_fraction / layer.node_count for layer in section.layers for _ in range(layer.node_count)
        ]
        midpoints = np.array(
            [end - fraction / 2 for end, fraction in zip(accumulate(fractions), fractions, strict=True)]
        )
        sections.append(warm + (cold - warm) * midpoints)
    return np.concatenate(([warm_temperature], sections[0], [ends[1]], sections[1], [ends[2]]))


def test_cooler_is_a_registered_link_type():
    assert LINK_TYPES["cooler"] is load_cooler


def test_a_cooler_table_defines_its_matrix_nodes(tmp_path):
    config = _config(tmp_path, tmp_path, ("PT", ("Stage 1", "Stage 2")))

    assert [node.name for node in load_link_nodes(config)] == _MATRIX_NAMES
    assert load_link_nodes({"type": "conduction"}) == ()


def test_nodes_run_along_each_head_chain_from_the_warm_flange(tmp_path):
    cooler = _cooler(tmp_path)

    assert cooler.nodes == ("Warm flange", *_MATRIX_NAMES[:2], "Stage 1", *_MATRIX_NAMES[2:], "Stage 2")


def test_matrix_nodes_hold_only_their_share_of_the_matrix(tmp_path):
    nodes = load_matrix_nodes(_config(tmp_path, tmp_path, ("PT", ("Stage 1", "Stage 2"))))
    masses = [0.9758, 1.112, 0.2212, 0.353, 0.3597, *[0.1791 / 8] * 8]
    materials = ["304 SS", "Pb", "Pb", "Er3Ni", "HoCu2", *["Gd2O2S"] * 8]

    assert [node.name for node in nodes] == _MATRIX_NAMES
    assert [len(node.composition) for node in nodes] == [1] * 13
    assert [node.composition[0][0].name for node in nodes] == materials
    np.testing.assert_allclose([node.composition[0][1] for node in nodes], masses)


def test_heads_on_one_manifold_each_have_their_own_matrix_nodes(tmp_path):
    heads = (("PT A", ("Stage 1", "Stage 2")), ("PT B", ("Stage 1", "Stage 2")))
    cooler = _cooler(tmp_path, *heads)
    names = [node.name for node in load_matrix_nodes(_config(tmp_path, tmp_path, *heads))]

    assert names == [name.replace("PT", "PT A") for name in _MATRIX_NAMES] + [
        name.replace("PT", "PT B") for name in _MATRIX_NAMES
    ]
    assert cooler.nodes == (
        "Warm flange",
        *names[:2],
        "Stage 1",
        *names[2:13],
        "Stage 2",
        *names[13:15],
        *names[15:],
    )


def test_reads_head_parameters_and_drive_from_their_files(tmp_path):
    cooler = _cooler(tmp_path, compressors=2)
    (mounted,) = cooler.heads

    head = mounted.head
    assert (head.model.name, head.packing, head.resistance, head.frequency) == (
        "PT1.5",
        _PACKING,
        _RESISTANCE,
        _FREQUENCY,
    )
    assert mounted.stages == ("Stage 1", "Stage 2")
    assert cooler.compressors == (Compressor(displacement_rate=_DISPLACEMENT_RATE, relief_swing=_RELIEF_SWING),) * 2
    assert cooler.warm_volume == _WARM_VOLUME
    # HEPAK's density at 295 K and 1.70 MPa times the warm volume and the head's pores, tubes, reservoirs, and dead
    # volume, from pt1_5.toml.
    head_volume = sum(_PORE_VOLUME) + 2.5077e-4 + 2.0086e-4 + 2 * 5.2800e-3 + 4.5e-5
    assert cooler.charge == pytest.approx(2.7524272578 * (_WARM_VOLUME + head_volume), rel=1e-4)


def test_relative_file_paths_resolve_from_the_repository_root(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    cooler = load_cooler(
        _config(
            "tests/integration/fixtures/cooler/drive.toml",
            "tests/integration/fixtures/cooler/pt.toml",
            ("PT", ("Stage 1", "Stage 2")),
        )
    )

    assert cooler.heads[0].head.model.name == "PT1.5"
    assert len(cooler.compressors) == 1


def test_heat_flow_at_rest_with_capped_swing_is_the_heads(tmp_path):
    # At rest, the manifold holds its charge at the fill pressure. Ample displacement drives the swing past the relief
    # point, where the relief valve holds it at the relief swing.
    cooler = _cooler(tmp_path, displacement_rate=1.0)
    temperature = np.full(len(cooler.nodes), _ROOM_TEMPERATURE)

    expected = calc_heat_flow(cooler.heads[0].head, temperature, _FILL_PRESSURE, _RELIEF_SWING)

    np.testing.assert_allclose(cooler.calc_energy_balance_terms(temperature).heat_flow, expected, rtol=1e-6, atol=1e-9)


def test_warm_flange_takes_no_heat(tmp_path):
    cooler = _cooler(tmp_path)

    assert cooler.calc_energy_balance_terms(_chain_temperature((60.0, 10.0))).heat_flow[0] == 0.0


def test_identical_heads_on_a_scaled_drive_each_match_one_head(tmp_path):
    temperature = _chain_temperature((60.0, 10.0))
    single = _cooler(tmp_path).calc_energy_balance_terms(temperature).heat_flow
    pair = _cooler(
        tmp_path,
        ("PT A", ("Stage 1 A", "Stage 2 A")),
        ("PT B", ("Stage 1 B", "Stage 2 B")),
        compressors=2,
        warm_volume=2 * _WARM_VOLUME,
    )

    flow = pair.calc_energy_balance_terms(np.concatenate((temperature, temperature[1:]))).heat_flow

    np.testing.assert_allclose(flow[: len(single)], single, rtol=1e-6)
    np.testing.assert_allclose(flow[len(single) :], single[1:], rtol=1e-6)


def test_heads_sharing_stages_add_their_stage_heat_flows(tmp_path):
    temperature = _chain_temperature((60.0, 10.0))
    single = _cooler(tmp_path).calc_energy_balance_terms(temperature).heat_flow
    pair = _cooler(
        tmp_path,
        ("PT A", ("Stage 1", "Stage 2")),
        ("PT B", ("Stage 1", "Stage 2")),
        compressors=2,
        warm_volume=2 * _WARM_VOLUME,
    )
    stage_1, stage_2 = 3, 15
    b_nodes = [*range(1, 3), *range(4, 15)]
    pair_temperature = np.concatenate((temperature, temperature[b_nodes]))

    flow = pair.calc_energy_balance_terms(pair_temperature).heat_flow

    np.testing.assert_allclose(flow[[stage_1, stage_2]], 2 * single[[stage_1, stage_2]], rtol=1e-6)
    np.testing.assert_allclose(flow[16:], single[b_nodes], rtol=1e-6)


def test_heat_capacity_at_rest_is_pore_helium_at_fill_pressure(tmp_path):
    cooler = _cooler(tmp_path)
    temperature = np.full(len(cooler.nodes), _ROOM_TEMPERATURE)
    volumetric = float(
        helium.density(_ROOM_TEMPERATURE, _FILL_PRESSURE)
        * helium.isobaric_specific_heat(_ROOM_TEMPERATURE, _FILL_PRESSURE)
    )
    expected = np.array([0.0, *_PORE_VOLUME[:2], 0.0, *_PORE_VOLUME[2:], 0.0]) * volumetric

    np.testing.assert_allclose(cooler.calc_energy_balance_terms(temperature).heat_capacity, expected, rtol=1e-8)


def test_heat_capacity_takes_the_running_mean_pressure(tmp_path):
    cooler = _cooler(tmp_path)
    model = load_head_model("PT1.5")
    temperature = _chain_temperature((40.0, 4.0))
    charge = calc_charge(
        [model], warm_volume=_WARM_VOLUME, fill_pressure=_FILL_PRESSURE, fill_temperature=_ROOM_TEMPERATURE
    )
    mean_pressure = calc_mean_pressure(
        lambda pressure: calc_gas_mass(model, temperature, pressure),
        charge,
        warm_volume=_WARM_VOLUME,
        warm_temperature=_ROOM_TEMPERATURE,
    )
    pore = np.array([0.0, *_PORE_VOLUME[:2], 0.0, *_PORE_VOLUME[2:], 0.0])
    expected = (
        pore * helium.density(temperature, mean_pressure) * helium.isobaric_specific_heat(temperature, mean_pressure)
    )

    assert mean_pressure < 0.95 * _FILL_PRESSURE
    np.testing.assert_allclose(cooler.calc_energy_balance_terms(temperature).heat_capacity, expected, rtol=1e-6)
