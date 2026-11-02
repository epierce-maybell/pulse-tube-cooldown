"""Tests for loading links from configuration."""

from functools import partial

import numpy as np
import pytest

from pulse_tube_cooldown.core.errors import ConfigurationError
from pulse_tube_cooldown.library.links import radiation
from pulse_tube_cooldown.library.links.conduction import Conduction, load_conduction
from pulse_tube_cooldown.library.links.joint_contact import JointContact, load_joint_contact
from pulse_tube_cooldown.library.links.load import LINK_TYPES, load_link
from pulse_tube_cooldown.library.links.multi_layer_insulation import (
    MultiLayerInsulation,
    load_multi_layer_insulation,
)
from pulse_tube_cooldown.library.links.radiation import Radiation, load_radiation
from pulse_tube_cooldown.library.materials.load import load_material
from pulse_tube_cooldown.library.surfaces.load import load_surface

_ROD = {
    "name": "Rod",
    "type": "conduction",
    "nodes": ["CFP", "STP"],
    "material": "304 SS",
    "area": 1.0e-4,
    "length": 0.1,
}


def test_registers_conduction():
    assert LINK_TYPES["conduction"] is load_conduction


def test_loads_conduction_with_its_material():
    link = load_link(_ROD)

    assert isinstance(link, Conduction)
    assert (link.name, link.nodes, link.area, link.length) == ("Rod", ("CFP", "STP"), 1.0e-4, 0.1)
    assert link.material.name == "304 SS"
    temperature = np.array([10.0, 50.0, 250.0])
    np.testing.assert_array_equal(
        link.material.thermal_conductivity(temperature), load_material("304 SS").thermal_conductivity(temperature)
    )


def test_unknown_type_names_the_link_and_known_types():
    with pytest.raises(ConfigurationError, match=r"(?s)Rod.*'strut'.*conduction"):
        load_link(_ROD | {"type": "strut"})


def test_missing_type_names_the_link():
    config = {key: value for key, value in _ROD.items() if key != "type"}

    with pytest.raises(ConfigurationError, match="Rod"):
        load_link(config)


@pytest.mark.parametrize("missing", ["name", "nodes", "material", "area", "length"])
def test_missing_parameter_is_named(missing):
    config = {key: value for key, value in _ROD.items() if key != missing}

    with pytest.raises(ConfigurationError, match=missing):
        load_link(config)


def test_unknown_material_names_the_material():
    with pytest.raises(ConfigurationError, match="Unobtainium"):
        load_link(_ROD | {"material": "Unobtainium"})


_SHIELD = {
    "name": "4K shield to still shield",
    "type": "radiation",
    "nodes": ["STP", "CFP"],
    "area": [0.119, 0.45],
    "surface": ["Bright Finish", "Dull Finish"],
}


@pytest.fixture
def surface_directory(tmp_path, monkeypatch):
    for name, emissivity in (("Bright Finish", 0.02), ("Dull Finish", 0.1)):
        (tmp_path / f"{name.split()[0].lower()}.toml").write_text(
            f'name = "{name}"\n\n[emissivity]\nform = "table"\n'
            f"table = [[4.0, {emissivity}], [300.0, {emissivity}]]\nrange = [4.0, 300.0]\n",
            encoding="utf-8",
        )
    monkeypatch.setattr(radiation, "load_surface", partial(load_surface, directory=tmp_path))


def test_registers_radiation():
    assert LINK_TYPES["radiation"] is load_radiation


def test_loads_radiation_with_its_surfaces(surface_directory):
    link = load_link(_SHIELD)

    assert isinstance(link, Radiation)
    assert (link.name, link.nodes, link.area) == ("4K shield to still shield", ("STP", "CFP"), (0.119, 0.45))
    assert [surface.name for surface in link.surface] == ["Bright Finish", "Dull Finish"]
    np.testing.assert_allclose([s.emissivity(np.array([77.0]))[0] for s in link.surface], [0.02, 0.1])


@pytest.mark.parametrize("missing", ["nodes", "area", "surface"])
def test_missing_radiation_parameter_is_named(surface_directory, missing):
    config = {key: value for key, value in _SHIELD.items() if key != missing}

    with pytest.raises(ConfigurationError, match=missing):
        load_link(config)


@pytest.mark.parametrize(
    "pair", [{"nodes": ["STP", "CFP", "PRP"]}, {"area": [0.119]}, {"surface": ["Bright Finish"] * 3}]
)
def test_radiation_pairs_must_hold_two_entries(surface_directory, pair):
    with pytest.raises(ValueError, match="unpack"):
        load_link(_SHIELD | pair)


def test_unknown_surface_names_the_surface(surface_directory):
    with pytest.raises(ConfigurationError, match="Black Paint"):
        load_link(_SHIELD | {"surface": ["Bright Finish", "Black Paint"]})


_BLANKET = {
    "name": "Stage 1 blanket",
    "type": "multi_layer_insulation",
    "nodes": ["Stage 1", "Chamber"],
    "area": 0.05,
    "layers": 4,
    "layer_density": 2500.0,
    "spacer_coefficient": 6.789e-13,
    "layer_density_exponent": 2.56,
    "radiation_coefficient": 5.39e-10,
    "emissivity": 0.03,
    "source": "Lockheed 1974, NTRS 20190030285",
}


def test_registers_multi_layer_insulation():
    assert LINK_TYPES["multi_layer_insulation"] is load_multi_layer_insulation


def test_loads_multi_layer_insulation():
    link = load_link(_BLANKET)

    assert isinstance(link, MultiLayerInsulation)
    assert (link.name, link.nodes) == ("Stage 1 blanket", ("Stage 1", "Chamber"))
    assert (
        link.area,
        link.layers,
        link.layer_density,
        link.spacer_coefficient,
        link.layer_density_exponent,
        link.radiation_coefficient,
        link.emissivity,
        link.source,
    ) == (0.05, 4, 2500.0, 6.789e-13, 2.56, 5.39e-10, 0.03, "Lockheed 1974, NTRS 20190030285")


@pytest.mark.parametrize(
    "missing",
    [
        "nodes",
        "area",
        "layers",
        "layer_density",
        "spacer_coefficient",
        "layer_density_exponent",
        "radiation_coefficient",
        "emissivity",
        "source",
    ],
)
def test_missing_multi_layer_insulation_parameter_is_named(missing):
    config = {key: value for key, value in _BLANKET.items() if key != missing}

    with pytest.raises(ConfigurationError, match=missing):
        load_link(config)


_JOINT = {
    "name": "Stage 2 flange",
    "type": "joint_contact",
    "nodes": ["PT stage 2", "CFP"],
    "force": 1.0e4,
    "conductance_per_force": 1.9e-3,
    "reference_temperature": 4.2,
    "temperature_exponent": 1.1,
    "source": "Dhuley 2019, Tables 2 and 5, Au-plated Cu",
}


def test_registers_joint_contact():
    assert LINK_TYPES["joint_contact"] is load_joint_contact


def test_loads_joint_contact():
    link = load_link(_JOINT)

    assert isinstance(link, JointContact)
    assert (link.name, link.nodes) == ("Stage 2 flange", ("PT stage 2", "CFP"))
    assert (
        link.force,
        link.conductance_per_force,
        link.reference_temperature,
        link.temperature_exponent,
        link.source,
    ) == (1.0e4, 1.9e-3, 4.2, 1.1, "Dhuley 2019, Tables 2 and 5, Au-plated Cu")


@pytest.mark.parametrize(
    "missing",
    ["nodes", "force", "conductance_per_force", "reference_temperature", "temperature_exponent", "source"],
)
def test_missing_joint_contact_parameter_is_named(missing):
    config = {key: value for key, value in _JOINT.items() if key != missing}

    with pytest.raises(ConfigurationError, match=missing):
        load_link(config)
