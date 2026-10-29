"""Tests for loading links from configuration."""

from functools import partial

import numpy as np
import pytest

from pulse_tube_cooldown.core.errors import ConfigurationError
from pulse_tube_cooldown.library.links import radiation
from pulse_tube_cooldown.library.links.conduction import Conduction, load_conduction
from pulse_tube_cooldown.library.links.load import LINK_TYPES, load_link
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
