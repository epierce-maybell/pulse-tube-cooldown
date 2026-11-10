"""Tests for loading links from configuration."""

from functools import partial

import numpy as np
import pytest

from pulse_tube_cooldown.core.errors import ConfigurationError
from pulse_tube_cooldown.library.links import radiation
from pulse_tube_cooldown.library.links.load import load_link
from pulse_tube_cooldown.library.surfaces.load import load_surface

_ROD = {
    "name": "Rod",
    "type": "conduction",
    "nodes": ["CFP", "STP"],
    "material": "304 SS",
    "area": 1.0e-4,
    "length": 0.1,
}


def test_loads_conduction_with_its_material():
    link = load_link(_ROD)

    assert (link.name, link.nodes) == ("Rod", ("CFP", "STP"))
    # Integrating the NIST 304 SS conductivity equation from 10 to 50 K gives 136.400814 W/m.
    flow = link.calc_energy_balance_terms(np.array([10.0, 50.0])).heat_flow
    np.testing.assert_allclose(flow, [0.136400814, -0.136400814], rtol=1e-7)


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


def test_loads_radiation_with_its_surfaces(surface_directory):
    link = load_link(_SHIELD)

    assert (link.name, link.nodes, link.area) == ("4K shield to still shield", ("STP", "CFP"), (0.119, 0.45))
    # The two gray surfaces have emissivities 0.02 and 0.1; their radiosity network carries 2.564784 mW.
    flow = link.calc_energy_balance_terms(np.array([45.0, 70.0])).heat_flow
    np.testing.assert_allclose(flow, [0.00256478421, -0.00256478421], rtol=1e-7)


@pytest.mark.parametrize("missing", ["nodes", "area", "surface"])
def test_missing_radiation_parameter_is_named(surface_directory, missing):
    config = {key: value for key, value in _SHIELD.items() if key != missing}

    with pytest.raises(ConfigurationError, match=missing):
        load_link(config)


def test_unknown_surface_names_the_surface(surface_directory):
    with pytest.raises(ConfigurationError, match="Black Paint"):
        load_link(_SHIELD | {"surface": ["Bright Finish", "Black Paint"]})


def test_radiation_rejects_inner_area_larger_than_outer_area(surface_directory):
    with pytest.raises(ConfigurationError, match="4K shield to still shield"):
        load_link(_SHIELD | {"area": [0.5, 0.1]})


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


def test_loads_multi_layer_insulation():
    link = load_link(_BLANKET)

    assert (link.name, link.nodes) == ("Stage 1 blanket", ("Stage 1", "Chamber"))
    # At 50 and 300 K, spacer and radiation terms contribute 0.14842610 and 0.07476105 W.
    flow = link.calc_energy_balance_terms(np.array([50.0, 300.0])).heat_flow
    np.testing.assert_allclose(flow, [0.223187155, -0.223187155], rtol=1e-7)


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


def test_loads_joint_contact():
    link = load_link(_JOINT)

    assert (link.name, link.nodes) == ("Stage 2 flange", ("PT stage 2", "CFP"))
    # Mean temperature 8.4 K: 19 W/K * 2^1.1 * 4 K = 162.909566 W.
    flow = link.calc_energy_balance_terms(np.array([6.4, 10.4])).heat_flow
    np.testing.assert_allclose(flow, [162.909566, -162.909566], rtol=1e-7)


@pytest.mark.parametrize(
    "missing",
    ["nodes", "force", "conductance_per_force", "reference_temperature", "temperature_exponent", "source"],
)
def test_missing_joint_contact_parameter_is_named(missing):
    config = {key: value for key, value in _JOINT.items() if key != missing}

    with pytest.raises(ConfigurationError, match=missing):
        load_link(config)
