"""Tests for loading links from configuration."""

import numpy as np
import pytest

from pulse_tube_cooldown.core.errors import ConfigurationError
from pulse_tube_cooldown.library.links.conduction import Conduction, load_conduction
from pulse_tube_cooldown.library.links.load import LINK_TYPES, load_link
from pulse_tube_cooldown.library.materials.load import load_material

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
