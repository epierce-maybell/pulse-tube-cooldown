"""Exchange-gas requirements: pressure from the vacuum-can cubic and gap conductance over facing surfaces."""

import warnings

import numpy as np
import pytest
from scipy.constants import R

from pulse_tube_cooldown.core.errors import ConfigurationError
from pulse_tube_cooldown.core.types import Gas
from pulse_tube_cooldown.library.gases.load import load_gas
from pulse_tube_cooldown.library.links.exchange_gas import ExchangeGas, FacingSurfaces
from pulse_tube_cooldown.library.links.load import load_link

# monitor_ovc.yaml's mbar coefficients times 100, in Pa / K^3, Pa / K^2, Pa / K, and Pa.
_MONITOR_OVC_COEFFICIENTS = (3.16030871e-07, -1.32989499e-04, 1.90126130e-02, -7.78545808e-01)


def _gas():
    # At 10 K the molecular conductance per area equals pressure numerically; conductivity is 0.02 W/(m K).
    return Gas(
        name="Worked example",
        molar_mass=R / (80 * np.pi),
        heat_capacity_ratio=3.0,
        accommodation_coeff=lambda temperature: np.full_like(temperature, 0.5),
        thermal_conductivity=lambda temperature: np.full_like(temperature, 0.02),
        saturation_pressure=lambda temperature: np.full_like(temperature, np.inf),
        adsorption_energy=0.0,
        adsorption_constant=0.0,
        adsorption_capacity=0.0,
    )


def _chain(pressure_coefficients, gas=None):
    return ExchangeGas(
        name="Exchange gas",
        nodes=("Can", "Middle", "Inner"),
        gas=gas or _gas(),
        pressure_node="Middle",
        pressure_coefficients=pressure_coefficients,
        surfaces=(
            FacingSurfaces(nodes=("Can", "Middle"), area=2.0, gap=0.01),
            FacingSurfaces(nodes=("Middle", "Inner"), area=3.0, gap=0.02),
        ),
    )


def test_pressure_is_the_cubic_at_the_pressure_node():
    # At 2 K: 1 * 8 - 2 * 4 + 3 * 2 - 4 = 2 Pa.
    link = _chain((1.0, -2.0, 3.0, -4.0))

    assert link.calc_pressure(np.array([50.0, 2.0, 70.0])) == pytest.approx(2.0)


def test_negative_cubic_clamps_pressure_to_zero():
    # At 1 K: 1 - 2 + 3 - 4 = -2 Pa.
    link = _chain((1.0, -2.0, 3.0, -4.0))

    assert link.calc_pressure(np.array([50.0, 1.0, 70.0])) == 0.0


@pytest.mark.parametrize(("pressure_temperature", "expected"), [(300.0, 1.48902), (69.0, 3.9805e-3), (67.0, 0.0)])
def test_monitor_ovc_cubic_shuts_off_near_68_kelvin(pressure_temperature, expected):
    # Evaluated term by term from monitor_ovc.yaml: 0.0148902 mbar at 300 K, 3.9805e-5 mbar at 69 K, -6.64e-5 mbar
    # at 67 K.
    link = _chain(_MONITOR_OVC_COEFFICIENTS)

    pressure = link.calc_pressure(np.array([300.0, pressure_temperature, 300.0]))

    assert pressure == pytest.approx(expected, rel=1e-4, abs=0.0)


@pytest.mark.parametrize(
    ("pressure", "expected"),
    [
        # Free-molecular: conductance = area * pressure, 2e-12 and 3e-12 W/K.
        (1e-12, [8e-12, -2e-11, 1.2e-11]),
        # Viscous: conductance = area * k / gap, 4 and 3 W/K.
        (1e12, [16.0, -28.0, 12.0]),
    ],
)
def test_flows_sum_gap_conductance_over_facing_pairs(pressure, expected):
    # Both pairs have a mean temperature of 10 K; differences are +4 K and -4 K across them.
    flow = _chain((0.0, 0.0, 0.0, pressure)).calc_heat_flow(np.array([8.0, 12.0, 8.0]))

    np.testing.assert_allclose(flow, expected, rtol=1e-9)
    assert flow.sum() == pytest.approx(0.0, abs=1e-12 * np.abs(flow).max())


def test_shut_off_gas_carries_no_heat_below_its_property_range():
    link = _chain(_MONITOR_OVC_COEFFICIENTS, gas=load_gas("N2"))

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        flow = link.calc_heat_flow(np.array([290.0, 60.0, 40.0]))

    np.testing.assert_array_equal(flow, [0.0, 0.0, 0.0])


_CONFIG = {
    "name": "Exchange gas",
    "type": "exchange_gas",
    "nodes": ["Vacuum can", "PRP", "CFP"],
    "gas": "N2",
    "pressure_node": "CFP",
    "pressure_coefficients": list(_MONITOR_OVC_COEFFICIENTS),
    "surfaces": [
        {"nodes": ["Vacuum can", "PRP"], "area": 0.5, "gap": 0.02},
        {"nodes": ["PRP", "CFP"], "area": 0.3, "gap": 0.01},
    ],
}


def test_loads_registered_exchange_gas_with_its_gas():
    link = load_link(_CONFIG)

    assert isinstance(link, ExchangeGas)
    assert (link.name, link.nodes, link.pressure_node) == ("Exchange gas", ("Vacuum can", "PRP", "CFP"), "CFP")
    assert link.gas.name == "N2"
    assert link.pressure_coefficients == _MONITOR_OVC_COEFFICIENTS
    assert link.surfaces == (
        FacingSurfaces(nodes=("Vacuum can", "PRP"), area=0.5, gap=0.02),
        FacingSurfaces(nodes=("PRP", "CFP"), area=0.3, gap=0.01),
    )


def test_pressure_node_outside_the_link_is_named():
    with pytest.raises(ConfigurationError, match=r"(?s)Exchange gas.*STP"):
        load_link(_CONFIG | {"pressure_node": "STP"})


def test_surface_node_outside_the_link_is_named():
    surfaces = [*_CONFIG["surfaces"], {"nodes": ["CFP", "STP"], "area": 0.2, "gap": 0.01}]

    with pytest.raises(ConfigurationError, match=r"(?s)Exchange gas.*STP"):
        load_link(_CONFIG | {"surfaces": surfaces})
