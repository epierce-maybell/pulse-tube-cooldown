"""Test the archetypal pulse-tube cooler against test-stand data and physical invariants."""

import pytest

from tests.acceptance.fixtures.pulse_tube_cooler.compile_records import compile_records


@pytest.fixture
def objectives():
    """Return one-compressor measurement objectives and absolute tolerances."""
    return compile_records()[1]
