"""Models cooling power from the pulse tube head."""

__all__ = ["PulseTubeHead"]

from dataclasses import dataclass

from _pytest.nodes import Node

@dataclass(frozen=True, slots=True, kw_only=True)
class HeadParameters:
    """Parameters fit for a specific pulse tube head."""


@dataclass(frozen=True, slots=True, kw_only=True)
class RegeneratorSegment:
    """One piece of the regenerator matrix, modeled as a Node."""
    node: Node
    span: tuple[float, float]

@dataclass(frozen=True, slots=True, kw_only=True)
class PulseTubeHead:
    """The pulse tube head, including both stage-1 and stage-2 pulse tubes and the regenerator.

    Attributes:
        name: Unique link identifier within the graph
    """
    name: str
    stages: tuple[str, str]
    segments: tuple[RegeneratorSegment, ...]
    clamp: str | None
    clamp_position: float
    contact_resistance: float
