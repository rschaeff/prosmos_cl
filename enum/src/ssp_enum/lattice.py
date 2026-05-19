"""Hexagonal lattice over Z = {-1, 0, 1} (Chitturi 2016 Fig. 1).

Each lattice point is a 3-tuple (x, y, z). The XY plane is hexagonal:
each point has six in-plane neighbors. The Z axis carries layer membership
({-1, 0, 1}) — the paper limits SSPs to "up to three layers".

A lattice point is *directed* once an SSE is placed at it: the SSE points
either up (+Z) or down (-Z). Paper convention: consecutive SSEs in the
sequence travel in opposite directions ("up-down topology").

This module is the geometry-only foundation: it defines points, neighbors,
adjacency, and the predicates that compactness checks (compactness.py) call.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class LatticePoint:
    """A point on the hexagonal lattice; immutable.

    Coordinates are integer triples; the hexagonal in-plane structure is
    encoded in the neighbor function rather than in the coordinates.
    """

    x: int
    y: int
    z: int

    def neighbors(self) -> tuple["LatticePoint", ...]:
        """Return the up-to-8 lattice neighbors (6 in-plane + up to 2 cross-layer).

        Not yet implemented. Cross-reference: paper Fig. 1a, panel (a)
        showing the 6 in-plane positions adjacent to an S2 skeleton.
        """
        raise NotImplementedError
