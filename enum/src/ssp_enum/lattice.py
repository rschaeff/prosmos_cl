"""Hexagonal lattice over Z = {-1, 0, 1} (Chitturi 2016 Fig. 1).

We use axial integer coordinates (q, r) for the hexagonal XY plane plus a
Z layer, yielding 3-tuples (q, r, z). With axial coords, the six in-plane
neighbors of (q, r) are:

    (q+1, r), (q-1, r), (q, r+1), (q, r-1), (q+1, r-1), (q-1, r+1)

The Z axis carries layer membership. The paper limits SSPs to "up to three
layers"; we use z in {-1, 0, 1}. Cross-layer adjacency: a point at (q, r, z)
is adjacent to (q, r, z+1) and (q, r, z-1) — vertically stacked positions in
adjacent layers. (The 2012 binary's `enumerateNeighborPositions(bool flag)`
suggests two neighbor flavors; we conservatively treat in-plane and
cross-layer as the only adjacency relations and revisit when validation
against CG-2012 outputs identifies discrepancies.)

Direction (up/down along Z) is a property of an SSE *placed* at a lattice
point, not of the lattice point itself. See `node.py`.
"""

from __future__ import annotations

from dataclasses import dataclass

# Axial-coord neighbor offsets in the hex XY plane.
_IN_PLANE_OFFSETS: tuple[tuple[int, int], ...] = (
    (1, 0), (-1, 0), (0, 1), (0, -1), (1, -1), (-1, 1),
)


@dataclass(frozen=True)
class LatticePoint:
    """A point on the hexagonal lattice; immutable.

    Axial XY coords plus layer index Z. The hexagonal in-plane structure is
    encoded in the neighbor function rather than in the coordinates.
    """

    q: int
    r: int
    z: int = 0

    def in_plane_neighbors(self) -> tuple["LatticePoint", ...]:
        """The 6 in-plane hex neighbors at the same Z layer."""
        return tuple(
            LatticePoint(self.q + dq, self.r + dr, self.z)
            for dq, dr in _IN_PLANE_OFFSETS
        )

    def cross_layer_neighbors(self) -> tuple["LatticePoint", ...]:
        """Vertically stacked neighbors in adjacent Z layers (±1)."""
        return tuple(
            LatticePoint(self.q, self.r, self.z + dz)
            for dz in (-1, 1)
            if -1 <= self.z + dz <= 1
        )

    def neighbors(self) -> tuple["LatticePoint", ...]:
        """All adjacent lattice points: 6 in-plane + 0–2 cross-layer."""
        return self.in_plane_neighbors() + self.cross_layer_neighbors()

    def is_adjacent(self, other: "LatticePoint") -> bool:
        return other in self.neighbors()
