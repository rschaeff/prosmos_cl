"""Skeleton = an ordered tuple of LatticePoints with adjacency-induced edges.

A skeleton represents the SSE arrangement of an SSP before SSE-type
(H/E) or interaction-type (parallel/antiparallel, H-bonded/not) assignment.
SSEs are numbered 1..N in protein sequence order; the same physical
spatial arrangement with different sequence labellings is a different
skeleton (paper Methods: "Nodes were assembled into skeletons by linking
numerically sequential SSEs in an antiparallel fashion").

Skeleton equality and hashing use the canonical-form induced by:
  - lattice symmetries (6-fold rotation, mirror) of the hex plane
  - layer-flip (z -> -z) — for *un*handed skeletons only; layer-flip is
    the operation that distinguishes L from R chirality, so chiral
    skeletons keep it out of the quotient.

so that "same arrangement up to lattice symmetry, same sequence labeling,
same chirality" deduplicates to a single skeleton.

Chirality
---------
Closed cycles (e.g., the S3 triangle) are intrinsically chiral: the
sequence walk around the cycle is CW or CCW, distinguishable by 2D
cross-product sign. Open chains have no inherent chirality and gain
L/R variants only when SSEs are vertically displaced across Z layers
(paper terminology: "handedness via Z-displacement"). In this
representation chirality is recorded as a `'L' | 'R' | None` label
that maps directly to CG-2012's `hand i j k L|R` IA.txt line and to
ProSMoS's `handedness i j k L|R` query line. The geometric realization
(specific Z assignments per node) is a separate concern, to be derived
from (lattice points + chirality) when SSE-type and direction
assignment land.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Optional

from .lattice import LatticePoint


@dataclass(frozen=True)
class Skeleton:
    """An ordered tuple of LatticePoints (sequence position -> lattice point).

    `points[i-1]` is the lattice point of SSE labeled `i` (1-based in the
    matrix format; 0-based in this tuple).

    `chirality` is one of ``None`` (unhanded / planar / no `hand` line in
    IA.txt), ``'L'``, or ``'R'``. Two skeletons with identical lattice
    points but different chirality are distinct SSPs.
    """

    points: tuple[LatticePoint, ...]
    chirality: Optional[str] = None

    @property
    def dim(self) -> int:
        return len(self.points)

    def adjacency_matrix(self) -> tuple[tuple[bool, ...], ...]:
        """Upper-triangular boolean: rows[i][j] = is_adjacent(i, j) for i<j."""
        n = self.dim
        rows: list[tuple[bool, ...]] = []
        for i in range(n):
            row = tuple(
                self.points[i].is_adjacent(self.points[j]) if j > i else False
                for j in range(n)
            )
            rows.append(row)
        return tuple(rows)

    def is_connected_by_sequence(self) -> bool:
        """Every consecutive sequence pair (i, i+1) must be lattice-adjacent.

        Paper Methods: "links model loops" between sequential SSEs, but the
        added SSE must be adjacent to ≥2 existing SSEs of the skeleton. The
        weaker "consecutive sequence are adjacent" condition is the minimum
        requirement; full compactness lives in `compactness.py`.
        """
        return all(
            self.points[i].is_adjacent(self.points[i + 1])
            for i in range(self.dim - 1)
        )
