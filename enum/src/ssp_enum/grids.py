"""Induced-grid canonical form + SCC-2 whitelists (Chitturi 2016 Appendix § "secondary compactness criterion 2").

An *induced grid* is the unlabeled set of lattice positions occupied by an
S_n skeleton (paper Appendix §1.1: "The set of points a skeleton induces
on the hexagonal grid is its induced grid"). Two skeletons whose underlying
lattice positions coincide *as a set* (modulo translation, hex rotation,
mirror, and z-flip) have the same induced grid. SCC criterion 2 says: a
candidate skeleton passes only if its induced grid is one of a
predefined set of "allowed grids" for that dimension (Appendix Fig. S3).

Allowed grid counts from Appendix:

  S3: 2 grids — 3-collinear and equilateral triangle.
  S4: 3 grids — Fig. S3(a, b, c). (a)=K4 with all-pair interactions;
      (b)=K4 minus one edge; (c)=4 nodes with only peripheral edges.
  S5: 5 grids — Fig. S3(d, e, f, g, h). The Appendix says Grid 5 (S3-h)
      has the same handedness specification as Grid 4 (S3-g); g and h
      share the same *unlabeled adjacency graph* (i.e., same set of
      edges up to vertex relabeling) but differ in their visual lattice
      embedding — a mirror-image pair that the lattice-symmetry group
      can't unify. Under the canonical form used here (unlabeled
      adjacency graph), g and h collapse to one signature, so we
      observe **4** S5 signatures from oracle data even though paper
      counts 5 grids.

Canonical-form: we represent an induced grid by the canonical
adjacency-graph signature (lex-min frozenset of edges over all N!
vertex relabelings). This is a coarser quotient than the paper's
lattice-embedding equivalence (it collapses g and h at S5), but it is
*sufficient* for SCC-2 — every skeleton with a g-or-h embedding has
the same adjacency graph, so accept/reject decisions are identical.

The whitelists below were derived empirically from the CG-2012 oracle
(`reference/IA-S{4,5}.txt`) and verified against the per-grid
handedness lists in the paper Appendix. They are paper-traceable but
the derivation script also lives in `tests/test_grids.py` so the
constants stay synchronized with the oracle.
"""

from __future__ import annotations

from itertools import permutations
from typing import Iterable

from .skeleton import Skeleton


GridSignature = tuple[tuple[int, int], ...]


def unlabeled_grid_signature(skel: Skeleton) -> GridSignature:
    """Canonical unlabeled adjacency-graph signature of `skel`.

    Two skeletons whose induced-grid adjacency structures are isomorphic
    (as undirected graphs, ignoring sequence labels) produce the same
    signature. The signature is the lex-smallest frozenset of (i, j)
    edges (i<j) over all N! relabelings — a standard graph
    canonical-form construction.
    """
    m = skel.adjacency_matrix()
    n = len(m)
    edges = [(i, j) for i in range(n) for j in range(i + 1, n) if m[i][j]]
    best: GridSignature | None = None
    for perm in permutations(range(n)):
        relabeled = tuple(
            sorted(tuple(sorted((perm[i], perm[j]))) for i, j in edges)
        )
        if best is None or relabeled < best:
            best = relabeled
    return best if best is not None else ()


# Whitelists keyed by dimension. Each value is a frozenset of canonical
# unlabeled-adjacency-graph signatures. A skeleton at dim N passes SCC-2
# iff its signature is in WHITELIST[N].
#
# S3 whitelist: derived from the four S3 base shapes (Linear, Triangle,
# Bent-A, Bent-B). Under unlabeled adjacency, Linear/Bent-A/Bent-B all
# collapse to the path P3 (two edges in a chain); Triangle is K3.
#
# S4 whitelist: derived from oracle IA-S4 (155 unique (skel_id, third)
# tuples → 3 signatures).
#
# S5 whitelist: derived from oracle IA-S5 (2377 unique (skel_id, third)
# tuples → 4 signatures). See module docstring re: the 4-vs-5 paper
# discrepancy (Appendix Fig. S3 grids g and h collapse).

WHITELIST_S3: frozenset[GridSignature] = frozenset({
    ((0, 1), (0, 2)),                # P3: linear / bent (any sequence labeling)
    ((0, 1), (0, 2), (1, 2)),        # K3: triangle
})

WHITELIST_S4: frozenset[GridSignature] = frozenset({
    # Fig. S3(c): 4 nodes with peripheral edges only (3 edges, P4 path).
    ((0, 1), (0, 2), (1, 3)),
    # Fig. S3(b): K4 minus one edge (5 edges).
    ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3)),
    # Fig. S3(a): K4 (all 6 pairs adjacent).
    ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)),
})

WHITELIST_S5: frozenset[GridSignature] = frozenset({
    # Sparsest: 4 edges, P5 path.
    ((0, 1), (0, 2), (1, 3), (2, 4)),
    # 7 edges (K5 minus 3): Fig. S3 most likely (f) given handedness density.
    ((0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (2, 4)),
    # K5 minus 1 edge (9 edges).
    ((0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (1, 4), (2, 3), (2, 4)),
    # K5 (all 10 pairs adjacent).
    ((0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (1, 4),
     (2, 3), (2, 4), (3, 4)),
})


WHITELIST: dict[int, frozenset[GridSignature]] = {
    3: WHITELIST_S3,
    4: WHITELIST_S4,
    5: WHITELIST_S5,
}


def is_in_whitelist(skel: Skeleton) -> bool:
    """True iff `skel`'s induced grid is in the SCC-2 whitelist for its dim.

    Dimensions outside {3, 4, 5} return True (no SCC-2 constraint defined
    yet). At dim ≤ 2, the lattice arrangement is trivially unique and
    no SCC-2 check applies.
    """
    n = skel.dim
    if n < 3:
        return True
    whitelist = WHITELIST.get(n)
    if whitelist is None:
        return True
    return unlabeled_grid_signature(skel) in whitelist
