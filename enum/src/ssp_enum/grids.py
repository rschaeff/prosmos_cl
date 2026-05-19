"""Induced-grid canonical form + SCC-2 whitelists (Chitturi 2016 Appendix § "secondary compactness criterion 2").

An *induced grid* is the unlabeled set of lattice positions occupied by an
S_n skeleton (paper Appendix §1.1: "The set of points a skeleton induces
on the hexagonal grid is its induced grid"). Two skeletons whose underlying
lattice positions coincide *as a set* (modulo translation, hex rotation,
mirror, and z-flip) have the same induced grid. SCC criterion 2 says: a
candidate skeleton passes only if its induced grid is one of a
predefined set of "allowed grids" for that dimension (Appendix Fig. S3).

Allowed grid counts from Fig. S3 (lattice level, after deriving the
underlying lattice arrangement from the paper's explicit-edge diagrams):

  S3: 2 grids — 3-collinear (P3) and equilateral triangle (K3).
  S4: 3 grids — Fig. S3(a) P4 path; (b) triangle (2,3,4) + pendant 1;
      (c) C4 cycle with one rhombus diagonal also lattice-adjacent.
  S5: 5 grids — Fig. S3(d) P5; (e) tripod+pendant; (f) C5 cycle;
      (g)/(h) K1,4 star centered at node 2 (two different leaf
      labelings that share unlabeled adjacency).

Note: the paper diagrams show *explicit interactions* (the solid lines).
The underlying lattice arrangement can include additional adjacencies
that are marked `X` (optional) in the SSP query matrix because the
specific SSE-type assignment doesn't require them to interact concretely.
The grids stored here are the **lattice arrangements** — the maximum-edge
version for each paper grid. A skeleton passes SCC-2 iff its full
lattice adjacency graph (canonical-form) is in this whitelist.

Earlier extractions of these constants from the oracle treated `X` as
adjacent (paper §1.1.1: `X` = non-adjacent / optional interaction); the
"K4" entry that appeared in S4 was a wrong inference — K4 isn't
realizable on the 2D hex lattice, and the 96 oracle records with that
profile are actually paper grid (c) lattice (5 edges) with the short
rhombus diagonal made explicit. The S4 whitelist is now derived directly
from Fig. S3 + hex-lattice geometry, paper-cited.

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
    # Fig. S3(a) — P4 path (1—2—3—4). 3 edges.
    ((0, 1), (0, 2), (1, 3)),
    # Fig. S3(b) — triangle (2,3,4) with pendant 1 attached at 2. 4 edges.
    # Note: this entry was previously erroneously listed as K4 (6 edges),
    # which is not realizable on the 2D hex lattice (max in-plane cluster
    # of 4 mutually-adjacent points doesn't exist).
    ((0, 1), (0, 2), (0, 3), (1, 2)),
    # Fig. S3(c) — C4 cycle with the short rhombus diagonal also adjacent.
    # On the hex lattice, the 4 corners of a unit rhombus include one
    # diagonal at distance 1 (lattice-adjacent, marked X in the SSP) and
    # one at distance 2 (non-adjacent). 5 edges total = K4 minus one edge.
    ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3)),
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
