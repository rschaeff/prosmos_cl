"""Model B (combinatorial / declared-edges) tests.

Phase 1 (`Skeleton.edges` field) is the load-bearing piece for downstream
design-target work — these tests pin the data-structure semantics.
Phase 2 (`enumerate_skeletons_combinatorial`) is experimental; we pin
its current counts as regression markers but they reflect known
over-enumeration that future iterations should narrow.
"""

from __future__ import annotations

from ssp_enum.combine_b import (
    combine_b_single_node,
    enumerate_skeletons_combinatorial,
)
from ssp_enum.lattice import LatticePoint
from ssp_enum.skeleton import Skeleton


def test_skeleton_edges_field_overrides_geometric():
    """When `edges` is set, adjacency_matrix uses it — not geometry."""
    pts = (LatticePoint(0, 0, 0), LatticePoint(1, 0, 0), LatticePoint(2, 0, 0))
    geometric = Skeleton(points=pts)  # edges=None → use lattice positions
    declared_sparse = Skeleton(points=pts, edges=frozenset({(0, 1)}))
    # Geometric: 1-2 adj, 2-3 adj, 1-3 not adj. So matrix row 0 = (F, T, F).
    assert geometric.adjacency_matrix()[0][1] is True   # geometric: lattice adj
    assert geometric.adjacency_matrix()[1][2] is True
    # Declared: only edge (0,1) — row 0 col 1 True, row 1 col 2 False
    assert declared_sparse.adjacency_matrix()[0][1] is True
    assert declared_sparse.adjacency_matrix()[1][2] is False


def test_skeleton_edges_default_is_none():
    """Backward compat: existing constructors yield edges=None."""
    s = Skeleton(points=(LatticePoint(0, 0, 0),))
    assert s.edges is None


def test_model_b_s3_count():
    """Phase 2 regression pin: Model B S3 = 19 — many lattice variants of
    the same declared P3 path graph. CG-2012 oracle S3 = 11. The
    over-count reflects per-lattice canonical-key dedup which is finer
    than CG-2012's graph-isomorphism dedup."""
    assert len(enumerate_skeletons_combinatorial(3)) == 19


def test_model_b_reaches_k14_star_at_s5():
    """Phase 2 reach check: Model B produces K1,4 star declared graphs,
    which Model A (combine.py) cannot. The full S5 count is over-
    enumerated (4842 vs oracle 648) but the declared-graph coverage
    matches what oracle has."""
    from itertools import permutations
    K14 = ((0, 1), (0, 2), (0, 3), (0, 4))

    def canon_edges(edges, n):
        best = None
        for perm in permutations(range(n)):
            rel = tuple(sorted(tuple(sorted((perm[i], perm[j]))) for i, j in edges))
            if best is None or rel < best:
                best = rel
        return best

    skels = enumerate_skeletons_combinatorial(5)
    k14_count = sum(1 for s in skels if canon_edges(s.edges, 5) == K14)
    assert k14_count > 0


def test_model_b_s5_count_regression():
    """Phase 2 regression pin: 4842 — used as a marker. Future iterations
    should bring this down via graph-isomorphism canonicalization."""
    assert len(enumerate_skeletons_combinatorial(5)) == 4842


def test_combine_b_single_node_declares_one_edge():
    """Each candidate from `combine_b_single_node` adds exactly one
    declared edge to s1's edge set."""
    s2 = Skeleton(
        points=(LatticePoint(0, 0, 0), LatticePoint(1, 0, 0)),
        edges=frozenset({(0, 1)}),
    )
    candidates = list(combine_b_single_node(s2))
    for c in candidates:
        n_new = len(c.edges) - len(s2.edges)
        assert n_new == 1, f"expected 1 new edge, got {n_new}: {c.edges}"
