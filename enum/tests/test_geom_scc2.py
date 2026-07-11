"""Geometric SCC-2 (Fig-S3 hex-congruence) — Step A of the root fix.

Locks the paper-faithful S5 count (140), the per-grid breakdown, the
subset relation to the graph-based whitelist, and the fact that the 58
dropped skeletons are exactly the grid-e-graph bent variants.
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from ssp_enum.compactness import set_scc2_mode
from ssp_enum.enumerate import enumerate_skeletons
from ssp_enum.combine import canonical_key
from ssp_enum.grids import unlabeled_grid_signature
from ssp_enum.geom_scc2 import (
    geometric_canonical, grid_of, FIG_S3_S5, GEOM_WHITELIST,
)


def test_figs3_grids_have_expected_edge_counts():
    from ssp_enum.lattice import LatticePoint
    expected = {"d": 4, "e": 5, "f": 7, "gh": 6}
    for name, pts in FIG_S3_S5.items():
        P = [LatticePoint(q, r, 0) for q, r in pts]
        e = sum(1 for i in range(len(P)) for j in range(i + 1, len(P))
                if P[i].is_adjacent(P[j]))
        assert e == expected[name], (name, e)


def test_four_grids_are_geometrically_distinct():
    assert len(GEOM_WHITELIST[5]) == 4


def test_geometric_canonical_is_isometry_invariant():
    # grid e, rotated + reflected, must canonicalize identically.
    base = FIG_S3_S5["e"]
    rot = [(-r, q + r) for q, r in base]
    refl = [(q + r, -r) for q, r in base]
    shifted = [(q + 7, r - 3) for q, r in base]
    c = geometric_canonical(base)
    assert geometric_canonical(rot) == c
    assert geometric_canonical(refl) == c
    assert geometric_canonical(shifted) == c


def _reset():
    set_scc2_mode("graph")


def test_graph_mode_is_198_geometric_is_140():
    try:
        set_scc2_mode("graph")
        assert len(enumerate_skeletons(5)) == 198
        set_scc2_mode("geometric")
        assert len(enumerate_skeletons(5)) == 140
    finally:
        _reset()


def test_geometric_per_grid_breakdown():
    try:
        set_scc2_mode("geometric")
        from collections import Counter
        by = Counter(grid_of(s) for s in enumerate_skeletons(5))
        assert dict(by) == {"f": 60, "e": 40, "gh": 28, "d": 12}
        assert None not in by  # every survivor maps to a Fig-S3 grid
    finally:
        _reset()


def test_geometric_is_strict_subset_and_drops_are_grid_e_graph():
    try:
        set_scc2_mode("graph")
        graph = {canonical_key(s): s for s in enumerate_skeletons(5)}
        set_scc2_mode("geometric")
        geom = {canonical_key(s) for s in enumerate_skeletons(5)}
        assert geom <= set(graph)
        dropped = [graph[k] for k in set(graph) - geom]
        assert len(dropped) == 58
        # every dropped skeleton lives in grid-e's 5-edge adjacency graph
        assert all(len(unlabeled_grid_signature(s)) == 5 for s in dropped)
        # ...and none is hex-congruent to a Fig-S3 grid
        assert all(grid_of(s) is None for s in dropped)
    finally:
        _reset()


def test_s3_s4_unaffected_by_mode():
    try:
        for n in (3, 4):
            set_scc2_mode("graph")
            a = len(enumerate_skeletons(n))
            set_scc2_mode("geometric")
            b = len(enumerate_skeletons(n))
            assert a == b, (n, a, b)
    finally:
        _reset()
