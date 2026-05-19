"""SCC-2 whitelist tests.

The whitelist constants in `grids.py` are derived empirically from the
CG-2012 oracle and cross-checked against Chitturi 2016 Appendix Fig. S3.
These tests pin both ends: every oracle record's induced grid must be
in the whitelist, and the whitelist size must match the paper-stated
count of allowed grids (modulo the documented g≡h collapse at S5).
"""

from __future__ import annotations

from pathlib import Path

import pytest

from ssp_enum.compactness import is_compact, passes_scc_2
from ssp_enum.enumerate import enumerate_dim
from ssp_enum.grids import (
    WHITELIST_S3,
    WHITELIST_S4,
    WHITELIST_S5,
    is_in_whitelist,
    unlabeled_grid_signature,
)
from ssp_enum.lattice import LatticePoint
from ssp_enum.oracle import parse
from ssp_enum.skeleton import Skeleton

REFERENCE = Path(__file__).parent.parent / "reference"


def _matrix_to_signature(matrix):
    """Reproduce unlabeled_grid_signature() from an oracle interaction matrix.

    Lattice adjacency = explicit interaction cells {c, t, u, v, C, T}.
    `X` (uppercase) and `-` are non-lattice-adjacent (paper §1.1.1:
    "Non adjacent SSEs can either interact optionally (X) or not interact
    at all (-)"). `*` is the diagonal.

    Note: oracle records only show the *explicit-interaction* projection
    of the lattice; the full lattice may have additional adjacencies
    marked X. So this signature is a *subgraph* of the true lattice
    grid — callers should account for that when comparing to WHITELIST.
    """
    from itertools import permutations
    n = len(matrix)
    edges = []
    for i in range(n):
        for j in range(i + 1, n):
            cell = matrix[i][j].strip()
            if cell and cell not in ("*", "-", "X"):
                edges.append((i, j))
    best = None
    for perm in permutations(range(n)):
        relabeled = tuple(
            sorted(tuple(sorted((perm[i], perm[j]))) for i, j in edges)
        )
        if best is None or relabeled < best:
            best = relabeled
    return best if best is not None else ()


def _is_subgraph_of_any(small_sig, whitelist):
    """True iff `small_sig`'s edge set is a subset of some `whitelist`
    entry under graph isomorphism (permute small_sig's labels to match)."""
    from itertools import permutations
    small_edges = set(small_sig)
    if not small_edges:
        return True
    n_small = max(max(e) for e in small_edges) + 1
    for big_sig in whitelist:
        big_edges = set(big_sig)
        # Distinct vertices in big_sig
        big_vertices = set()
        for e in big_edges:
            big_vertices.update(e)
        if len(big_vertices) < n_small:
            continue
        # Try all relabelings of small_sig's vertices into big_sig's vertex set
        big_vlist = sorted(big_vertices)
        # Iterate over injective maps from {0..n_small-1} to big_vertices
        from itertools import permutations as _p
        for chosen in _p(big_vlist, n_small):
            relabeled = {tuple(sorted((chosen[i], chosen[j]))) for i, j in small_edges}
            if relabeled.issubset(big_edges):
                return True
    return False


def _oracle_signatures(path: Path):
    if not path.exists() or path.stat().st_size == 0:
        pytest.skip(f"oracle file {path.name} missing/empty in this checkout")
    return [_matrix_to_signature(r.matrix) for r in parse(path)]


def test_whitelist_s3_sizes():
    """Paper Appendix § Size 3: two grids — 3 collinear, equilateral triangle."""
    assert len(WHITELIST_S3) == 2


def test_whitelist_s4_size():
    """Paper Appendix Fig. S3(a-c): 3 induced grids at S4."""
    assert len(WHITELIST_S4) == 3


def test_whitelist_s5_size():
    """Paper Appendix Fig. S3(d-h) lists 5 grids; (g) and (h) share the
    same unlabeled adjacency graph (Appendix lists Grid 5 handedness as
    'Same as Grid 4 (S3-g)'), so under unlabeled-graph canonical form
    they collapse to 1 signature → 4 total at S5."""
    assert len(WHITELIST_S5) == 4


def test_oracle_s4_every_record_subgraph_of_whitelist():
    """Every S4 oracle record's *explicit-interaction* grid (cells in
    {c,t,u,v,C,T}, excluding X) must be a subgraph of some WHITELIST_S4
    entry (which represents the full lattice grid including X-marked
    optional interactions). This is the right semantic since oracle
    records show SSPs (explicit) and WHITELIST stores lattice grids."""
    sigs = _oracle_signatures(REFERENCE / "IA-S4.txt")
    off = [s for s in sigs if not _is_subgraph_of_any(s, WHITELIST_S4)]
    assert not off, (
        f"{len(off)} S4 oracle records have explicit-edge grids that don't "
        f"sit under any WHITELIST_S4 lattice grid as a subgraph"
    )


def test_oracle_s5_every_record_subgraph_of_whitelist():
    """Same semantic at S5. NOTE: WHITELIST_S5 still encodes the
    *X-included* extraction from oracle (Phase A vintage) — the
    Fig. S3-derived lattice whitelist for S5 is pending. This test
    currently passes by happy accident for the simpler grids."""
    sigs = _oracle_signatures(REFERENCE / "IA-S5.txt")
    off = [s for s in sigs if not _is_subgraph_of_any(s, WHITELIST_S5)]
    # Allow up to a small number of off-whitelist records as TODO marker
    # for the pending WHITELIST_S5 correction.
    assert len(off) <= 0, (
        f"{len(off)} S5 oracle records have explicit-edge grids that don't "
        f"sit under any WHITELIST_S5 lattice grid"
    )


def test_s3_enumeration_passes_scc_2():
    """All 11 S3 SSPs we produce must satisfy SCC-2 (P3 or K3 grid)."""
    for ssp in enumerate_dim(3):
        assert passes_scc_2(ssp), f"{ssp} fails SCC-2"
        assert is_compact(ssp), f"{ssp} fails full compactness"


def test_s4_off_whitelist_arrangement_rejected():
    """A 4-node skeleton whose adjacency graph isn't in WHITELIST_S4
    must fail SCC-2. We construct one by taking three K3-arranged points
    plus a 4th point adjacent to only one of them — passes PCC, doesn't
    extend a line, doesn't get ≥2 neighbors. (This will also fail SCC-1,
    but the SCC-2 check is what we want to exercise.)"""
    skel = Skeleton(points=(
        LatticePoint(0, 0, 0),
        LatticePoint(1, 0, 0),
        LatticePoint(0, 1, 0),
        LatticePoint(5, 5, 0),  # disconnected — won't appear in any S4 grid
    ))
    # The unlabeled adjacency graph of this skeleton: edges (0,1), (0,2),
    # (1,2) only (the disconnected 4th has no edges). That's K3 + isolated,
    # which is not in WHITELIST_S4 (all S4 grids are connected).
    assert not passes_scc_2(skel)
