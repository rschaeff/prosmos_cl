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

    Adjacency = any non-blank, non-'*', non-'-' cell on the upper triangle.
    """
    from itertools import permutations
    n = len(matrix)
    edges = []
    for i in range(n):
        for j in range(i + 1, n):
            cell = matrix[i][j].strip()
            if cell and cell not in ("*", "-"):
                edges.append((i, j))
    best = None
    for perm in permutations(range(n)):
        relabeled = tuple(
            sorted(tuple(sorted((perm[i], perm[j]))) for i, j in edges)
        )
        if best is None or relabeled < best:
            best = relabeled
    return best if best is not None else ()


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


def test_oracle_s4_every_record_in_whitelist():
    sigs = _oracle_signatures(REFERENCE / "IA-S4.txt")
    off = [s for s in sigs if s not in WHITELIST_S4]
    assert not off, f"{len(off)} S4 oracle records have off-whitelist grids"
    # And the empirical set of signatures matches the whitelist exactly
    assert set(sigs) == WHITELIST_S4


def test_oracle_s5_every_record_in_whitelist():
    sigs = _oracle_signatures(REFERENCE / "IA-S5.txt")
    off = [s for s in sigs if s not in WHITELIST_S5]
    assert not off, f"{len(off)} S5 oracle records have off-whitelist grids"
    assert set(sigs) == WHITELIST_S5


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
