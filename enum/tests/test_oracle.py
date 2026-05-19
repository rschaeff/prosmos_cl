"""Tests against Chalam's CG-2012 IA.txt outputs.

These are the green/red signal for our enumeration's correctness.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from ssp_enum.oracle import parse, SSPRecord

REFERENCE = Path(__file__).parent.parent / "reference"


def _records(name: str) -> list[SSPRecord]:
    path = REFERENCE / name
    if not path.exists() or path.stat().st_size == 0:
        pytest.skip(f"oracle file {name} is missing or empty in this checkout")
    return list(parse(path))


def test_s4_count_matches_cg2012():
    """CG-2012 produces 203 IA.txt blocks at dimension 4 (per `grep -cE` count)."""
    records = _records("IA-S4.txt")
    assert len(records) == 203


def test_s5_count_matches_cg2012():
    """CG-2012 produces 2807 IA.txt blocks at dimension 5 (1472 unique skeletons)."""
    records = _records("IA-S5.txt")
    assert len(records) == 2807


def test_s5_first_record_shape():
    """The first S5 record should be [5-0-0] [0-0], all-E, antiparallel meander."""
    records = _records("IA-S5.txt")
    first = records[0]
    assert (first.dim, first.skeleton_id, first.third_idx) == (5, 0, 0)
    assert (first.sub_first, first.sub_second) == (0, 0)
    assert first.sse_types == ("E", "E", "E", "E", "E")
    # The matrix has 't' on the four sequence-adjacent pairs (β-meander)
    m = first.matrix
    assert m[0][1] == "t"  # 1-2 antiparallel H-bond
    assert m[1][2] == "t"  # 2-3
    assert m[2][3] == "t"  # 3-4
    assert m[3][4] == "t"  # 4-5
    # All five strands in one sheet
    assert first.same_sheet == ((1, 2, 3, 4, 5),)


def test_s5_design_target_141_7_7():
    """The 14th design target (5-141-7-7) is a group-vi pattern (3 H + 2 E)."""
    records = _records("IA-S5.txt")
    target = next(
        r for r in records
        if (r.skeleton_id, r.third_idx, r.sub_first, r.sub_second) == (141, 7, 7, 1004)
    )
    assert target.sse_types == ("H", "E", "H", "H", "E")
    assert target.same_sheet == ((2, 5),)  # β-hairpin only
    # The two β-strands are antiparallel H-bonded
    assert target.matrix[1][4] == "t"
