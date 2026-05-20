"""Tests for the ProSMoS query writer."""

from __future__ import annotations

from pathlib import Path

import pytest

from ssp_enum.oracle import parse
from ssp_enum.prosmos import find_record, write_query

REFERENCE = Path(__file__).parent.parent / "reference"


def _records():
    p = REFERENCE / "IA-S5.txt"
    if not p.exists() or p.stat().st_size == 0:
        pytest.skip("oracle IA-S5.txt missing in this checkout")
    return list(parse(p))


def test_find_record_returns_matching():
    rs = _records()
    r = find_record(rs, "5-141-7-7")
    assert r.dim == 5
    assert r.skeleton_id == 141
    assert r.third_idx == 7
    assert r.sub_first == 7


def test_find_record_raises_for_missing():
    rs = _records()
    with pytest.raises(KeyError):
        find_record(rs, "5-99999-0-0")


def test_write_query_header_format():
    """First two lines are SSE indices and types."""
    rs = _records()
    r = find_record(rs, "5-141-7-7")
    out = write_query(r)
    lines = out.splitlines()
    assert lines[0] == "1 2 3 4 5"
    assert lines[1] == "H E H H E"  # matches r.sse_types


def test_write_query_matrix_upper_triangular_with_indent():
    """Row i has 2*i leading spaces."""
    rs = _records()
    r = find_record(rs, "5-141-7-7")
    out = write_query(r)
    lines = out.splitlines()
    # Matrix rows are lines 2..6 (0-indexed: 2, 3, 4, 5, 6)
    for i in range(5):
        row = lines[2 + i]
        # Leading spaces: 2*i
        assert row[: 2 * i] == " " * (2 * i), f"row {i}: {row!r}"
        # First non-space character is '*'
        assert row.lstrip()[0] == "*"


def test_write_query_length_per_sse():
    """Every SSE has a length directive with the right min/max."""
    rs = _records()
    r = find_record(rs, "5-141-7-7")
    out = write_query(r)
    lengths = [ln for ln in out.splitlines() if ln.startswith("length ")]
    assert len(lengths) == r.dim
    # Spot-check: helix gets min 8, strand gets min 5.
    for ln in lengths:
        # length <i> <T> <min> <max>
        parts = ln.split()
        idx = int(parts[1])
        typ = parts[2]
        mn = int(parts[3])
        expected_type = r.sse_types[idx - 1]
        assert typ == expected_type
        assert mn == (8 if typ == "H" else 5)


def test_write_query_sheet_directives():
    """5-141-7-7 has β-hairpin (E2-E5) → one sheetS line."""
    rs = _records()
    r = find_record(rs, "5-141-7-7")
    out = write_query(r)
    sheet_lines = [ln for ln in out.splitlines() if ln.startswith("sheet")]
    # Exactly one sheetS line for {2, 5}.
    assert sheet_lines == ["sheetS 2 5"]


def test_write_query_handedness_directives():
    """5-141-7-7 has 9 handedness triples — all should appear."""
    rs = _records()
    r = find_record(rs, "5-141-7-7")
    out = write_query(r)
    hand_lines = [ln for ln in out.splitlines() if ln.startswith("handedness ")]
    assert len(hand_lines) == 9
    # Spot-check one: oracle has (1, 2, 4, 'R') based on prior inspection
    # of the same record; assert that line is present.
    assert "handedness 1 2 4 R" in hand_lines


def test_write_query_group_i_target_5_269_0_0():
    """Sanity check on a group-(i) design target: all-strand β-sheet."""
    rs = _records()
    r = find_record(rs, "5-269-0-0")
    out = write_query(r)
    lines = out.splitlines()
    assert lines[0] == "1 2 3 4 5"
    assert lines[1] == "E E E E E"
    # Length lines: all should be E with min 5.
    length_lines = [ln for ln in lines if ln.startswith("length ")]
    assert len(length_lines) == 5
    for ln in length_lines:
        parts = ln.split()
        assert parts[2] == "E"
        assert int(parts[3]) == 5


def test_oracle_design_targets_reachable():
    """9 of 14 design targets are in CG-2012's IA-S5.txt. The other 5
    are from a later (post-2012) enumeration on the website. Pin which."""
    rs = _records()
    reachable = []
    unreachable = []
    targets = [
        "5-269-0-0", "5-311-0-0", "5-289-0-0", "5-288-0-0",
        "5-280-0-0", "5-282-0-0", "5-306-0-0", "5-309-0-0",
        "5-283-1-2", "5-307-1-2", "5-243-1-2",
        "5-265-7-7", "5-234-7-7", "5-141-7-7",
    ]
    for t in targets:
        try:
            find_record(rs, t)
            reachable.append(t)
        except KeyError:
            unreachable.append(t)
    assert len(reachable) == 9
    assert len(unreachable) == 5
    # The 5 unreachable are the ones with `1-2` or `7-7` suffix (except 5-141-7-7).
    assert "5-283-1-2" in unreachable
    assert "5-141-7-7" in reachable
