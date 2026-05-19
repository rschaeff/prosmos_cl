"""Smoke tests against Chalam's CG-2012 IA.txt outputs.

These tests are the green/red signal for our enumeration's correctness.
Until enumerate.py exists, they're red and that's expected.
"""

from __future__ import annotations

from pathlib import Path

import pytest

REFERENCE = Path(__file__).parent.parent / "reference"


@pytest.mark.xfail(reason="oracle.parse() not yet implemented")
def test_s3_count_matches_cg2012():
    """CG-2012 produces 11 SSPs at dimension 3 (per its Info.txt)."""
    from ssp_enum.oracle import parse

    records = list(parse(REFERENCE / "IA-S3.txt"))
    assert len(records) == 11


@pytest.mark.xfail(reason="oracle.parse() not yet implemented")
def test_s4_count_matches_cg2012():
    """CG-2012 produces 203 IA.txt blocks at dimension 4."""
    from ssp_enum.oracle import parse

    records = list(parse(REFERENCE / "IA-S4.txt"))
    assert len(records) == 203


@pytest.mark.xfail(reason="oracle.parse() not yet implemented")
def test_s5_count_matches_cg2012():
    """CG-2012 produces 2807 IA.txt blocks at dimension 5 (1472 unique SSPs)."""
    from ssp_enum.oracle import parse

    records = list(parse(REFERENCE / "IA-S5.txt"))
    assert len(records) == 2807
