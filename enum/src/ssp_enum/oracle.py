"""Parser for Chalam's CG-2012 IA.txt format.

Used solely as a validation oracle. Each block in IA.txt is one SSP record:

    [5-269-0]      <- outer header: dim-skeletonId-thirdIdx
    [0-1911]       <- inner header: sub-first - sub-second (purpose: serial id)
     1 2 3 4 5     <- SSE index header
     E E E E E     <- SSE type per index (E=strand, H=helix, Z=unset)
     * - c - c     <- upper-triangular interaction matrix
       * t c -
         * - -
           * -
             *
    sS 1 3 5 2 4   <- sheetS: which SSEs share a sheet (order = spatial, R-to-L)
    sD 1 2         <- sheetD: SSE pairs explicitly in different sheets
    hand 1 2 5 L   <- handedness of SSE triple (L = left, R = right)

This module reads that format into structured Python records so we can:
  - count SSPs per dimension (validation against CG-2012 Info.txt totals)
  - look up specific SSPs by id (cross-check against our fresh enumeration)
  - never accept IA.txt as input to anything downstream (use ProSMoS queries)

This file is the *only* place in the package that consumes IA.txt.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator

_OUTER_RE = re.compile(r"^\[(\d+)-(\d+)-(\d+)\]$")
_INNER_RE = re.compile(r"^\[(\d+)-(\d+)\]$")


@dataclass
class SSPRecord:
    """One SSP entry parsed from IA.txt. Pure data, no semantics."""

    dim: int
    skeleton_id: int
    third_idx: int
    sub_first: int
    sub_second: int
    sse_types: tuple[str, ...]
    matrix: tuple[tuple[str, ...], ...]
    same_sheet: tuple[tuple[int, ...], ...] = field(default_factory=tuple)
    diff_sheet: tuple[tuple[int, ...], ...] = field(default_factory=tuple)
    handedness: tuple[tuple[int, int, int, str], ...] = field(default_factory=tuple)

    @property
    def block_id(self) -> str:
        return f"{self.dim}-{self.skeleton_id}-{self.third_idx}"

    @property
    def full_id(self) -> str:
        return f"{self.block_id}.{self.sub_first}-{self.sub_second}"


def parse(path: str | Path) -> Iterator[SSPRecord]:
    """Yield one SSPRecord per entry in an IA.txt file.

    The IA.txt format is line-oriented; we drive a small state machine.
    Records are delimited by the outer `[dim-skel-third]` header. Within
    a record, the inner `[sub_first-sub_second]` header is required and
    is followed by the SSE-index header, the types line, and the
    upper-triangular interaction matrix (one row per SSE, with leading
    indent padding the diagonal). Optional constraint lines (`sS`/`sD`/
    `hand`) follow before the next outer header or EOF.
    """
    lines = Path(path).read_text().splitlines()
    i = 0
    n = len(lines)

    while i < n:
        # Skip blanks / advance until we find an outer header
        m_outer = _OUTER_RE.match(lines[i].strip()) if lines[i].strip() else None
        if not m_outer:
            i += 1
            continue

        dim, skel, third = (int(g) for g in m_outer.groups())
        i += 1

        # Next non-blank line must be the inner header
        while i < n and not lines[i].strip():
            i += 1
        if i >= n:
            return
        m_inner = _INNER_RE.match(lines[i].strip())
        if not m_inner:
            # Malformed entry; skip
            continue
        sub_first, sub_second = (int(g) for g in m_inner.groups())
        i += 1

        # SSE index line — sanity check it's a row of integers 1..dim
        idx_line = lines[i].split() if i < n else []
        if [int(x) for x in idx_line if x.isdigit()] != list(range(1, dim + 1)):
            # Skip malformed
            continue
        i += 1

        # SSE types line
        sse_types = tuple(lines[i].split())
        if len(sse_types) != dim:
            continue
        i += 1

        # Matrix: dim rows, each upper-triangular. We read dim lines and
        # reconstruct the (dim x dim) symmetric matrix as upper-triangular
        # only (lower triangle stored as empty string by convention).
        matrix: list[list[str]] = [["" for _ in range(dim)] for _ in range(dim)]
        for row in range(dim):
            cells = lines[i + row].split()
            # cells[0] is the diagonal '*', then upper-triangular entries
            for col_offset, cell in enumerate(cells):
                col = row + col_offset
                if col < dim:
                    matrix[row][col] = cell
        i += dim

        # Constraint lines until next outer header or blank/EOF
        same_sheet: list[tuple[int, ...]] = []
        diff_sheet: list[tuple[int, ...]] = []
        handedness: list[tuple[int, int, int, str]] = []
        while i < n and lines[i].strip() and not _OUTER_RE.match(lines[i].strip()):
            line = lines[i].strip()
            tokens = line.split()
            kw = tokens[0]
            if kw == "sS":
                same_sheet.append(tuple(int(t) for t in tokens[1:]))
            elif kw == "sD":
                diff_sheet.append(tuple(int(t) for t in tokens[1:]))
            elif kw == "hand":
                if len(tokens) == 5:
                    handedness.append(
                        (int(tokens[1]), int(tokens[2]), int(tokens[3]), tokens[4])
                    )
            i += 1

        yield SSPRecord(
            dim=dim,
            skeleton_id=skel,
            third_idx=third,
            sub_first=sub_first,
            sub_second=sub_second,
            sse_types=sse_types,
            matrix=tuple(tuple(row) for row in matrix),
            same_sheet=tuple(same_sheet),
            diff_sheet=tuple(diff_sheet),
            handedness=tuple(handedness),
        )
