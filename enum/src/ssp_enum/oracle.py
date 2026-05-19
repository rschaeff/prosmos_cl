"""Parser for Chalam's CG-2012 IA.txt format.

Used solely as a validation oracle. Each block in IA.txt is one SSP record:

    [5-269-0]      <- outer header: dim-skeletonId-thirdIdx
    [0-1911]       <- inner header: sub-block id (purpose unclear; serial id)
     1 2 3 4 5     <- SSE index header
     E E E E E     <- SSE type per index (E=strand, H=helix)
     * - c - c     <- upper-triangular interaction matrix
       * t c -
         * - -
           * -
             *
    sS 1 3 5 2 4   <- sheetS: which SSEs share a sheet (order: spatial, R-to-L)
    sD 1 2         <- sheetD: SSE pairs explicitly in different sheets
    hand 1 2 5 L   <- handedness of SSE triple (L = left, R = right)

This module reads that format into structured Python records so we can:
  - count SSPs per dimension (validation against CG-2012 Info.txt totals)
  - look up specific SSPs by id (cross-check against our fresh enumeration)
  - never accept IA.txt as input to anything downstream (use ProSMoS queries)

This file is the *only* place in the package that consumes IA.txt.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator


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
    """Yield one SSPRecord per entry in an IA.txt file. Not yet implemented."""
    raise NotImplementedError("parse() pending; first commit is scaffold-only")
