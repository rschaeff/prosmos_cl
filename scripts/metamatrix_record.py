#!/usr/bin/env python3
"""Read generateMatrix records: elements, sheets and the full N x N interaction matrix.

Works on a single generateMatrix output file or a whole metamatricesDB (records are
simply concatenated). Import it rather than writing another parser:

    from metamatrix_record import iter_records
    for rec in iter_records("metamatricesDB.clean"):
        rec.name, rec.elements[0].start, rec.matrix[0][2], rec.sheets

or print records from the command line:

    metamatrix_record.py <file_or_db> [NAME ...]     # NAME with or without .ssd

Record layout (generateMatrix/src/control.h, proInteractionMatr)
----------------------------------------------------------------
    header   name %-32s | N %4d | N element blocks of 68 bytes each:
               type(1) chain(1) start %5s "--" end %5s " " length %4d " "
               then 6 x %8.3f: first-point x y z, last-point x y z
    sheets   "sheet <id> <n_strands> <element numbers...>" lines, or "sheet    0"
    matrix   the upper triangle of the N x N matrix, row by row, diagonal ('*')
             included, so N(N+1)/2 characters; an empty line when N == 0

The header is FIXED WIDTH: slice it, never split it or regex it. Adjacent
coordinates abut when a value fills all 8 columns ("  43.6281005.842"), the
second character is the chain id (not always 'A'), and start/end are strings
that can carry insertion codes ("  1B") and minus signs. A record whose fields
have overflowed their widths, or whose matrix does not fit its element count, is
skipped with a note on stderr, as searchmatrix skips it (MalformedRecord is
raised instead with skip_malformed=False).

Matrix codes (readme.original, "Parameters of ProSMoS meta-matrices"):
    c  parallel strands, > 2 H-bonds        t  antiparallel strands, > 2 H-bonds
    u  interact, angle < 85 deg             v  interact, angle >= 95 deg
    N  interact, 85 <= angle < 95 deg       -  no interaction      * diagonal

Element numbers in `sheets` are 1-based, as written; `elements` and `matrix`
are ordinary 0-based Python lists, so sheet member k is elements[k - 1].
"""
from __future__ import annotations

import sys
from dataclasses import dataclass
from typing import Iterator, List, Tuple

NAME_WIDTH = 32
COUNT_WIDTH = 4
FIRST_ELEMENT = NAME_WIDTH + COUNT_WIDTH   # 36
ELEMENT_WIDTH = 68


class MalformedRecord(ValueError):
    pass


@dataclass
class Element:
    type: str                              # 'E' strand or 'H' helix
    chain: str
    start: str                             # residue ids as written: may be '1B', '-3'
    end: str
    length: int
    first: Tuple[float, float, float]      # axis end points
    last: Tuple[float, float, float]


@dataclass
class Sheet:
    id: int
    members: List[int]                     # 1-based element numbers


@dataclass
class Record:
    name: str                              # as written, e.g. 'ecod_1715949.ssd'
    elements: List[Element]
    sheets: List[Sheet]
    matrix: List[List[str]]                # N x N, symmetric, '*' on the diagonal

    @property
    def uid(self) -> str:
        return self.name[:-4] if self.name.endswith(".ssd") else self.name

    def code(self, i: int, j: int) -> str:
        """Interaction code between 1-based elements i and j, as in queries and sheets."""
        return self.matrix[i - 1][j - 1]


def parse_record(header: str, sheet_lines: List[str], matrix_line: str) -> Record:
    header = header.rstrip("\n")
    name = header[:NAME_WIDTH].strip()
    try:
        n = int(header[NAME_WIDTH:FIRST_ELEMENT])
    except ValueError:
        raise MalformedRecord(f"{name}: element count unreadable: {header[NAME_WIDTH:FIRST_ELEMENT]!r}")
    if len(header) != FIRST_ELEMENT + ELEMENT_WIDTH * n:
        raise MalformedRecord(f"{name}: header is {len(header)} bytes, expected "
                              f"{FIRST_ELEMENT + ELEMENT_WIDTH * n} for {n} elements "
                              f"(a field overflowed its width)")
    elements = []
    for k in range(n):
        b = header[FIRST_ELEMENT + ELEMENT_WIDTH * k: FIRST_ELEMENT + ELEMENT_WIDTH * (k + 1)]
        if b[0] not in "EH" or b[7:9] != "--":
            raise MalformedRecord(f"{name}: element {k + 1} misaligned: {b[:20]!r}")
        xyz = [float(b[20 + 8 * c: 28 + 8 * c]) for c in range(6)]
        elements.append(Element(b[0], b[1], b[2:7].strip(), b[9:14].strip(),
                                int(b[15:19]), tuple(xyz[:3]), tuple(xyz[3:])))
    sheets = []
    for line in sheet_lines:
        f = line.split()
        if len(f) > 2:                     # "sheet    0" means no sheets
            sheets.append(Sheet(int(f[1]), [int(x) for x in f[3:]]))
    tri = matrix_line.rstrip("\n")
    if len(tri) != n * (n + 1) // 2:
        raise MalformedRecord(f"{name}: matrix line has {len(tri)} codes, "
                              f"expected {n * (n + 1) // 2} for {n} elements")
    matrix = [[""] * n for _ in range(n)]
    k = 0
    for i in range(n):
        for j in range(i, n):
            matrix[i][j] = matrix[j][i] = tri[k]
            k += 1
    return Record(name, elements, sheets, matrix)


def _is_header(line: str) -> bool:
    return ".ssd" in line[:40] and not line.startswith(("*", "sheet"))


def iter_records(path: str, skip_malformed: bool = True) -> Iterator[Record]:
    """Yield every record in a generateMatrix output file or metamatricesDB.

    Malformed records are reported on stderr and skipped, as searchmatrix skips them;
    pass skip_malformed=False to raise MalformedRecord instead. Even .clean DBs hold
    a few: ecod_db_exp_full/metamatricesDB.clean has 75 of 1,400,111, all headers
    with overflowed fields (residue '2--0', coordinates like 35049.429).
    """
    with open(path, errors="replace") as fh:
        lines = iter(fh)
        line = next(lines, None)
        while line is not None:
            if not _is_header(line):
                msg = f"{path}: expected a record header, got {line[:40]!r}"
                if not skip_malformed:
                    raise MalformedRecord(msg)
                print(f"skipped: {msg}", file=sys.stderr)
                line = next(lines, None)
                continue
            header, sheet_lines = line, []
            line = next(lines, None)
            while line is not None and line.startswith("sheet"):
                sheet_lines.append(line)
                line = next(lines, None)
            matrix_line = "" if line is None else line   # N == 0 writes an empty line
            line = next(lines, None)
            try:
                yield parse_record(header, sheet_lines, matrix_line)
            except MalformedRecord as e:
                if not skip_malformed:
                    raise
                print(f"skipped: {e}", file=sys.stderr)


def format_record(rec: Record) -> str:
    out = [f"{rec.name}  {len(rec.elements)} elements"]
    for i, e in enumerate(rec.elements, 1):
        out.append(f"{i:4d}  {e.type}{e.chain} {e.start:>5}-{e.end:<5} len {e.length}")
    for s in rec.sheets:
        out.append(f"sheet {s.id}: elements {' '.join(map(str, s.members))}")
    if rec.matrix:
        n = len(rec.matrix)
        out.append("      " + " ".join(f"{j:>2}" for j in range(1, n + 1)))
        for i, row in enumerate(rec.matrix, 1):
            out.append(f"{i:4d}  " + " ".join(f"{c:>2}" for c in row))
    return "\n".join(out)


def main(argv: List[str]) -> int:
    if not argv or argv[0] in ("-h", "--help"):
        print(__doc__.split("Record layout")[0].rstrip())
        return 0 if argv else 2
    path, wanted = argv[0], {w if w.endswith(".ssd") else w + ".ssd" for w in argv[1:]}
    shown = 0
    for rec in iter_records(path):
        if wanted and rec.name not in wanted:
            continue
        print(format_record(rec) + "\n")
        shown += 1
    if wanted and shown < len(wanted):
        print(f"found {shown} of {len(wanted)} requested records", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
