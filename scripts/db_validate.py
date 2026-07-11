#!/usr/bin/env python3
"""Validate (and optionally clean) a ProSMoS metamatricesDB.

Why this exists
---------------
searchmatrix's DB reader (searchControl.h:277-324) parses records by an
even/odd line counter: non-`s` lines alternate HEADER (must contain `.ssd`) and
MATRIX (starts with `*`); `s`-prefixed sheet lines are skipped. It assumes every
record is exactly `[header] [0+ sheet lines] [one matrix line]`. A record that
violates this — most commonly an **orphan `sheet`+`matrix` block with no header**
(a generateMatrix glitch) — desyncs the alternation, and the engine's recovery
is incomplete: it silently mis-parses and drops hits across a large downstream
region. On the 4.92M-record AFDB DB this undercounts S5 hits ~40x (see
enum/docs/searchmatrix_fulldb_undercount.md).

This tool simulates that same parser, reports malformed records, and (with
--clean) re-emits ONLY well-formed records so searchmatrix parses the whole DB
in sync.

A well-formed record:
  - exactly one header line (contains `.ssd` within the first 40 chars),
  - then zero or more sheet lines (line[0] == 's'),
  - then exactly one matrix line (line[0] == '*', no `.ssd`),
  - and nothing else before the next header.

Usage:
  db_validate.py <db>                 # report only
  db_validate.py <db> --clean <out>   # write cleaned DB, dropping malformed records
  db_validate.py <db> --show 10       # print context for the first N malformed records
"""
from __future__ import annotations
import argparse
import sys


# searchmatrix parses the header at FIXED byte offsets (searchControl.h:1208-1310):
# name %-32s | count %4d @32 | then per SSE, 68 bytes each starting @36:
#   type,chain (2) | begin %5s (5) | "--" (2) | end %5s (5) | " " | len %4d (4)
#   | " " | 6x coord %8.3f (48).
# So for SSE k the "--" separator must sit at offset 43 + 68*k. If any field
# overflows its width (name >32, residue >5 digits, |coord| >= 10000 or <= -1000,
# count/len >9999) every following field shifts and the "--" moves. This check
# detects that misalignment directly. On the current AFDB DB, zero records fail
# it (max name 24, residues <=4 digits, no coord overflow) -- it's a forward
# guard against a future generateMatrix run that overflows a column.
SSE_STRIDE = 68
FIRST_SSE = 36


def header_field_overflow(line: str) -> bool:
    """True if `line` (a header) is misaligned at searchmatrix's fixed offsets
    (i.e. some field overflowed its fixed width)."""
    name_ssd = line[:line.find(".ssd") + 4] if ".ssd" in line else ""
    if len(name_ssd) > 32:
        return True
    try:
        count = int(line[32:36])
    except ValueError:
        return True
    for k in range(count):
        sep = FIRST_SSE + SSE_STRIDE * k + 7
        if line[sep:sep + 2] != "--":
            return True
    return False


def line_kind(line: str) -> str:
    if not line:
        return "blank"
    if line[0] == "s":
        return "sheet"
    if ".ssd" in line[:40]:
        return "header"
    if line[0] == "*":
        return "matrix"
    return "other"


def validate(path: str, clean_out: str | None = None, show: int = 0):
    total = wellformed = 0
    malformed = {"orphan_block": 0, "header_no_matrix": 0, "extra_matrix": 0, "stray": 0}
    field_overflow = 0
    shown = 0
    out = open(clean_out, "w") if clean_out else None

    header = None          # current record's header line
    sheets: list[str] = []
    have_header = False

    def flush_ok(matrix_line: str):
        nonlocal wellformed
        wellformed += 1
        if out:
            out.write(header)
            for s in sheets:
                out.write(s)
            out.write(matrix_line)

    def report(kind: str, ctx: list[str]):
        nonlocal shown
        malformed[kind] += 1
        if show and shown < show:
            shown += 1
            sys.stderr.write(f"--- malformed [{kind}] ---\n")
            for L in ctx[-6:]:
                sys.stderr.write(f"   {line_kind(L):8s} {L[:78].rstrip()}\n")

    recent: list[str] = []
    with open(path) as f:
        for line in f:
            recent.append(line)
            if len(recent) > 8:
                recent.pop(0)
            k = line_kind(line)
            if k == "blank":
                continue
            if k == "header":
                if have_header:
                    # previous header never reached a matrix
                    report("header_no_matrix", recent)
                total += 1
                if header_field_overflow(line):
                    field_overflow += 1
                    if show and shown < show:
                        shown += 1
                        sys.stderr.write(f"--- field-overflow (misaligned at fixed offsets) ---\n   {line[:78].rstrip()}\n")
                header, sheets, have_header = line, [], True
            elif k == "sheet":
                if have_header:
                    sheets.append(line)
                # else: sheet with no open header -> part of an orphan block (ignored)
            elif k == "matrix":
                if have_header:
                    flush_ok(line)
                    have_header = False
                    header, sheets = None, []
                else:
                    report("orphan_block", recent)   # matrix where a header was expected
            else:  # "other"
                report("stray", recent)
        if have_header:
            report("header_no_matrix", recent)
    if out:
        out.close()

    bad = sum(malformed.values())
    print(f"records (headers) seen : {total}")
    print(f"well-formed            : {wellformed}")
    print(f"malformed              : {bad}")
    for k, v in malformed.items():
        if v:
            print(f"    {k:18s}: {v}")
    print(f"field-overflow (fixed-column misalignment): {field_overflow}")
    if clean_out:
        print(f"cleaned DB written     : {clean_out}  ({wellformed} records)")
    return bad


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("db")
    ap.add_argument("--clean", metavar="OUT", help="write cleaned DB here")
    ap.add_argument("--show", type=int, default=0, help="print context for first N malformed")
    a = ap.parse_args()
    bad = validate(a.db, a.clean, a.show)
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
