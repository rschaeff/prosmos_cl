#!/usr/bin/env python3
"""Post-ProSMoS acceptance rule: turn searchmatrix hits into accepted occurrences.

searchmatrix returns *candidates*: elements whose pairwise codes match the query and that
share one PALSSE sheet (the query's `sheetS` line). The four-strand census accepts a
candidate only if it also passes four clauses, run here as four sequential stages:

    hits      searchmatrix output  ->  one row per candidate instance
    sieve     1. disjointness      the four residue ranges share no residue
              2. axial overlap     on every interior strand, both flanks overlap it by
                                   more than 3 A along its PALSSE axis
    pairing   3. local pairing     on every interior strand, >= 2 consecutive residues whose
                                   nearest flank Ca atoms lie within 6.0 A of Ca(i) and on
                                   opposite sides of the plane Ca(i-1), Ca(i), Ca(i+1)
    bridge    4. DSSP bridges      every pair the query marks c/t carries >= 1 DSSP bridge
    census    count accepted (query, record) pairs at three sheet-size conventions

Every stage reads the previous stage's table and writes every row back out. A row that
fails a clause keeps its place with `status` set to `reject:<clause>`; a row that cannot be
evaluated gets `fail:<reason>`. Nothing is dropped silently, and `census` prints the full
status breakdown. Only rows whose status is `ok` are processed by later stages.

Rejects and failures are different things. `reject:*` is the rule working. `fail:*` means
the input could not be judged -- no structure file, DSSP could not assign the structure,
a record name searchmatrix truncated beyond recovery -- and such a row is NOT accepted.

Usage (see RUNBOOK_GREY.md for a worked example):
    acceptance.py hits    --hits DIR_OR_TSV [...] --queries QDIR --db DB.clean -o t1.tsv
    acceptance.py sieve   -i t1.tsv --db DB.clean --queries QDIR -o t2.tsv
    acceptance.py pairing -i t2.tsv --structures TEMPLATE --queries QDIR -o t3.tsv
    acceptance.py bridge  -i t3.tsv --structures TEMPLATE --queries QDIR --dssp BIN -o t4.tsv
    acceptance.py census  -i t4.tsv --queries QDIR -o census.json

`--shard K/N` on pairing and bridge keeps only records with index K of every N (0-based),
for SLURM arrays; the next stage takes all shard outputs as its -i list.

Exit status: 0 on success, 1 on an input problem (missing file, no hits at all, a malformed
query), 2 on bad arguments. Row-level failures do not change the exit status; read the
breakdown.

Provenance: this is a portable re-implementation of the rule behind census_plane.json in
the ProSMoS 2026 paper (ruczinski/bridge_prep.py with the sign clause off, nvg_worker.py,
bridge_worker.py + bridge_rules.py, census_plane.py). It reproduces those per-instance
verdicts except where the originals were silently lossy; the differences are listed in
scripts/acceptance/README.md.
"""
from __future__ import annotations

import argparse
import csv
import glob
import json
import math
import os
import re
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from metamatrix_record import MalformedRecord, _is_header, parse_record  # noqa: E402

OV_MIN = 3.0          # axial overlap floor, A
PAIR_D = 6.0          # local pairing distance, A
PAIR_RUN = 2          # consecutive qualifying residues
BRIDGE_MIN = 1        # DSSP bridges per required pair
CONVENTIONS = (("complete", 4), ("core", 6), ("raw", None))
VACUOUS = 999.0

COLS = ["query", "record", "sheet_size", "positions", "ranges", "status", "detail",
        "overlap_min", "pairing_dstar", "min_pair", "min_sand"]


class InputError(Exception):
    """An input problem the user has to fix: reported on stderr, exit 1."""


# ----------------------------------------------------------------------------- queries

class Query:
    def __init__(self, name, types, codes):
        self.name, self.types, self.codes = name, types, codes     # codes[(i, j)], 1-based
        n = len(types)
        adj = {p: [] for p in range(1, n + 1)}
        self.required = []
        for (i, j), c in sorted(codes.items()):
            if c in "ct":
                self.required.append((i, j))
                adj[i].append(j); adj[j].append(i)
        self.interiors = [(p, tuple(adj[p])) for p in range(1, n + 1) if len(adj[p]) == 2]


def read_query(path):
    lines = [ln.rstrip("\n") for ln in open(path)]
    try:
        ids = lines[0].split()
        types = lines[1].split()
        n = len(ids)
        if n < 2 or len(types) != n:
            raise ValueError("element and type lines disagree")
        codes = {}
        for i in range(n):
            row = lines[2 + i].split()
            if not row or row[0] != "*":
                raise ValueError(f"matrix row {i + 1} does not start with '*'")
            for k, c in enumerate(row[1:]):
                codes[(i + 1, i + 2 + k)] = c
        if len(codes) != n * (n - 1) // 2:
            raise ValueError("matrix is not a full upper triangle")
    except IndexError:
        raise InputError(f"{path}: malformed query (file ends inside the element/matrix lines)")
    except ValueError as e:
        raise InputError(f"{path}: malformed query ({e})")
    return Query(Path(path).stem, types, codes)


def load_queries(qdir):
    files = sorted(glob.glob(os.path.join(qdir, "*.query")))
    if not files:
        raise InputError(f"no *.query files in {qdir}")
    return {q.name: q for q in map(read_query, files)}


# ----------------------------------------------------------------------------- table io

def read_table(paths):
    rows = []
    for p in paths:
        if not os.path.exists(p):
            raise InputError(f"input table {p} does not exist")
        with open(p, newline="") as fh:
            rd = csv.DictReader(fh, delimiter="\t")
            if rd.fieldnames != COLS:
                raise InputError(f"{p}: not an acceptance table (header {rd.fieldnames})")
            rows.extend(rd)
    return rows


def write_table(rows, path):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "") for c in COLS})


def report(stage, rows, path):
    c = Counter(r["status"] for r in rows)
    parts = ", ".join(f"{k} {v:,}" for k, v in sorted(c.items(), key=lambda kv: -kv[1]))
    print(f"{stage}: {len(rows):,} rows -> {path}  ({parts})", file=sys.stderr)


def shard_filter(rows, shard):
    if not shard:
        return rows
    k, n = shard
    recs = sorted({r["record"] for r in rows})
    mine = set(recs[k::n])
    return [r for r in rows if r["record"] in mine]


def parse_shard(s):
    m = re.fullmatch(r"(\d+)/(\d+)", s or "")
    if not m or int(m.group(1)) >= int(m.group(2)):
        raise argparse.ArgumentTypeError("--shard must be K/N with 0 <= K < N")
    return int(m.group(1)), int(m.group(2))


# ----------------------------------------------------------------------------- metamatrix

def scan_db(db, keep):
    """Parse the records `keep(name)` selects. Returns (records by name, malformed names)."""
    if not os.path.exists(db):
        raise InputError(f"metamatrix DB {db} does not exist")
    recs, bad = {}, set()
    with open(db, errors="replace") as fh:
        lines = iter(fh)
        line = next(lines, None)
        while line is not None:
            if not _is_header(line):
                line = next(lines, None)
                continue
            header, sheets = line, []
            line = next(lines, None)
            while line is not None and line.startswith("sheet"):
                sheets.append(line)
                line = next(lines, None)
            matrix = "" if line is None else line
            line = next(lines, None)
            name = header[:32].strip()
            uid = name[:-4] if name.endswith(".ssd") else name
            if not keep(uid):
                continue
            try:
                recs[uid] = parse_record(header, sheets, matrix)
            except MalformedRecord:
                bad.add(uid)
    return recs, bad


RESID = re.compile(r"^(-?\d+)([A-Za-z]?)$")


def resid(s):
    """'12' -> (12, ''), '12A' -> (12, 'A'), '-3' -> (-3, ''). None if unparseable."""
    m = RESID.match(s.strip())
    return (int(m.group(1)), m.group(2)) if m else None


# ----------------------------------------------------------------------------- stage: hits

# Ranges can carry insertion codes ("82A") and minus signs; archive_hits.py's `\d+` dropped
# such segments, which is why an archived row can hold fewer segments than the query has.
SEG = re.compile(r"segment-Type:\s*(\S+)\s+Position:\s*(\d+)\s+"
                 r"Range:\s*(-?\d+[A-Za-z]?)\s*--\s*(-?\d+[A-Za-z]?)")


def parse_hit_file(path):
    """searchmatrix hit file -> (record name as written, [[(type, position, lo, hi)], ...])."""
    blocks, cur, name = [], None, None
    with open(path, errors="replace") as fh:
        for ln in fh:
            if name is None and ln.strip():
                name = ln.strip()
            if ln.startswith("MOTIF:"):
                if cur:
                    blocks.append(cur)
                cur = []
            elif ln.startswith("END"):
                if cur:
                    blocks.append(cur)
                cur = None
            elif cur is not None:
                m = SEG.search(ln)
                if m:
                    cur.append((m.group(1), int(m.group(2)), resid(m.group(3))[0], resid(m.group(4))[0]))
    if cur:
        blocks.append(cur)
    return name, blocks


def iter_hits(paths, single_query):
    """Yield (query, record-as-written, segments) from searchmatrix dirs or archived TSVs."""
    found = False
    for p in paths:
        if os.path.isdir(p):
            files = glob.glob(os.path.join(p, "pdb*.txt"))
            if files:                                    # one query's output directory
                if not single_query:
                    raise InputError(f"{p} holds hit files directly; name its query with --query")
                layout = [(single_query, files)]
            else:                                        # manifest layout: <p>/<query>/pdb*.txt
                layout = [(Path(d).name, glob.glob(os.path.join(d, "pdb*.txt")))
                          for d in sorted(glob.glob(os.path.join(p, "*"))) if os.path.isdir(d)]
            for q, fs in layout:
                for f in sorted(fs):
                    name, blocks = parse_hit_file(f)
                    rec = name[3:] if name and name.startswith("pdb") else Path(f).stem[3:]
                    for b in blocks:
                        found = True
                        yield q, rec, b
        elif os.path.isfile(p):                          # archive_hits.py TSV, optionally zstd
            if p.endswith(".zst"):
                proc = subprocess.Popen(["zstd", "-dc", p], stdout=subprocess.PIPE, text=True)
                fh = proc.stdout
            else:
                proc, fh = None, open(p)
            for ln in fh:
                t = ln.rstrip("\n").split("\t")
                if len(t) < 4 or t[0] == "query":
                    continue
                segs = []
                for s in t[3].split(","):
                    typ, pos, rng = s.split(":")
                    m = re.fullmatch(r"(-?\d+[A-Za-z]?)-(-?\d+[A-Za-z]?)", rng)
                    segs.append((typ, int(pos), resid(m.group(1))[0] if m else None,
                                 resid(m.group(2))[0] if m else None))
                found = True
                yield t[0], t[1], segs
            if proc:
                proc.wait()
            else:
                fh.close()
        else:
            raise InputError(f"hits path {p} does not exist")
    if not found:
        raise InputError("no hits found in any --hits path (an empty searchmatrix output is "
                         "almost always an error upstream, not a real zero)")


def smallest_sheet(rec, positions):
    ps = set(positions)
    sizes = [len(s.members) for s in rec.sheets if ps <= set(s.members)]
    return min(sizes) if sizes else None


def cmd_hits(a):
    queries = load_queries(a.queries)
    raw = list(iter_hits(a.hits, a.query))
    if a.records:
        wanted = {ln.strip() for ln in open(a.records) if ln.strip()}
        raw = [h for h in raw if h[1] in wanted]
    names = {h[1] for h in raw}
    lengths = sorted({len(n) for n in names})
    candidates = defaultdict(list)

    def keep(uid):
        if uid in names:
            return True
        hit = False
        for L in lengths:
            if L < len(uid) and uid[:L] in names:
                candidates[uid[:L]].append(uid)
                hit = True
        return hit

    recs, bad = scan_db(a.db, keep)
    rows = []
    for q, name, segs in raw:
        row = {"query": q, "record": name, "positions": ",".join(str(s[1]) for s in segs),
               "status": "ok", "detail": ""}
        rows.append(row)
        query = queries.get(q)
        if query is None:
            row.update(status="fail:unknown_query", detail=f"no {q}.query in {a.queries}")
            continue
        if len(segs) != len(query.types):
            row.update(status="fail:segment_count",
                       detail=f"{len(segs)} segments for a {len(query.types)}-element query; an "
                              "archived hit table drops segments whose range has an insertion "
                              "code -- rerun the search, or read the hit tree, to recover it")
            continue
        # resolve the record: exact name, else a unique truncation match on element ranges
        pool = [name] if name in recs or name in bad else candidates.get(name, [])
        if not pool:
            row.update(status="fail:record_not_in_db")
            continue
        fits = []
        for uid in pool:
            if uid in bad:
                continue
            r = recs[uid]
            ok = all(1 <= p <= len(r.elements) for _, p, _, _ in segs)
            if ok:
                for _, p, lo, hi in segs:
                    e = r.elements[p - 1]
                    s, t = resid(e.start), resid(e.end)
                    if lo is not None and (s is None or t is None or (s[0], t[0]) != (lo, hi)):
                        ok = False
                        break
            if ok:
                fits.append(uid)
        if not fits:
            if all(u in bad for u in pool):
                row.update(status="fail:malformed_record",
                           detail="metamatrix record has overflowed fields")
            else:
                row.update(status="fail:range_mismatch",
                           detail=f"hit ranges match no element set in {len(pool)} candidate(s)")
            continue
        if len(fits) > 1:
            row.update(status="fail:ambiguous_record",
                       detail=f"truncated name matches {len(fits)} records: {','.join(fits[:4])}")
            continue
        uid = fits[0]
        r = recs[uid]
        pos = [s[1] for s in segs]
        row["record"] = uid
        if uid != name:
            row["detail"] = f"resolved from truncated name {name}"
        els = [r.elements[p - 1] for p in pos]
        row["ranges"] = ",".join(f"{e.chain}:{e.start}-{e.end}" for e in els)
        wrong = [(i, j) for (i, j) in query.required
                 if r.code(pos[i - 1], pos[j - 1]) != query.codes[(i, j)]]
        if wrong:
            row.update(status="fail:slot_mismatch",
                       detail=f"record codes disagree with query on pairs {wrong}")
            continue
        sz = smallest_sheet(r, pos)
        if sz is None:
            row.update(status="fail:no_sheet", detail="no PALSSE sheet holds all four elements")
            continue
        row["sheet_size"] = sz
    write_table(rows, a.out)
    report("hits", rows, a.out)


# ----------------------------------------------------------------------------- stage: sieve

def axial_overlap(I, A, B):
    """Overlap (A) of flanks A and B along interior element I's axis; None if I is degenerate."""
    e1 = I.first
    u = [I.last[k] - e1[k] for k in range(3)]
    nu = math.sqrt(sum(x * x for x in u))
    if nu < 1e-6:
        return None
    u = [x / nu for x in u]
    proj = lambda p: sum((p[k] - e1[k]) * u[k] for k in range(3))
    pa = [proj(A.first), proj(A.last)]
    pb = [proj(B.first), proj(B.last)]
    return min(max(pa), max(pb)) - max(min(pa), min(pb))


def cmd_sieve(a):
    queries = load_queries(a.queries)
    rows = read_table(a.inp)
    need = {r["record"] for r in rows if r["status"] == "ok"}
    recs, bad = scan_db(a.db, lambda uid: uid in need)
    for row in rows:
        if row["status"] != "ok":
            continue
        uid = row["record"]
        if uid not in recs:
            row.update(status="fail:record_not_in_db" if uid not in bad else "fail:malformed_record")
            continue
        r, q = recs[uid], queries[row["query"]]
        pos = [int(x) for x in row["positions"].split(",")]
        els = [r.elements[p - 1] for p in pos]
        # 1. disjointness
        spans = []
        for e in els:
            s, t = resid(e.start), resid(e.end)
            if s is None or t is None:
                spans = None
                break
            spans.append((s, t))
        if spans is None:
            row.update(status="fail:unparsed_residue", detail=row["ranges"])
            continue
        spans.sort()
        if not all(spans[i][1] < spans[i + 1][0] for i in range(len(spans) - 1)):
            row["status"] = "reject:disjoint"
            continue
        # 2. axial overlap on every interior strand
        ovs, degenerate = [], False
        for ip, (na, nb) in q.interiors:
            ov = axial_overlap(els[ip - 1], els[na - 1], els[nb - 1])
            if ov is None:
                degenerate = True
                break
            ovs.append(ov)
        if degenerate:
            row.update(status="fail:degenerate_axis", detail="interior element has zero length")
            continue
        row["overlap_min"] = f"{min(ovs):.2f}" if ovs else ""
        if ovs and min(ovs) <= OV_MIN:
            row["status"] = "reject:overlap"
    write_table(rows, a.out)
    report("sieve", rows, a.out)


# ----------------------------------------------------------------------------- structures

def check_structures(template):
    if not os.path.isdir(template) and "{rec}" not in template:
        raise InputError(f"--structures {template} is neither a directory nor a template with {{rec}}")


def none_found(by_rec, template):
    """Every record missing its structure is a path problem, not a data problem."""
    if by_rec and not any(os.path.exists(structure_path(template, u)) for u in by_rec):
        eg = structure_path(template, next(iter(sorted(by_rec))))
        raise InputError(f"no structure found for any of {len(by_rec):,} records "
                         f"(first expected at {eg}); check --structures")


def structure_path(template, uid):
    """Fill a path template. {rec} is the record name; {bucket} is int(rec)//100 as %05d,
    for ECOD's /data/ecod/af2_pdb_domain_data/{bucket}/{rec}/{rec}.pdb layout."""
    fields = {"rec": uid}
    if "{bucket}" in template:
        fields["bucket"] = f"{int(uid) // 100:05d}" if uid.isdigit() else ""
    if os.path.isdir(template):
        return os.path.join(template, f"{uid}.pdb")
    return template.format(**fields)


def read_ca(path):
    """(chain, resnum, icode) -> Ca xyz, ATOM records only, altloc ' ' or 'A'.
    Also returns the atom names seen, for diagnosing structures DSSP cannot assign."""
    ca, names, left = {}, set(), False
    with open(path, errors="replace") as fh:
        for ln in fh:
            if not ln.startswith("ATOM"):
                continue
            nm = ln[12:16]
            names.add(nm.strip())
            if nm[0] != " " and len(nm.strip()) < 4 and not nm[0].isdigit():
                left = True
            if nm.strip() != "CA" or ln[16] not in (" ", "A"):
                continue
            try:
                key = (ln[21], int(ln[22:26]), ln[26].strip())
                ca[key] = (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
            except ValueError:
                continue
    return ca, names, left


def element_keys(ca, chain, span):
    """Ca keys of one element, in residue order. Falls back to all chains when the
    element's chain letter is absent from the file (single-chain domain files often
    relabel)."""
    s, t = span
    chains = {k[0] for k in ca}
    use = chain if chain in chains else None
    keys = [k for k in ca if (use is None or k[0] == use) and s <= (k[1], k[2]) <= t]
    return sorted(keys, key=lambda k: (k[0], k[1], k[2]))


def consecutive(k1, k2):
    return k1[0] == k2[0] and (k2[1] - k1[1] == 1 or (k2[1] == k1[1] and k2[2] != k1[2]))


def parse_ranges(ranges):
    out = []
    for g in ranges.split(","):
        chain, rng = g.split(":", 1)
        m = re.fullmatch(r"(-?\d+[A-Za-z]?)-(-?\d+[A-Za-z]?)", rng)
        out.append((chain, (resid(m.group(1)), resid(m.group(2)))))
    return out


# ----------------------------------------------------------------------------- stage: pairing

def sub(p, q): return (p[0] - q[0], p[1] - q[1], p[2] - q[2])
def dot(p, q): return p[0] * q[0] + p[1] * q[1] + p[2] * q[2]
def cross(p, q): return (p[1] * q[2] - p[2] * q[1], p[2] * q[0] - p[0] * q[2], p[0] * q[1] - p[1] * q[0])
def dist(p, q): d = sub(p, q); return math.sqrt(dot(d, d))


def dstar(ca, B, A, C, run):
    """Smallest cutoff at which interior strand B has `run` consecutive qualifying residues.

    Per non-terminal residue of B: the plane through Ca(i-1), Ca(i), Ca(i+1); each flank's
    nearest Ca must lie on opposite sides of it; the residue's value is the larger of the
    two nearest distances (VACUOUS if the sides test fails or a neighbour Ca is missing).
    """
    if len(B) < 3 or not A or not C:
        return VACUOUS
    Ax, Cx = [ca[k] for k in A], [ca[k] for k in C]
    per = []
    for i in range(1, len(B) - 1):
        k = B[i]
        p = ca[k]
        prev, nxt = B[i - 1], B[i + 1]
        if not (consecutive(prev, k) and consecutive(k, nxt)):
            per.append((k, VACUOUS)); continue
        n = cross(sub(p, ca[prev]), sub(ca[nxt], ca[prev]))
        if dot(n, n) < 1e-12:
            per.append((k, VACUOUS)); continue
        da, qa = min(((dist(p, q), q) for q in Ax), key=lambda t: t[0])
        dc, qc = min(((dist(p, q), q) for q in Cx), key=lambda t: t[0])
        if dot(sub(qa, p), n) * dot(sub(qc, p), n) >= 0:
            per.append((k, VACUOUS)); continue
        per.append((k, max(da, dc)))
    best = VACUOUS
    for i in range(len(per) - run + 1):
        w = per[i:i + run]
        if all(consecutive(w[j][0], w[j + 1][0]) for j in range(run - 1)):
            best = min(best, max(x[1] for x in w))
    return best


def cmd_pairing(a):
    queries = load_queries(a.queries)
    check_structures(a.structures)
    rows = shard_filter(read_table(a.inp), a.shard)
    by_rec = defaultdict(list)
    for r in rows:
        if r["status"] == "ok":
            by_rec[r["record"]].append(r)
    none_found(by_rec, a.structures)
    for uid, rs in sorted(by_rec.items()):
        path = structure_path(a.structures, uid)
        if not os.path.exists(path):
            for r in rs:
                r.update(status="fail:no_structure", detail=path)
            continue
        ca, _, _ = read_ca(path)
        if not ca:
            for r in rs:
                r.update(status="fail:no_ca_atoms", detail=path)
            continue
        for r in rs:
            q = queries[r["query"]]
            els = parse_ranges(r["ranges"])
            keys = [element_keys(ca, ch, sp) for ch, sp in els]
            ds = [dstar(ca, keys[ip - 1], keys[na - 1], keys[nb - 1], a.run)
                  for ip, (na, nb) in q.interiors]
            r["pairing_dstar"] = ",".join(f"{d:.2f}" for d in ds)
            if ds and max(ds) > a.distance:
                r["status"] = "reject:pairing"
    write_table(rows, a.out)
    report(f"pairing (D {a.distance} A, run {a.run})", rows, a.out)


# ----------------------------------------------------------------------------- stage: bridge

def run_dssp(dssp, pdb):
    """Classic DSSP text output (dsspcmbi or mkdssp --output-format dssp) -> residue rows."""
    with tempfile.TemporaryDirectory() as td:
        out = os.path.join(td, "out.dssp")
        cmd = [dssp, pdb, out]
        if "mkdssp" in os.path.basename(dssp):
            cmd = [dssp, "--output-format", "dssp", pdb, out]
        subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        if not os.path.exists(out):
            return None
        lines = Path(out).read_text(errors="replace").splitlines()
    st = [i for i, l in enumerate(lines) if l.startswith("  #  RESIDUE")]
    if not st:
        return None
    rows = []
    for ln in lines[st[0] + 1:]:
        if len(ln) < 38 or ln[13] == "!":
            continue
        try:
            rows.append((int(ln[0:5]), (ln[11], int(ln[5:10]), ln[10].strip()),
                         int(ln[25:29]), int(ln[29:33])))
        except ValueError:
            continue
    return rows or None


def dssp_failure_reason(path):
    _, names, left = read_ca(path)
    if names <= {"CA"}:
        return "fail:dssp_ca_only"
    if left:
        return "fail:dssp_atom_names"      # atom names left-justified in columns 13-16
    if "O" not in names:
        return "fail:dssp_no_carbonyl"
    return "fail:dssp_no_output"


def cmd_bridge(a):
    queries = load_queries(a.queries)
    if not (os.path.isfile(a.dssp) and os.access(a.dssp, os.X_OK)):
        raise InputError(f"DSSP binary {a.dssp} is not an executable file")
    check_structures(a.structures)
    rows = shard_filter(read_table(a.inp), a.shard)
    by_rec = defaultdict(list)
    for r in rows:
        if r["status"] == "ok":
            by_rec[r["record"]].append(r)
    none_found(by_rec, a.structures)
    for uid, rs in sorted(by_rec.items()):
        path = structure_path(a.structures, uid)
        if not os.path.exists(path):
            for r in rs:
                r.update(status="fail:no_structure", detail=path)
            continue
        drows = run_dssp(a.dssp, path)
        if drows is None:
            why = dssp_failure_reason(path)
            for r in rs:
                r.update(status=why, detail="DSSP assigned no residues")
            continue
        n2k = {n: k for n, k, _, _ in drows}
        partners = defaultdict(set)
        for _, k, b1, b2 in drows:
            for b in (b1, b2):
                if b and b in n2k:
                    partners[k].add(n2k[b])
        chains = {k[0] for _, k, _, _ in drows}
        for r in rs:
            q = queries[r["query"]]
            spans = parse_ranges(r["ranges"])

            def inel(k, slot):
                ch, (s, t) = spans[slot - 1]
                return s <= (k[1], k[2]) <= t and (k[0] == ch or ch not in chains)

            res = {slot: [k for _, k, _, _ in drows if inel(k, slot)]
                   for slot in range(1, len(spans) + 1)}
            min_pair = int(VACUOUS)
            for x, y in q.required:
                n = sum(1 for k in res[x] for p in partners.get(k, ()) if inel(p, y))
                min_pair = min(min_pair, n)
            min_sand = int(VACUOUS)
            for ip, (na, nb) in q.interiors:
                n = sum(1 for k in res[ip]
                        if any(inel(p, na) for p in partners.get(k, ()))
                        and any(inel(p, nb) for p in partners.get(k, ())))
                min_sand = min(min_sand, n)
            r["min_pair"], r["min_sand"] = str(min_pair), str(min_sand)
            if min_pair < a.min_bridges:
                r["status"] = "reject:bridge"
    write_table(rows, a.out)
    report(f"bridge (>= {a.min_bridges} per pair)", rows, a.out)


# ----------------------------------------------------------------------------- stage: census

def cmd_census(a):
    queries = load_queries(a.queries)
    rows = read_table(a.inp)
    best = defaultdict(dict)                   # query -> record -> smallest accepted sheet
    for r in rows:
        if r["status"] != "ok":
            continue
        if not r["min_pair"]:
            raise InputError("census input has ok rows that never went through the bridge "
                             "stage; run hits -> sieve -> pairing -> bridge first")
        sz = int(r["sheet_size"])
        d = best[r["query"]]
        if r["record"] not in d or sz < d[r["record"]]:
            d[r["record"]] = sz
    per = []
    for q in sorted(queries):
        d = best.get(q, {})
        e = {name: sum(1 for s in d.values() if cap is None or s <= cap) for name, cap in CONVENTIONS}
        per.append({"query": q, **e})
    realized = {name: sum(1 for p in per if p[name] > 0) for name, _ in CONVENTIONS}
    status = Counter(r["status"] for r in rows)
    fails = {k: v for k, v in status.items() if k.startswith("fail:")}
    out = {"queries": len(queries), "realized": realized,
           "status": dict(sorted(status.items())), "per_query": per}
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    json.dump(out, open(a.out, "w"), indent=1)
    print(f"census: {len(rows):,} rows; realized " +
          ", ".join(f"{k} {v}/{len(queries)}" for k, v in realized.items()) + f" -> {a.out}")
    for k, v in sorted(status.items(), key=lambda kv: -kv[1]):
        print(f"  {k:<28}{v:>12,}")
    if fails:
        print(f"  {sum(fails.values()):,} rows could not be judged and were NOT accepted "
              f"(fail:*); see the detail column", file=sys.stderr)


# ----------------------------------------------------------------------------- cli

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest="cmd", required=True)

    def common(p, needs_in=True):
        p.add_argument("--queries", required=True, help="directory of *.query files")
        if needs_in:
            p.add_argument("-i", "--inp", nargs="+", required=True, help="input table(s)")
        p.add_argument("-o", "--out", required=True, help="output path")

    p = sp.add_parser("hits", help="searchmatrix output -> candidate table")
    common(p, needs_in=False)
    p.add_argument("--hits", nargs="+", required=True,
                   help="searchmatrix output dir(s) or archive_hits.py TSV(.zst) file(s)")
    p.add_argument("--db", required=True, help="the metamatricesDB(.clean) that was searched")
    p.add_argument("--query", help="query name, when --hits is one query's output directory")
    p.add_argument("--records", help="file of record names to keep (testing, subsetting)")
    p.set_defaults(fn=cmd_hits)

    p = sp.add_parser("sieve", help="clauses 1-2: disjointness and axial overlap")
    common(p)
    p.add_argument("--db", required=True)
    p.set_defaults(fn=cmd_sieve)

    for name, fn, help_ in (("pairing", cmd_pairing, "clause 3: local pairing"),
                            ("bridge", cmd_bridge, "clause 4: DSSP bridge per required pair")):
        p = sp.add_parser(name, help=help_)
        common(p)
        p.add_argument("--structures", required=True,
                       help="directory holding <record>.pdb, or a path template with {rec} "
                            "(and {bucket} for the ECOD store)")
        p.add_argument("--shard", type=parse_shard, help="K/N: process records K, K+N, ...")
        if name == "bridge":
            p.add_argument("--dssp", required=True, help="dsspcmbi or mkdssp binary")
            p.add_argument("--min-bridges", type=int, default=BRIDGE_MIN)
        else:
            p.add_argument("--distance", type=float, default=PAIR_D)
            p.add_argument("--run", type=int, default=PAIR_RUN)
        p.set_defaults(fn=fn)

    p = sp.add_parser("census", help="count accepted occurrences")
    common(p)
    p.set_defaults(fn=cmd_census)

    a = ap.parse_args(argv)
    try:
        a.fn(a)
    except InputError as e:
        print(f"acceptance.py {a.cmd}: {e}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
