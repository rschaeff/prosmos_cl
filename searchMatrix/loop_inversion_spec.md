# searchmatrix loop-inversion — spec

## Problem

Today searchmatrix takes **one query** and streams the **whole DB**, parsing every
record, then matching that single query. The proteome sweep runs it **once per
query** (6,336 S5 queries as 6,336 SLURM tasks). So the 2.6 GB / 4.92M-record DB
is **parsed 6,336 times** — ~16 TB of redundant reads and 6,336× redundant record
parsing (header field-splitting, per-element `new`/`delete`, matrix build). The
*matching* work is irreducible (record × query), but the parsing and I/O are not.

**Goal:** parse the DB **once**, test **all queries** against each record. Parsing
and I/O drop ~6,336×; matching is unchanged. Expected ~2× total-compute reduction
(parsing is a large fraction of per-record cost after the `mkdir` fix) plus the
elimination of 16 TB of redundant DB reads.

## Current shape (searchControl.h `oneprocess`, searchMatrix.cpp `main`)

```
main(a1=query, a2=DB, a3=out): control.oneprocess(a1,a2,a3)
oneprocess:
  # ---- load the ONE query (lines ~249-266) ----
  qrow = quNumOfele(quline, qElecol)
  quryMatrix = new[qrow][qrow]
  formqMatrix(quryMatrix, qrow, inputfile, totalhandqury, totalrequire)
  # ---- stream the DB (lines ~277-408), per record ----
  while getline(intMline):
     header (.ssd) -> intMnumofele -> mrow, IntMnumEl, intElecol; alloc interActionM
     sheet (s...)  -> accumulate totalsheet
     matrix        -> getInterActionM(interActionM,...)
                      checksheetH(interActionM, mrow, totalsheet, pid, intElecol)   # record-dependent
                      searchM(intElecol, interActionM, mrow, qElecol, quryMatrix, qrow, totalpass, IntMnumEl)
                      selectMatrix(IntMnumEl, totalhandqury, totalrequire, totalpass, totalsheet, intElecol, interActionM)
                      printOuptfile(totalpass, IntMnumEl, pid, a3)   # writes a3/pdb<pid>.txt
```

Key fact: `checksheetH` and the record parse are **record-dependent only**;
`searchM` / `selectMatrix` / `printOuptfile` are the only **query-dependent** steps.

## Target shape

```
main(a1=query-manifest, a2=DB, a3=out-root)
oneprocess:
  # ---- load ALL queries once ----
  vector<QuerySpec> Q
  for path in read_manifest(a1):
     q.name  = basename(path, ".query")
     q.qrow  = quNumOfele(quline, q.qElecol)
     q.quryMatrix = new[qrow][qrow]; formqMatrix(q.quryMatrix, q.qrow, f, q.hand, q.require)
     mkdir a3/q.name              # once
     Q.push_back(q)
  minq = min(q.qrow for q in Q)   # for the cheap size early-exit
  # ---- stream the DB ONCE ----
  while getline(intMline):
     header/sheet parsing unchanged -> per record: mrow, IntMnumEl, intElecol, interActionM, totalsheet, pid
     on matrix line (record complete):
        if mrow < minq: free; continue          # nothing can match
        checksheetH(interActionM, mrow, totalsheet, pid, intElecol)   # ONCE per record
        for q in Q:
           if q.qrow > mrow: continue            # size early-exit (bonus)
           totalpass.clear()
           searchM(intElecol, interActionM, mrow, q.qElecol, q.quryMatrix, q.qrow, totalpass, IntMnumEl)
           selectMatrix(IntMnumEl, q.hand, q.require, totalpass, totalsheet, intElecol, interActionM)
           printOuptfile(totalpass, IntMnumEl, pid, a3/q.name/)       # per-query output dir
        free interActionM
```

## Concrete changes

1. **`QuerySpec` struct** — bundle the per-query state that is currently loose
   locals in `oneprocess`: `char** quryMatrix; int qrow; vector<elecol> qElecol;
   vector<handness> hand; vector<require> require; string name;`.

2. **`loadQuery(path) -> QuerySpec`** — extract oneprocess lines ~243-266 verbatim
   into a function returning a `QuerySpec` (open file, read 2 header lines,
   `quNumOfele`, alloc `quryMatrix`, `formqMatrix`). No logic change.

3. **Manifest input** — arg1 becomes a *manifest* (one query path per line) or a
   directory. **Backward compatible:** if arg1 ends in `.query` / is a single
   file, wrap it as a 1-element manifest (existing callers keep working).

4. **Refactor the reader (277-408)** — split into (a) *parse record* (header +
   sheets + matrix → the record state) and (b) *the query loop*. The record parse
   is exactly today's code; the inner search block (327-354 today) moves into
   `for q in Q`. `checksheetH` is called **once per record**, before the loop.

5. **Per-query output** — `printOuptfile` already takes a `path`; pass
   `a3 + "/" + q.name + "/"`. Create each `a3/q.name/` once at load. Mirrors the
   sweep's current `$OUT/hits/$NAME/` layout, so downstream tooling is unchanged.

6. **Size early-exit (bonus)** — `if (mrow < q.qrow) continue;` skips queries a
   record is too small to satisfy; `if (mrow < minq) continue;` skips the whole
   record (≈43% of the DB for S5). This is the disabled `a.size()<b.size()` guard
   from `searchM` (line ~1691), applied at the caller where it actually saves the
   `checksheetH` + matrix work.

7. **Memory** — `QuerySpec.quryMatrix` allocated once at load, freed at end.
   `interActionM` allocated per record, freed after the query loop. 6,336 queries
   × ~(qrow² + small vectors) ≈ a few MB resident — negligible.

## Parallelism

The hardening already makes any scan desync-safe, so chunking is no longer needed
for *correctness* — only for *throughput*. Recommended sweep driver:

- Partition the DB into `C` chunks (`db_validate`-clean boundaries, or just record
  counts). Run **C SLURM tasks**, each `searchmatrix <all-queries-manifest>
  <chunk_c> <out>`; union the per-query hit dirs. Total DB parse = **one** pass
  (distributed across C tasks) instead of 6,336. Pick `C` ≈ available cores.

This replaces the current "6,336 tasks × full-DB parse" with "C tasks × (DB/C
parse) × all-queries", i.e. the DB is parsed once total rather than 6,336 times.

## Validation

- **Hit-identity:** on the uniform-5000 sample and on 1tim, the inverted binary's
  per-query hit set must be **byte-identical** to the current per-query runs
  (diff `out/<qname>/` trees). This is the acceptance gate.
- **Determinism:** query order must not affect any query's hits (each query writes
  its own dir; no shared mutable state between queries except the read-only record
  — verify `searchM`/`selectMatrix` take the record by const/no-mutate).
- **Regression:** the existing `run_tests.sh` T1-T3 (byte-identical genmat, known
  hit, 1tim/1ubq end-to-end) must still pass with the inverted engine behind the
  single-query back-compat path.

## Effort / risk

- ~1 focused day. Mechanical refactor (struct + function extraction + loop swap);
  no algorithm change. Main risk is shared mutable state leaking between queries —
  mitigated by the hit-identity gate. Keep the single-query path as the default
  build; gate the manifest path so the change is opt-in until validated.

## Expected payoff

- **I/O:** 16 TB → 2.6 GB of DB reads for a full S5 sweep.
- **Compute:** parsing done once instead of 6,336× (~2× total, workload-dependent),
  plus ~43% of records size-skipped for S5 before any matching.
- Stacks with the committed `mkdir` (12×) and grep-removal fixes.
