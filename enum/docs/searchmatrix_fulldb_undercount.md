# searchmatrix undercounts on the full 4.92M-record DB (parser desync)

Triggered by Qian's report: profiling on her cluster gave S3/S4/S5 hit rates of
**58% / 27% / 11%**, ~20-40× above our published AFDB S5 rate of **0.31%**
(15,327 hitting domains / 4.92M, from `s5_full_afdb/`). Investigation of whether
this is (a) a leda/node issue or (b) a bug in our search.

## Answer: (b), but NOT the cwd-race bug — searchmatrix under-scans the large DB

The current search harness (`scripts/slurm_search/array.sbatch`) is correctly
isolated: each query runs in its own `/tmp/sm_${LINE}_$$/run` cwd with its own
`sheetbug/`. No worker-collision race. Our current engine on **small** DBs agrees
with Qian (see below). The defect is in scanning the **full 2.6 GB / 4.92M-record**
DB.

### The airtight proof (same query, same structures, clean vs dirty DB)

Same current binary (`qian_package/.../bin/searchmatrix`), same graph-198 query
`s5-0096-0000`, over the same 4.92M structures:

| DB scanned | hits |
|---|---|
| **dirty** full DB, single scan | **125** |
| dirty full DB split into 20 chunks, summed | **2,098** |
| **cleaned** full DB (5 orphan records removed), single scan | **2,178** |

The only difference between the 125 and the 2,178 is the **5 orphan records**
`db_validate` removed — so those 5 records cost **~94%** of this query's hits.
Chunking (fresh parse per chunk) independently recovers the same ~2,098. The
full-DB count is deterministic (re-runs reproduce 125; `s5-0132-0000` reproduces
its recorded 7,372 exactly, 871–1049 s, clean `rc=0`), so the original
`s5_full_afdb` numbers are faithfully reproduced — and are undercounts.

Even sharper — the per-chunk breakdown: **chunk 0 alone (records 0–246k) = 125**,
*exactly* the whole dirty-DB count. The scan dies at the first orphan (record
~67k) and silently drops chunks 1–19 (**95% of the DB**).

(An earlier "subset of 424 out-hits the full DB's 125" framing was retracted: the
424 came from a `geom-140`-numbered query, a *different* skeleton than the
`graph-198` `s5-0096`. The clean-vs-dirty comparison above is same-query and
needs no cross-set comparison.)

### Mechanism: parser desync from a HANDFUL of orphan records

searchmatrix's reader (searchControl.h:277-324) parses records by an **even/odd
line counter**: non-`s` lines alternate HEADER (must contain `.ssd`) and MATRIX
(`*...`); `s`-prefixed sheet lines are skipped. It assumes every record is
exactly `[header] [0+ sheet lines] [one matrix line]`. The recovery on a detected
orphan is incomplete.

Validating the DB (`scripts/db_validate.py`, which simulates this exact parser)
finds only **5 truly malformed records** in 4.92M — orphan `sheet`+`matrix`
blocks with **no header line** (a generateMatrix glitch, first cluster at records
67,224-68,309). Those 5 desync the alternation, and the incomplete recovery drops
hits across the large downstream region. Loss is region/query-dependent:
`s5-0132` kept 7,372; `s5-0096` fell to 125. Chunking (fresh parse per chunk)
recovers most (`s5-0096`: 125 → 2,098 over 20 chunks; finer chunks recover more).
**5 bad records in 4.92M silently crippled the entire scan.**

Separately, `db_validate` flags ~6,231 benign 0-SSE records (count field 0, a
blank line serves as the empty matrix — searchmatrix stays in sync on these).

### A second, deeper vulnerability: fixed-width column overflow

The `.ssd`/matrix format is **fixed-column**. When a residue number reaches 4
digits (≥1000) or an SSE count reaches 2 digits, the value overflows its field
and **packs against the adjacent token** with no delimiter: `8HA1163 --1173`
(residue 1163) or `62HA` (count 62). **67,020 records (1.36%)** show this packing.
searchmatrix parses SSE elements at fixed offsets, so packed records are at risk
of mis-parse. (This packing also silently inflated an earlier "32k 0-SSE" count
to 5× its true value — a regex was fooled by the same overflow.) This is the
"vulnerable to truncation" concern: a validator catches it, but the format itself
has no framing/delimiters to prevent it.

## The true population S5 hit rate

Current engine on a **uniform-random 5,000-rep** sample (seed 20260711;
representative — 56.6% have ≥5 SSEs, matching the population sieve), scanning
desync-free small DBs:

| query set | hitters / 5000 | hit % | dark % |
|---|---|---|---|
| all-typings (old basis, 198≡140) | 1229 | **24.6** | 75.4 |
| paper-faithful (all 3 rules) | 635 | **12.7** | 87.3 |

Consistent with Qian's independent 11% (her test set differs). The published
**0.31%** is a searchmatrix full-DB-scan artifact, undercounted **~40-80×**. The
paper-faithful rules still halve hits (24.6% → 12.7%), reproducing the Phase-4
finding on the correct sample.

## Consequences

- **The darkness fraction is wrong quantitatively.** "99.7% dark / 0.55% hit" is
  a scan artifact. True S5 darkness is the **majority (~83-89%)** but not
  near-total, and ~15,327 "hitting domains" is really **hundreds of thousands**.
- **The negative-space (800/800 zero-hit) result must be rechecked** — those
  queries were run against the full DB and could have lost hits to desync.
- **The D1/D2, dark-gallery, and rarefaction analyses** all sit on full-DB search
  outputs and need redoing on de-synced-safe scans.
- **The geometric-SCC-2 / paper-faithful work is unaffected** — it is about the
  query set, not the DB scan. The Phase-4 relative decomposition used small-DB
  sample scans (desync-free) and stands.

## The fix (implemented)

**(A) Harden searchmatrix's reader — landed.** `searchControl.h` now dispatches
each line by **type** (`.ssd` → header, `s…` → sheet, else → matrix) instead of
the fragile `readLincon % 2` parity, guards the matrix branch on a pending header
(so an orphan matrix is *skipped*, not blindly consumed), and discards an
incomplete pending record when a new header arrives. A single orphan can no
longer desync the rest of the scan. Regression: identical 38-query hit-set on
1tim vs the old engine; on the **dirty** full DB the hardened engine recovers
s5-0096 from 125 → ~2,178 (≈ the cleaned-DB count) — no DB cleaning required.
Rebuild: `cd searchMatrix/src && g++ -DSILENT -O0 -w searchMatrix.cpp -o
../build/searchmatrix`.

**(B) Validator + clean — landed** (`scripts/db_validate.py`). Flags the 5
orphan records (and benign 0-SSE), and `--clean` re-emits a well-formed DB in
11 s. Either the hardened engine *or* a cleaned DB fixes the undercount; the
hardened engine is preferred (no preprocessing, robust to future bad records).

**(C) generateMatrix output — NO change needed.** Investigated the suspected
fixed-width "column packing": it is a **non-issue for this DB**. generateMatrix
emits a fixed 68 bytes/SSE (`%c%c%5s--%5s %4d` + 6×`%8.3f`) that *exactly* matches
searchmatrix's 68-byte parse stride, and **no field overflows** its width —
`db_validate` (with the new alignment guard) reports **0** field-overflow records
(max `NAME.ssd` = 24 ≤ 32; residues ≤ 4 digits fit `%5s`; 0 coordinates outside
`%8.3f` range). The earlier "67,020 packed records" was a regex artifact of the
trailing-space field format, which reads correctly at the fixed offsets. The
latent risk (a future run with names > 32 chars, residues > 5 digits, or
`|coord| ≥ 10000`) is now caught by the validator's `field-overflow` count, so
generateMatrix needn't be churned unless that guard ever fires.

## Reproduce

```
# hardened engine on the dirty DB (no cleaning) recovers the true count:
searchMatrix/build/searchmatrix \
  queries_graph198_alltypings/s5/s5-0096-0000.query afdb_db/metamatricesDB.clean out/   # ~2178
# validate / clean:
scripts/db_validate.py afdb_db/metamatricesDB.clean            # 5 orphan_block, 0 field-overflow
scripts/db_validate.py afdb_db/metamatricesDB.clean --clean cleaned.db
```
