# Session handoff — 2026-07-11 (geometric SCC-2 + searchmatrix undercount + corrected sweep)

Branch: `s5-matrix-temporal-promiscuity`. All work committed except the live SLURM
sweep (running). Read this first after a compaction.

## LIVE STATE — a full-scale search is RUNNING (needs follow-up on completion)

- **SLURM jobs:** `622256 622257 622258 622259 622260 622261 622262` (7 chained
  array chunks, 1000 tasks each, %200 concurrency, `--time=08:00:00`) + merge
  `622263`. Watch: `squeue -j 622256 ...`.
- **What:** the 6,336 graph-198 S5 queries (`queries_graph198_alltypings/s5`,
  byte-identical to the original `s5_full_afdb`) re-run on the **hardened +
  fast** searchmatrix against the full 4.92M-record DB
  (`work/prosmos_2026/afdb_db/metamatricesDB.clean`). Output →
  `work/prosmos_2026/s5_full_afdb_hardened/`.
- **Why:** the original `s5_full_afdb` (0.31% hitters / 15,327 domains) was a
  searchmatrix **undercount** (see below). This is the corrected sweep.
- **Speed/ETA:** ~400–500s/query (node-local DB staging), ~46 concurrent
  (cluster contention), **~15h**, 0 failures so far. 8h/task limit → no timeout
  truncation.
- **Completion handler:** `/tmp/sweep_done.sh` (log `/tmp/sweep_done.log`) waits
  for `s5_full_afdb_hardened/summary.tsv`, then writes:
  - `distinct_hitters.txt` — distinct hitting structures (= corrected dark fraction),
  - `failed_queries.txt` — any `rc!=0` (timeouts) to retry.
  It prints `SWEEP_ANALYSIS_DONE`. (If the handler died in the compaction, re-run
  the same analysis by hand from `summary.tsv` + `hits/`.)

### On completion — DO THIS
1. **Corrected headline:** `wc -l distinct_hitters.txt` = true hitter count.
   Compare to original 15,327 (0.31%). Expected ~12–25% (my uniform-5000 probe
   gave 24.6% all-typings). Dark fraction = 1 − hitters/4,921,931.
2. **Retry pass:** if `failed_queries.txt` non-empty, re-run those with 8h via
   `scripts/slurm_search/submit.sh` (build a QDIR of just those `.query` files).
3. **Negative space:** the 800/800 zero-hit negspace result (`enum/negspace_queries/`,
   memory `project_negspace_finding`) was computed on the BUGGY engine — recheck
   which enumerated skeletons still never hit on the corrected sweep.
4. Aborted NFS-bound run's partial output is at `s5_full_afdb_hardened_nfsbound_aborted/`
   (can delete).

## What was done this session (commits, newest first)

| commit | what |
|---|---|
| `3784521` | slurm_search: **stage DB to node-local /tmp** (5–6× faster; NFS 16TB→130GB) |
| `2a06356` | searchmatrix: **loop-inversion spec** (`searchMatrix/loop_inversion_spec.md`) |
| `c1eec66` | searchmatrix: drop per-query full-DB grep sanity count |
| `74a5bcb` | searchmatrix: **12× speedup** — hoist per-record `mkdir` out of checksheetH |
| `65b7485` | searchmatrix: **harden DB reader** (line-type dispatch) + validator overflow guard |
| `0900f74` | db: fix undercount proof to same-query clean-vs-dirty (retract confounded subset claim) |
| `985c339` | db: **validator** (`scripts/db_validate.py`) + full-DB undercount diagnosis |
| `8528368` | enum: Phase 4 doc — fresh graph-198 queries == qian queries byte-identical |
| `dc2044b` | enum: **Phase 4 impact probe** (`scc2_phase4_impact.md`) |
| `6ce5590` | enum: **Phase 3** paper-faithful query generation (labeling map) |
| `dda980c` | enum: read S5 **minimal interaction** sets off Fig S3 |
| `e8bc425` | enum: read S5 **minimal handedness** lists (oracle) |
| `7d3fc48` | enum: **geometric SCC-2** (Fig-S3 hex-congruence) — paper-faithful S5 = **140** |

## Key results (durable)

### Geometric SCC-2 root fix (enum/)
- The old graph-based SCC-2 admitted bent variants; the paper's is **hex-congruence
  to Fig-S3 grids**. Fixing it: **S5 198 → 140** (opt-in `set_scc2_mode("geometric")`
  in `compactness.py`; default stays graph). Dropped 58 = bent grid-e variants, 0
  added. Code: `enum/src/ssp_enum/geom_scc2.py`.
- **Minimal handedness** (from oracle IA-S5.txt): d=0, e=2 `{(1,4,5),(2,3,5)}`,
  g/h=6 (full), f=7 + 2 conditional `{(1,4,5),(3,4,5)}` firing on typing HHHEE.
- **Minimal interactions** (off Fig S3): d/e/f/g-h required + grid-e `≥1 mandatory`
  broken pair `(2,5)‖(3,5)`. Verified partition of each grid's lattice edges.
- **Phase 3:** `skeletons_to_records(..., paper_faithful=True)` → **708 records**
  (grid-e emits 2 variants; f conditional). `reference_labeling()` maps survivors
  to Fig-S3 labels. Residual: 28 survivors (21 f, 7 g/h) have a mandatory triple
  coplanar under every labeling (sequence-parity degeneracy) → flagged
  `mandatory_handedness_gap`; likely an RCC/handedness-equivalence dedup artifact
  (Phase 5). Docs: `enum/docs/scc2_*_result.md`, `scc2_geometric_plan.md`.

### Phase 4 impact probe
- **Geometric SCC-2 is hit-neutral** (graph-198 == geom-140 hit sets, on 1tim and
  a sample). So the 198→140 basis change doesn't move darkness.
- The **paper's typing rule** (Rule 1: "≥4 collinear ⇒ β-sheet", Chitturi SI) is
  what drives paper-faithful darkness — it **darkens even 1tim** (38 S5 hits → 0),
  because 1tim matches grid-e queries typed with a helix on the collinear-4 run.
  Separate "geometric darkness" from "rule-induced darkness".

### searchmatrix full-DB UNDERCOUNT (the big one — see `enum/docs/searchmatrix_fulldb_undercount.md`)
- Trigger: Qian's cluster gave S5 hit rate **11%** vs our published **0.31%**.
- **Not** leda, **not** the cwd-race (harness isolated), **not** the queries
  (fresh == qian byte-identical). It's the engine: the reader alternated
  header/matrix by `readLincon%2` and **5 orphan records** (a `sheet`+`matrix`
  block with no header, first at DB record ~67k) desynced it; incomplete recovery
  **silently dropped ~95% of the DB**.
- Proof (same query, same DB): s5-0096 = **125** dirty vs **2,178** cleaned; chunk 0
  alone = 125 = whole dirty scan.
- **Fixes:** hardened reader (line-type dispatch, `65b7485`) — no single bad record
  can desync; `scripts/db_validate.py` validates/`--clean`s a DB. generateMatrix
  needs **no** change (format self-consistent: 68-byte/SSE stride matches; 0 field
  overflows; added a forward guard).
- **True S5 rate** ≈ 12–25% (uniform-5000 probe: 24.6% all-typings / 12.7%
  paper-faithful), consistent with Qian's 11%. The 0.31% is undercounted ~40–80×.
- **CONSEQUENCE:** the darkness fraction, D1/D2, dark-gallery, negspace, and
  rarefaction numbers all sit on the buggy full-DB scan and are quantitatively
  wrong (darkness overstated). The live sweep is the correction.

### Performance (searchmatrix)
- `checksheetH` ran `system("mkdir")` **per record** (fork+exec ×4.92M) = ~92% of
  runtime → hoisted to once/process = **12×** (`74a5bcb`).
- Removed per-query full-DB `grep` sanity count (`c1eec66`).
- **Loop-inversion spec** (`searchMatrix/loop_inversion_spec.md`): parse DB once,
  all queries per record → DB parsed 1× not 6,336×; ~2× more + a size early-skip.
  Not yet implemented; the biggest remaining lever.
- Node-local DB staging (`3784521`) removed NFS read contention (5–6×).

## Environment / how-to
- Enum python: `/sw/apps/Anaconda3-2023.09-0/bin/python`, `sys.path.insert(0,'enum/src')`.
- Hardened searchmatrix build: `cd searchMatrix/src && g++ -DSILENT -O0 -w
  searchMatrix.cpp -o ../build/searchmatrix` (`-O0` required; `-O2` risks
  miscompile and only ~9%). Production binary already rebuilt with all fixes.
- Search harness: `scripts/slurm_search/{submit.sh,array.sbatch}` — `QDIR OUT DB
  CONCURRENCY TIMELIMIT ./submit.sh`. array.sbatch now stages DB to node /tmp.
- Probe harness (small DBs): `work/prosmos_2026/run_probe.sh <qdir> <db> <out>`
  (query paths MUST be absolute — it cd's into an isolated dir).
- Query sets in `work/prosmos_2026/`: `queries_graph198_alltypings/s5` (6336, =orig),
  `queries_geom140_alltypings/s5` (4480), `queries_paper_faithful_140/s5` (708).
- Uniform sample DB: `work/prosmos_2026/sample_afdb_uniform_5000.db` (seed 20260711).

## Pending / next
1. **Live sweep completion** (above) — corrected darkness numbers + retry + negspace recheck.
2. **Phase 5** (enum): RCC/handedness-equivalence dedup to resolve the 28-residual and reconcile the 648 oracle.
3. **Loop inversion** (searchmatrix): implement the spec — the real perf/IO win.
4. Update memories `project_negspace_finding` (result now suspect — undercounted engine) and `project_s5_oracle_gap` (S5 now 140 geometric, was 198) once corrected numbers land.
