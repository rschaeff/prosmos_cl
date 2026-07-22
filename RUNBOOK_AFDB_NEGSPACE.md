# RUNBOOK: AFDB negative-space sweep on HGD

End-to-end instructions for running the 800-query ProSMoS negative-space
sweep against AFDB v4 structures on the HGD cluster. Owner: Gray.

Quick context: we identified 25 S5 skeletons × 32 H/E typings = **800
queries** that find zero matches across 19,015 ECOD manual-rep PDB
domains (and so far ~387 of ~440 also find zero in 707k F70 PDB domains,
sweep ongoing on leda). The question this sweep answers: **which of these
PDB-absent motifs are realized in AFDB-predicted structures?**

Background reading (in this repo):
- `enum/docs/coverage_gaps.md` — why these 800 queries
- `enum/docs/positive_controls.md` — methodology validation (FHB + RLM)
- `results/v4_ecod_manual_reps_summary.tsv` — the manual-rep sweep results

The 800 queries are pre-staged at **`enum/negspace_queries/`** in this
repo — no need to regenerate.

---

## What you'll be running

Three sequential pipelines:

1. **Build dependencies** — compile `searchmatrix`, `generateMatrix`, set
   up the `palsse_cl` Python package.
2. **Build the AFDB metamatricesDB** — run PALSSE + generateMatrix over
   AFDB structures, concatenate per-structure matrices into one DB file.
3. **Run the sweep** — submit 800 queries as a SLURM array against the
   AFDB metamatricesDB, merge per-task TSVs into a summary.

The sweep itself (step 3) is the fast step (hours). The DB build (step 2)
is the heavy step and its cost scales linearly with how many AFDB
structures you include.

---

## Sizing decision (read this first)

AFDB v4 has ~214M predicted structures. Building a full-AFDB
metamatricesDB at our observed ~2-3 sec/structure on 8-way parallel
tasks would be ~1.7M task-hours — not realistic.

Choose one:

- **AFDB cluster reps** (recommended): one representative per AFDB
  sequence cluster. ~52M reps for AFDB clusterRep30 or ~5M for
  clusterRep50; final choice depends on what's pre-indexed on HGD.
- **Lab-curated AFDB subset** (e.g. existing `afdb_pyecod_mini_bench`
  set on leda is a starting template — confirm equivalent on HGD).
- **Taxonomy-subset** (e.g. archaea-only, ~600k structures) if you want
  a fast first pass before committing to broader sampling.

Whichever you pick, **document the chosen subset's size + selection
criteria** in the OUT dir (`README.md` next to the manifest). The
negative-space interpretation depends on knowing what fraction of AFDB
you actually searched.

---

## Step 1 — Build dependencies on HGD

Both repos are at `/data/ECOD0/html/distributions/prosmos_afdb_negspace/`
(also on the leda side at `~rschaeff/dev/{prosmos_cl,palsse_cl}` if you
need to grab updates).

```bash
# On HGD, in a scratch workdir of your choice:
cd ~/your_workdir
tar xzf /path/to/prosmos_cl.tar.gz
tar xzf /path/to/palsse_cl.tar.gz
```

### 1a — searchmatrix and generateMatrix

```bash
cd prosmos_cl/searchMatrix/src
g++ -O2 -DSILENT searchMatrix.cpp -o ../build/searchmatrix
# Warnings only (format-string + unused-result), no errors expected.

cd ../../generateMatrix/src
# generateMatrix is MPI-linked. Use the system mpicxx, NOT the conda one.
which mpicxx       # confirm it points to /usr/bin/mpicxx or HGD equivalent
mpicxx generMatrix.cpp -o ../build/generateMatrix
```

`-DSILENT` is **essential** for searchmatrix — without it, the binary
writes tens of millions of debug-print lines per query, which crushes
any NFS log path. See commit `80700c1`.

### 1b — palsse_cl Python package

```bash
cd ../../../palsse_cl
python3 -m venv .venv
source .venv/bin/activate
pip install biopython numpy
pip install -e .
# Sanity:
python3 -c "import palsse; print(palsse.__file__)"
```

The PALSSE Python implementation is a byte-for-byte port of the Grishin
lab's original Python 2 oracle. The 4 commits at the tip of palsse_cl
(`4b30449`, `81ef0d8`, `79ba278`, `6fe1391`, `9224daf`) fix specific
divergences from the oracle; the codebase is stable as of this handoff.

---

## Step 2 — Build the AFDB metamatricesDB

Pattern is the same as the leda ECOD builds — see
`scripts/ecod_db_build/{submit.sh,array.sbatch,process_chunk.sh}` for
reference. The leda scripts pull a manifest from the ECOD postgres DB;
on HGD you'll build the manifest from your chosen AFDB subset directly
(no DB pull needed).

### 2a — Build the manifest

```bash
OUT=$HOME/work/afdb_negspace_db
mkdir -p $OUT/{manifest,chunks,fragments,logs,work}

# Replace this with however you enumerate your AFDB subset paths.
find /path/to/afdb/subset -name '*.pdb' -type f | sort > $OUT/manifest/paths.txt
N_PATHS=$(wc -l < $OUT/manifest/paths.txt)
echo "$N_PATHS structures"
```

### 2b — Split into chunks and submit

The `scripts/ecod_db_build/process_chunk.sh` script handles one chunk:
for each PDB path in the chunk, run PALSSE to get SSE definitions, then
run generateMatrix to produce the per-structure interaction matrix, then
append to a fragment file. The chunks are merged into one
metamatricesDB at the end.

```bash
# Adapt scripts/ecod_db_build/submit.sh — the DB_FILTER / PG* env vars
# are for the ECOD pull on leda. On HGD just set:
#   - SKIP_DB_PULL=1 (or just remove the psql block in your local copy)
#   - the manifest must already be at $OUT/manifest/paths.txt
#
# Or write a thin HGD-specific submit wrapper. Either way you want:
#   - N_CHUNKS chosen so each chunk is ~5-10k structures (so each task
#     finishes within --time=02:00:00)
#   - --cpus-per-task=8 (process_chunk.sh forks 8 PALSSE+generateMatrix
#     workers)
#   - --exclude any bad-MPI nodes on HGD (see leda exclude pattern in
#     submit.sh — generateMatrix's MPI_Init fails on cluster nodes built
#     without --with-pmi)

bash scripts/ecod_db_build/submit.sh  # adapted as above
```

The leda exclude pattern was `leda23,leda25,...,leda45` (odd-numbered
older hardware generation). Test one chunk on a few HGD nodes before
launching the full sweep; if generateMatrix dies in MPI_Init on some
nodes, add them to `--exclude`.

### 2c — Filter post-merge

generateMatrix occasionally emits malformed entries when input PDBs have
multi-character chain IDs (e.g. "40") — these produce 5-digit
coordinate fields that searchmatrix later crashes on with
`std::out_of_range` in `basic_string::assign`. The leda v3 build hit
~84 such entries in 707k inputs (~0.01%).

After the array+merge job lands, filter the DB:

```bash
# Find offending entries: any row with >4-digit residue numbers in the
# coord block. The leda v3 .clean filter pattern (awk skip) is in the
# git history at commit 120f836; copy it here.

awk '...' $OUT/metamatricesDB > $OUT/metamatricesDB.clean
grep -c '\.ssd' $OUT/metamatricesDB.clean
```

**Always use the `.clean` version for the sweep.** Otherwise expect
100% rc=134 task failures partway through every sweep.

### 2d — Sanity-check the DB

Before the 800-query sweep, run **one positive control** to confirm the
DB is sound. Pick the RLM α-variant winner:

```bash
mkdir -p /tmp/sanity_hits
$PROSMOS_HOME/searchMatrix/build/searchmatrix \
    enum/negspace_queries/../queries_typed/s5/s5-0098-0021.query \
    $OUT/metamatricesDB.clean \
    /tmp/sanity_hits > /tmp/sanity.log 2>&1
echo "hits: $(find /tmp/sanity_hits -type f | wc -l)"
# Expect: thousands of hits on AFDB (RLM is ubiquitous).
# Zero hits = your DB build is broken; investigate before continuing.
```

(You'll need to regenerate the full corpus for `s5-0098-0021` via
`python3 enum/scripts/generate_enumerated_queries.py --dims 5`, or just
copy the file directly from leda at
`~rschaeff/dev/prosmos_cl/example/ssp_enumerated/queries_typed/s5/s5-0098-0021.query`.)

---

## Step 3 — Run the 800-query sweep

```bash
QDIR=$PWD/enum/negspace_queries \
DB=$OUT/metamatricesDB.clean \
OUT=$HOME/work/afdb_negspace_results \
CONCURRENCY=128 \
TIMELIMIT=06:00:00 \
EXCLUDE_NODES="<comma-separated bad-node list from step 2b, if any>" \
bash scripts/slurm_search/submit.sh
```

What this does:
- Submits one SLURM array job with 800 tasks (well under the typical
  MaxJobCount=10000)
- Each task picks one query from `QDIR`, runs searchmatrix against the
  DB, writes hit files to `$OUT/hits/<query>/`, writes a per-task TSV
  to `$OUT/parts/<line>.tsv`
- A dependent merge job concatenates the part TSVs into
  `$OUT/summary.tsv`

Monitor:
```bash
# One-shot snapshot
bash scripts/slurm_search/watch_failures.sh $OUT

# Background loop logging every 5 min
nohup bash -c "while true; do bash scripts/slurm_search/watch_failures.sh $OUT; echo; sleep 300; done" >> $OUT/failures.log 2>&1 &
```

Expected runtime: at the leda observed rate of 26 min mean per query
against a 707k DB, 800 queries × 26 min / 128 concurrent ≈ **3 hours
wall**. AFDB DB will likely be larger, scaling linearly with entry
count.

### If any tasks fail (rc != 0)

`watch_failures.sh` will surface them with the searchmatrix log tail.
Most likely cause: a malformed DB entry that survived the post-merge
filter. Identify the offending entry, add to the filter, regenerate
the `.clean` DB, re-run the failed line numbers (see leda
`ecod_search_v4_retry` pattern).

---

## Step 4 — Interpret results

The headline question for each query: **was it zero-hit?**

```bash
# Zero-hit queries — the AFDB-confirmed PDB-negative-space:
awk -F'\t' 'NR>1 && $4==0' $OUT/summary.tsv > $OUT/still_zero.tsv
wc -l $OUT/still_zero.tsv

# Lit-up queries — motifs PDB doesn't have but AFDB does:
awk -F'\t' 'NR>1 && $4>0 {print}' $OUT/summary.tsv | sort -t$'\t' -k4,4 -rn > $OUT/lit_up.tsv
head $OUT/lit_up.tsv
```

Cross-reference against the leda PDB-side results to identify the four
interesting categories:

| Category | PDB hits | AFDB hits | Meaning |
|---|---|---|---|
| Truly absent | 0 | 0 | Motif unrealized in nature even at AFDB scale |
| **AFDB-only** | **0** | **>0** | **Headline result — motifs evolution made but PDB missed** |
| AFDB-loss | >0 | 0 | Shouldn't happen if AFDB ⊃ PDB content; flag for diagnosis |
| Common | >0 | >0 | (expected for any negspace queries that lit up in the leda 707k sweep) |

The leda zero-hit set after the 707k sweep is in
`/data/ECOD0/html/distributions/prosmos_afdb_negspace/leda_zero_hit_queries.txt`
(updated as the leda sweep finishes — currently partial, will be final
within ~2 weeks).

---

## Files of interest in this repo

- `scripts/slurm_search/submit.sh` — sweep launcher
- `scripts/slurm_search/array.sbatch` — per-task wrapper (uses local
  /tmp for hot path — see commit `821db38`)
- `scripts/slurm_search/watch_failures.sh` — sweep monitor
- `scripts/slurm_search/recover_s6.sh` — throttled launcher pattern,
  use if AFDB sweep needs to be larger than MaxJobCount limit
- `scripts/ecod_db_build/{submit.sh,array.sbatch,process_chunk.sh}` —
  metamatricesDB build pipeline (ECOD-flavored on leda; adapt to AFDB)
- `enum/docs/{coverage_gaps.md,positive_controls.md}` — methodology
  context
- `searchMatrix/src/searchMatrix.cpp` — the SILENT patch is committed
  here; rebuild from this source, not the legacy `Linux/` binary

## Contacts

- Dustin Schaeffer (dustin.schaeffer@gmail.com) for question on
  methodology, enum design, or what counts as a "real" PDB absence
- The leda-side in-flight stage 1 sweep state is in
  `~rschaeff/work/prosmos_2026/ecod_search_v4_negspace{,_part2}/` — ask if you need
  the survivor list before mine completes
