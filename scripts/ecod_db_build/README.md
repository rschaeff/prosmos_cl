# ECOD metamatricesDB build

Build a ProSMoS-format `metamatricesDB` from ECOD domain PDB structures.
The build runs as a SLURM array — one task per chunk of domains —
followed by a dependent merge job that concatenates per-task fragments
into the final DB.

## Pipeline per domain

```
ECOD .pdb  ──►  PALSSE  ──►  .ssd  ──►  generateMatrix -os  ──►  .out
```

Per-chunk:
```
N .out files  ──►  cat  ──►  fragment_<task_id>.frag
```

Final:
```
N fragments  ──►  cat (sorted by task_id)  ──►  metamatricesDB
```

## Files

| File | Role |
|---|---|
| `submit.sh` | Pulls manifest from `ecod_commons.derived_files`, splits into chunks, submits the array + merge. The orchestration entry point. |
| `array.sbatch` | Per-task script. Reads `$SLURM_ARRAY_TASK_ID`-th chunk file and runs `process_chunk.sh`. |
| `process_chunk.sh` | The actual per-chunk worker: parallel PALSSE → parallel generateMatrix → concatenate. Not SLURM-specific; runnable standalone for local testing. |
| `README.md` | This file. |

## Usage

```bash
OUT=$HOME/work/ecod_db \
N_CHUNKS=200 \
PARALLEL=8 \
TIMELIMIT=02:00:00 \
./submit.sh
```

`$OUT` must be NFS-visible from compute nodes.

After the merge job completes, `$OUT/metamatricesDB` is ready for
`searchmatrix`.

## Environment knobs

| Var | Default | Note |
|---|---|---|
| `OUT` | required | Output root. NFS-visible. |
| `N_CHUNKS` | 200 | Number of array tasks. With ~1.77M ECOD domains and 200 tasks, each task processes ~9K domains. |
| `PARALLEL` | 8 | xargs concurrency per task (also `--cpus-per-task` for SLURM). |
| `TIMELIMIT` | `02:00:00` | Per-task `--time`. 2h is generous for 9K domains at ~0.7s each. |
| `DB_FILTER` | complete PDB-source | SQL WHERE clause for manifest selection. See below for subsets. |
| `CONCURRENCY` | 64 | Max concurrent array tasks (`%N` in `--array`). |

## Subsetting the input

The default `DB_FILTER` selects **all complete PDB-source domains** (~1.77M
at v294.2). To build a smaller DB for testing or representative work,
override `DB_FILTER`. Examples:

### F40 representatives only (~150K domains)
```bash
DB_FILTER="
  file_type_id=2 AND status='complete' AND domain_source_type='pdb'
  AND ecod_uid IN (
    SELECT ecod_uid FROM ecod_commons.domains
    WHERE f40_representative = true
  )
"
```

(Adjust column names to match the actual schema — verify with
`\d ecod_commons.domains` first.)

### A random sample for smoke testing
```bash
DB_FILTER="
  file_type_id=2 AND status='complete' AND domain_source_type='pdb'
  AND random() < 0.01
"
```
(About 17K random domains.)

### Including AlphaFold-predicted structures
```bash
DB_FILTER="file_type_id=2 AND status='complete'"   # no domain_source_type filter
```

## Sizing notes

Empirical timing from a 100-domain smoke run:
- PALSSE: ~85ms per domain (single-threaded), wall-time ~10s with 8 workers
- generateMatrix: ~560ms per domain (single-threaded), wall-time ~7s with 8 workers
- **Total per domain (8-way parallel): ~85ms wall-time**

Extrapolation:
- 1.77M domains / 1600 cores (200 tasks × 8 workers): ~95K core-seconds = ~16 hours wall-time per task. **Too long for the default 2h `TIMELIMIT`.**
- 1.77M domains, 8 cores per task, 9K domains per task: 9000 × 0.7s / 8 = 13 min wall-time per task. **OK at 2h limit.**

The default settings (200 tasks × 9K domains × 8-way parallel × 2h limit)
give comfortable headroom. Each task should finish in 15-30 min.

## Output format

The final `metamatricesDB` is a concatenation of per-domain entries:

```
<basename>.ssd <NSSE><SSE_line_1><SSE_line_2>...<SSE_line_N>
sheet <id> <count> <indices>...
sheet <id> <count> <indices>...
*<upper-triangular interaction matrix as packed string>*
```

(SSE lines are concatenated on the same line as the filename; only sheet
records and the matrix line have their own newlines.)

`searchmatrix` reads this format directly. Concatenation order is sorted
by task ID then by filename within each task — reproducible across runs
given the same manifest.

## Validation after build

After the merge job lands a `metamatricesDB`:

```bash
# Entry count (should match manifest line count, minus any PALSSE failures)
grep -c '\.ssd' $OUT/metamatricesDB

# Quick search-binary sanity check
mkdir -p $OUT/sanity
$REPO/searchMatrix/build/searchmatrix \
    $REPO/example/ssp_design_targets/queries_enum/02-5-311-0-0.query \
    $OUT/metamatricesDB \
    $OUT/sanity/ \
    2>&1 | tail
ls $OUT/sanity/ | wc -l   # hit count
```

## Re-running failed tasks

If some array tasks failed (timeouts, node issues, etc.), the
`afterany` merge runs anyway and includes only the successful fragments.
To fill in gaps:

```bash
# Find which task IDs failed (no fragment produced)
seq 1 $N_CHUNKS | while read i; do
  [ -f $OUT/fragments/fragment_${i}.frag ] || echo $i
done > $OUT/manifest/failed_tasks.txt

# Re-submit only those, e.g. as --array=3,17,42 — the chunk files are
# already on disk so the per-task script just re-runs them. Bump
# TIMELIMIT first if the failure was a timeout.
sbatch --array=3,17,42 \
    --time=04:00:00 \
    --export=CHUNK_DIR=$OUT/chunks,WORK_ROOT=$OUT/work,FRAG_DIR=$OUT/fragments,PARALLEL=8 \
    array.sbatch
```

Then re-cat fragments → metamatricesDB once the retries finish.

## Known caveats

- **`/tmp` on leda is not NFS-shared**; the array.sbatch uses `WORK_ROOT`
  (under `$OUT/work`) for intermediate `.ssd`/`.out` files. Per-task
  cleanup removes these at task exit.
- **PALSSE's HELIX/SHEET column positions** had an off-by-one vs the PDB
  spec (chain at col 18 instead of col 19) in commit-pre-95e8bcd. Older
  PALSSE checkouts will produce 0-SSE matrix entries — `pip install -e .`
  the latest `~/dev/palsse_cl/` before running this build.
- **generateMatrix's `-ds` directory mode** scales pathologically (>4
  min for 100 files). We use `-os` per-file in a parallel xargs pool
  instead — much faster.
- **AlphaFold-predicted vs PDB-deposited**: by default this build uses
  only PDB-source domains. AlphaFold-source domains are also in
  `ecod_commons.derived_files` but represent inferred structures that
  ProSMoS searches may not interpret the same way as crystal/cryo-EM
  inputs. Override `DB_FILTER` to include them if needed.
