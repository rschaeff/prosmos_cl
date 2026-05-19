# SLURM array runner for ProSMoS searchmatrix

Submit a ProSMoS motif search as an embarrassingly-parallel SLURM array — one
array task per query.

## Files

| File | Role |
|---|---|
| `array.sbatch` | The per-task script. Reads `$SLURM_ARRAY_TASK_ID`-th line of a query manifest, runs `searchmatrix`, writes a one-line TSV record. Not called directly. |
| `submit.sh` | Builds the manifest from a directory of `*.query` files and submits the array. Also queues a dependent merge job that produces a single `summary.tsv`. |
| `README.md` | This file. |

## Usage

```
QDIR=/home/rschaeff/dev/prosmos_cl/example/ssp_design_targets/queries \
OUT=$HOME/work/ssp_smoke_2010db \
./submit.sh
```

`$OUT` must be NFS-visible from compute nodes (do **not** use `/tmp` on leda —
its `/tmp` is not shared). `$OUT/hits/<query_name>/` will receive the per-hit
output files (PDB-id-keyed text files written by searchmatrix). `$OUT/summary.tsv`
will be populated by the merge job after the array completes.

Environment knobs:

| Var | Default | Note |
|---|---|---|
| `QDIR` | (required) | Directory of `*.query` files. All `*.query` are submitted. |
| `OUT` | (required) | Output root. NFS-visible, writable. |
| `DB`  | `/home/rschaeff/src/Prosmos/ProSMoS/metamatrixdb/metamatricesDB` | The 2010-era 306 MB metamatricesDB. Override for a freshly built DB. |
| `CONCURRENCY` | `64` | Max concurrent array tasks. Cap on busy clusters — for the full 1488-SSP rerun do not omit. |
| `SEARCH` | repo `searchMatrix/build/searchmatrix` | Path to the binary, picked up by `array.sbatch`. |

## What each array task does

1. Reads the `SLURM_ARRAY_TASK_ID`-th query path from `$QUERY_LIST`.
2. Creates a per-task work dir at `$OUT/work/<task_id>/run/` and `cd`s into it.
   This is **load-bearing**: `searchmatrix` unconditionally writes debug output
   to `<cwd>/../sheetbug/` and silently `exit(0)`s if it can't open
   `<cwd>/../sheetbug/total.txt`. Per-task work dirs avoid both the silent
   failure and cross-task collisions.
3. Runs `searchmatrix <query> <DB> $OUT/hits/<name>/`, logging stdout/stderr
   to `$OUT/logs/<name>.stdout`.
4. Records `name\truntime_sec\texit_code\thit_count` to
   `$OUT/parts/<task_id>.tsv`. The merge job concatenates these parts into
   `$OUT/summary.tsv` after the array finishes.
5. Cleans up the per-task `sheetbug/` debug output (large, useless after run).

## Why this exists (vs. running `searchmatrix` in a shell loop)

- 11 queries against the 2010 DB take ~10–60 min each sequentially. 1488 SSPs
  (full Chitturi-2016 rerun scope) would take days in a single shell.
- Each search is single-threaded and reads the same DB → embarrassingly
  parallel across queries. SLURM array is the natural shape.
- Per-task isolation also gets us around the sheetbug `exit(0)` foot-gun in
  one place rather than every caller's loop.

## Resource sizing

Defaults (`array.sbatch`): `--time=01:00:00`, `--mem=2G`, `--cpus-per-task=1`.

- `--time`: the slowest individual query we've measured exceeded 10 min on the
  306 MB DB; an hour is generous but cheap. Bump for a much larger DB.
- `--mem`: the DB is 306 MB; searchmatrix's resident set is well under 1 GB.
  2 GB is a comfortable ceiling. (We saw `bad_alloc` in the 32-bit 2010 binary
  but not in the rebuilt 64-bit one.)
- `--cpus-per-task=1`: searchmatrix has no internal parallelism.

## Known caveats inherited from `searchmatrix`

- `exit(0)` on parse errors and on the sheetbug-fopen-failure path. The
  array-task script ignores `$?` and treats `hit_count` as the only reliable
  signal of "the search did something." Inspect
  `$OUT/logs/<name>.stdout` for "is wrong" / "can't open" lines if a task
  reports `hits=0` and you suspect a translation bug rather than a true miss.
- The query format docs in the top-level README list `sheet S` / `chain S` —
  the parser actually wants `sheetS` / `chainS` as single tokens. Same for
  `sheetD` / `chainD`. The queries under `example/ssp_design_targets/queries/`
  already use the correct form.

## Merging by hand (if the dependent merge job didn't run)

```
{ printf 'query\truntime_sec\texit_code\thits\n'; cat $OUT/parts/*.tsv | sort; } > $OUT/summary.tsv
```
