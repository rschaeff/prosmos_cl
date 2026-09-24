# ProSMoS

**Pro**tein **S**econdary structure **Mo**tif **S**earch — search the PDB for proteins whose secondary-structure topology matches a user-defined motif. Originally released ~2010 by Shuoyong Shi (UT Southwestern); free for academic use.

This tree is a hardened working copy of the upstream source; the original release notes are
preserved as [`readme.original`](./readme.original).

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22075370.svg)](https://doi.org/10.5281/zenodo.22075370)

---

## Citing this work

This repository provides the **search engine** used by:

> Schaeffer RD, Guo R, Cong Q, Grishin NV. *A topology census of four-strand β-sheets in
> experimental and predicted protein structures.* [TODO: journal, year, DOI]

The search engine itself is archived on Zenodo. Cite the **version** DOI in a paper's
methods; the **concept** DOI below always resolves to the latest version:

> Schaeffer RD. *ProSMoS — protein secondary-structure motif search engine.* Zenodo.
> https://doi.org/10.5281/zenodo.22075370

The **analysis** for that paper — the enumeration, both censuses, every control and
sensitivity analysis, the tables and the figure renderers — is deposited separately, with a
per-file manifest and checksums:

> [TODO: Zenodo DOI]

Please cite both if you use the pipeline end to end. The original method is:

> Shi S, Zhong Y, Majumdar I, Sri Krishna S, Grishin NV. Searching for three-dimensional
> secondary structural patterns in proteins with ProSMoS. *Bioinformatics* 2007;23(11):1331–8.
> Shi S, Chitturi B, Grishin NV. ProSMoS server. *Nucleic Acids Res* 2009;37:W526–31.

### What is and is not part of the paper

| | |
|---|---|
| `searchMatrix/`, `generateMatrix/` | the pipeline the census ran on — **in scope** |
| `enum/docs/positive_controls.md` | recall controls quoted in the paper's SI — **in scope** |
| `archive/` | exploratory work, explicitly **not** part of the census; see [`archive/README.md`](./archive/README.md) |
| `enum/` (rest), `scripts/` | working code and session notes, not paper artifacts |

The paper's scope is deliberately narrow: ProSMoS is used as an implementation of an
already-defined 96-state description of the four-strand β-sheet. Neither the paper nor this
repository claims ProSMoS is an adequate general representation of protein fold space.

### Corrections in this tree

Running the legacy code at whole-database scale required fixing defects that are load-bearing
for the census. The most consequential is `1e5546c`: `intMnumofele()` bounded its
element-reading loop by input line length rather than by the declared element count, so
structures with five-digit coordinates overflowed an eight-character field, injected phantom
SSEs, and could abort a chunk — while the batch wrapper swallowed the non-zero exit. Also
here: a `Makefile` that refuses to sweep against a stale binary (`fbd3d0f`), buffer-overflow
fixes in the query-path and output-path handling, and a hardened database reader.

## What it does

ProSMoS reduces each protein structure to a square interaction matrix over its secondary-structure elements (SSEs). The matrix captures, for every pair of SSEs, whether they interact and how — parallel/antiparallel β-strand H-bonding, helix-helix contact angle, sheet/chain membership, handedness. A user-supplied query matrix (a small motif) is then matched against a precomputed database of these matrices over the PDB.

```
PDB + PALSSE SSE definitions
        │
        ▼  generateMatrix  (C++ / MPI)
  per-structure N×N interaction matrix
        │
        ▼  cat *.out > metamatricesDB
  metamatricesDB (~306 MB)
        │
   query.txt ──► searchMatrix (C++)
        │
        ▼
   ranked hits  (PDB id + matched SSE ranges)
```

The two-stage design (expensive offline matrix generation, fast online search) is what makes whole-PDB motif scans tractable.

## Layout

```
generateMatrix/         offline: build interaction matrices from PALSSE output
  src/                  C++ + MPI sources, original Makefile
  Linux/generateMatrix  32-bit i386 binary (2010)  -- needs libstdc++.so.5
  build/generateMatrix  64-bit binary rebuilt from src (OpenMPI 4.1, g++ 11)
searchMatrix/           online: motif search against metamatricesDB
  src/                  C++ sources
  Linux/searchmatrix    32-bit i386 binary (2010)  -- loads but crashes (bad_alloc) on the real DB
  build/searchmatrix    64-bit binary rebuilt from src (g++ 11)
scripts/
  acceptance/           post-search acceptance rule for the four-strand census
  fetchmatrix/          perl helpers for inspecting individual matrices
  generateinsightIIlog/ build InsightII visualization log from hits
  generatemolscript/    build MolScript input from PDB hits
example/
  query.txt             β-grasp motif query
  pdbhitslist           sample hit list
readme.original         upstream readme (full original docs)
```

What was deliberately **not** copied from the upstream tree:
- `metamatrixdb/metamatricesDB` (306 MB precomputed DB) — at `/home/rschaeff/src/Prosmos/ProSMoS/metamatrixdb/`
- `generateMatrix/Linux/*.pmm` (60k+ per-PDB matrix files, ~249 MB)
- `scripts/fetchmatrix/transgenmatrix/` and large `*.tar.gz` / sample `*.list` files
- `example/testlog.tar.gz` (121 MB), `pdbhits.tar`, `super1.71.tar`
- SunOS binary (legacy)

If you need any of those, pull them directly from `/home/rschaeff/src/Prosmos/ProSMoS/`.

## Building from source

`generateMatrix` (MPI, needed only if you regenerate the DB):
```
mkdir -p generateMatrix/build      # not tracked; a fresh clone lacks it
/usr/bin/mpicxx generateMatrix/src/generMatrix.cpp -o generateMatrix/build/generateMatrix
```
Requires system OpenMPI (`/usr/bin/mpicxx`, the `libopenmpi-dev` package) and g++.
Builds clean against system OpenMPI 4.1 / g++ 11 on this host (warnings only — uninitialized return paths and `%d` vs `size_t` mismatches; cosmetic). The result is a 64-bit native binary at `generateMatrix/build/generateMatrix`, kept separate from the 2010 `Linux/` binary.

Note: the conda `mpicxx` at `/sw/apps/Anaconda3-2023.09-0/bin/mpicxx` is broken on this host (its `x86_64-conda-linux-gnu-c++` wrapper isn't on PATH) — use `/usr/bin/mpicxx` explicitly. The binary links system OpenMPI, so never launch it with conda's `mpirun` either (see the invocation notes: you should not need `mpirun` at all).

`searchMatrix` (no MPI):
```
make -C searchMatrix          # g++ -O2 -DSILENT; `make stale` exits 1 if build/ is older than src/
```
Builds clean against g++ 11 (warnings only — `%d` vs pointer/`size_t` format mismatches; cosmetic). The result is a 64-bit native binary at `searchMatrix/build/searchmatrix`.

`-DSILENT` redirects stdout to `/dev/null` at startup. The original source has 200+ leftover debug `cout`/`printf` calls inside the per-DB-entry match loop; on the 710k-entry F70 DB that's ~85 MB of stdout per query (measured), which crippled v1/v2 sweeps via NFS I/O contention before we moved logs to compute-node local `/tmp`. The flag eliminates the write entirely without touching the call sites. Errors still surface via `cerr`. Drop the flag if you want the original verbose debug output for a one-off invocation. Use the Makefile rather than a hand-written command: it builds at `-O2`, which is about 4x faster than `-O0`, and hit sets are byte-identical at every optimisation level.

The source as released does not compile on modern g++ until one duplicate parameter name is fixed: `searchControl.h:30` declared `searchM(...)` with two parameters both named `a`, which older g++ tolerated but g++ ≥ ~6 rejects as a conflicting declaration. The second `a` (the `vector<matrixElment>&` one) has been renamed to `totalele` to match the existing definition at line 1628.

The upstream Makefile in `generateMatrix/src/` only does `g++ -c` and references `/usr1/local/include`; treat it as a hint, not a working build script.

### The shipped Linux binaries

Both prebuilt binaries are **32-bit i386 ELF** from 2010.

| binary | runs here? | why |
|---|---|---|
| `searchMatrix/Linux/searchmatrix` | no | loads (32-bit deps resolve via `/lib32/`) but aborts with `std::bad_alloc` on the real DB — use `searchMatrix/build/searchmatrix` instead |
| `searchMatrix/build/searchmatrix` | yes | rebuilt from `src/` against g++ 11 after the duplicate-parameter fix (see above); end-to-end smoke test against the 306 MB metamatricesDB produces hits |
| `generateMatrix/Linux/generateMatrix` | no | links `libstdc++.so.5` — not installed |
| `generateMatrix/build/generateMatrix` | yes | rebuilt from `src/` against system OpenMPI |

To run the prebuilt `generateMatrix`, either use `generateMatrix/build/generateMatrix` (rebuilt from source here), install the legacy compat lib for the 2010 binary (Debian/Ubuntu: `libstdc++5:i386`), or rebuild yourself as below.

## Usage

### searchMatrix
```
searchmatrix <query.query> <metamatricesDB.clean> <output_dir>      # one query
searchmatrix <manifest>    <metamatricesDB.clean> <output_dir>      # many, one DB pass
```
Writes one file per matching DB record, each listing the matched elements (type, position,
residue range, chain, length). With a single query the hit files go straight into
`<output_dir>`. A manifest lists query files, one path per line, and each query's hits go into
`<output_dir>/<query name>/`; the DB is parsed once for all of them, which is what the SLURM
sweeps (`scripts/slurm_search/`) run. The output directory is created if missing, and its
trailing slash is optional. `searchmatrix --help` prints the summary.

Always search the `.clean` DB written by `scripts/db_validate.py --clean`: a few malformed
`generateMatrix` records can desynchronise the reader.

Exit status: 0 on success; 1 for an unreadable query, manifest entry or DB, an unwritable
output directory, or a malformed query or DB record (the message is on stderr); 2 for bad
arguments. Until 2026-09-22 every one of these exited 0: the production build sends stdout to
`/dev/null` and the errors went to stdout, so a typo in a path produced a clean-looking run with
no hits. The same silent zero happened when the working directory's parent was not writable
(the engine appends debug output to `../sheetbug/`, which is now optional). The search itself is
unchanged: `scripts/db_validation/searchmatrix_args_regression.sh` compares hit trees against
the previous build byte for byte.

### generateMatrix
One PALSSE file per call, no `mpirun`:
```
generateMatrix <file.ssd> <output_file>                  # -os is the default
generateMatrix -os <name.ssd> <palsse_dir> <output_file>  # original form, still accepted

for f in <palsse_dir>/*.ssd; do
  b=$(basename "$f" .ssd)
  generateMatrix "$f" <out_dir>/"$b".out
done
cat <out_dir>/*.out > metamatricesDB
```
`-os` (the default) takes sheets from the PALSSE file's SHEET records, which is how the
census DBs were built; `-o` recomputes them (legacy). The record is named after the
`.ssd` file's basename. The trailing slash on `<palsse_dir>` is optional. Run
`generateMatrix --help` for the summary.

It runs as an OpenMPI singleton, so it needs no `mpirun`. Parallelise with `xargs -P` or a
SLURM array, as `scripts/*_db_build/process_chunk.sh` does.

Errors exit nonzero with a message on stderr: an unreadable input (checked before the
output is created, so no empty file is left behind), an unwritable output, or bad
arguments (2). A structure with no helix or strand elements still writes its (empty)
record, with a warning.

**The batch modes `-ds <palsse_dir/>` and `-fs <listfile>` are disabled** (exit 2). They
never worked in this build: run directly they hang (rank 0 is manager-only and waits for
workers that do not exist), and under `/usr/bin/mpirun -np N` they exit 0 but write every
record with **zero elements**, because the worker path never calls `prepareIndex()`.

Until 2026-09-22 the front end failed silently: a directory without the trailing slash was
concatenated straight onto the file name, and a missing input exited 0 after creating an
empty output. The records themselves are unchanged. `scripts/db_validation/genmat_args_regression.sh`
checks the new front end against the v1.0.0 binary byte for byte, in every accepted spelling.

Generation depends on **PALSSE** SSE definitions as input — see the companion working copy at `~/dev/palsse_cl/` (or upstream: http://prodata.swmed.edu/palsse/).

## Query format

Plain text. First line: element numbers. Second: element types (`E` strand, `H` helix,
`X` any). Then the upper triangle of the query matrix, one row per element, `*` on the
diagonal. Then optional constraint lines. The smallest useful query, two strands paired
antiparallel in one sheet:
```
1 2
E E
* t
  *
sheetS 1 2
length 1 E 5 1000
length 2 E 5 1000
```
For a parallel pair change `t` to `c`. No single code means "paired, either direction",
so the two orientations are two queries.

### Matrix codes

`generateMatrix` writes six codes into the database. A query may use those six, which
match only themselves, or four wildcards (`searchControl.h`, `notequal()`):

| query code | matches DB code | meaning |
|---|---|---|
| `t` | `t` | strands paired antiparallel (see below) |
| `c` | `c` | strands paired parallel |
| `u` | `u` | in contact, axes at < 85° (pointing the same way) |
| `v` | `v` | in contact, axes at ≥ 95° (pointing opposite ways) |
| `N` | `N` | in contact, 85°–95° (roughly perpendicular) |
| `-` | `-` | no contact (but see below for strands in one sheet) |
| `T` | `v` or `t` | in contact, opposite-pointing, paired or not |
| `C` | `u` or `c` | in contact, same-pointing, paired or not |
| `x` | anything but `-` | in contact, angle not checked |
| `X` | anything | not checked |

"In contact" means the element axes overlap by more than 2.5 Å and lie within 11 Å
(`OVERLAP_DEFAULT`, `DISTANCE_DEFAULT` in `generateMatrix/src/external.h`).

Which codes actually occur depends on the element types. Counted over every element pair
in a 5,000-domain AFDB sample:

| pair | codes seen (most to least common) |
|---|---|
| E–E | `-` `t` `v` `u` `c` `N` |
| H–E | `-` `v` `u` `N`; `t`/`c` in 230 of ~42,000 pairs |
| H–H | `-` `v` `u` `N`; `t` in 6 of ~46,800 pairs |

So on a helix–strand pair `T` behaves as `v` and `C` as `u`.

### What `t` and `c` mean: PALSSE pairing, not hydrogen bonds

The upstream readme describes `t`/`c` as "more than 2 H-bonds". Nothing in the pipeline
counts hydrogen bonds:

1. **PALSSE pairs residues from Cα positions only.** It scores quadruplets (two residues on
   each of two strands) against empirical tables of Cα–Cα distances, virtual-bond geometry
   and torsions, chains the accepted ones into ladders, and writes every paired residue as a
   `PAIRS` record in the `.ssd`. Its `bond_score` is the virtual Cα–Cα bond. No N–H···O
   geometry or energy is used.
2. **`generateMatrix` counts those pairs** (`h_bond_E()` in `control.h`, despite the name).
   For elements A and B it counts A's PALSSE partners that fall inside B. **Two or more**
   makes the pair `t` or `c`: `t` if the first and last partners run backwards along B,
   otherwise `c`.

Consequences:

- A tie (every counted partner is the same residue of B) is coded `c`, parallel.
- Two strands that are in contact and in the same PALSSE sheet but not directly paired get
  `-`, not an angle code. Between strands of one sheet, `-` means "not ladder neighbours".
- The pairing count also runs for elements that are *not* in contact, and for any element
  types. That is where the rare helix `t`/`c` codes come from.
- A `t`/`c` pair is not always in one PALSSE sheet: in the AFDB sample 237 of 7,731 `t` and
  685 of 2,573 `c` pairs are not. Add `sheetS` when you mean one sheet.

If a result depends on real hydrogen bonds, check it against DSSP, which assigns bridges
from backbone H-bond energies. The four-strand census does this: every `t`/`c` pair in a
hit must carry a DSSP bridge.

### Constraint lines

- `length <i> <type> <min> [<max>]`: element-length bounds, inclusive. Give a generous
  maximum (`1000`); a low one silently drops long elements.
- `sheetS <i> <j> ...`: all in one PALSSE sheet. `sheetD <i> <j> ...`: not all in one sheet.
- `chainS` / `chainD`: same or different chain, likewise.
- `parallel <i> <j>` / `antiparallel <i> <j>`: orientation of strands in one sheet that are
  not paired with each other.
- `handedness <i> <j> <k> R|L`: chirality of a three-element unit. It is not reliable
  enough to support a claim of left-handedness.

`sheetS`, `sheetD`, `chainS` and `chainD` are single words. The upstream readme writes
`sheet S`, which the parser does not accept.

Example: β-grasp ([`example/query.txt`](./example/query.txt)):
```
1 2 3 4 5
E E H E E
* t C - c
  * T - -
    * x C
      * t
        *
handedness 2 3 4 R
length 3 H 8 1000
```

### Reading database records

`scripts/metamatrix_record.py` parses a `generateMatrix` output file or a whole
metamatricesDB into elements, sheets and the full N×N matrix:
```
scripts/metamatrix_record.py metamatricesDB.clean <name>       # print records
python3 -c "from metamatrix_record import iter_records"         # or import it
```
The record header is fixed width. Slice it; never split it on whitespace or regex it.

## Post-processing (`scripts/`)

- **`acceptance/`** — the four-strand census's acceptance rule, run on `searchmatrix` hits: disjointness, axial overlap, local pairing and DSSP bridges, as sequential stages that record why each candidate was rejected or could not be judged. See `scripts/acceptance/README.md`.
- **`fetchmatrix/`** — perl helpers to extract and pretty-print a single PDB's matrix from a hit list, useful when designing queries.
- **`generateinsightIIlog/`** — InsightII visualization log generator.
- **`generatemolscript/`** — MolScript input generator for hit structures.

## References

Shi S, Zhong Y, Majumdar I, Sri Krishna S, Grishin NV. *Searching for three-dimensional secondary structural patterns in proteins with ProSMoS.* Bioinformatics 23(11):1331–8 (2007).

Majumdar I, Krishna SS, Grishin NV. *PALSSE: A program to delineate linear secondary structural elements from protein structures.* BMC Bioinformatics 6:202 (2005).

Contact (upstream): shuoyong.shi@UTsouthwestern.edu
