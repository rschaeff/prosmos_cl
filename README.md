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
cd generateMatrix/src
/usr/bin/mpicxx generMatrix.cpp -o ../build/generateMatrix
```
Builds clean against system OpenMPI 4.1 / g++ 11 on this host (warnings only — uninitialized return paths and `%d` vs `size_t` mismatches; cosmetic). The result is a 64-bit native binary at `generateMatrix/build/generateMatrix`, kept separate from the 2010 `Linux/` binary.

Note: the conda `mpicxx` at `/sw/apps/Anaconda3-2023.09-0/bin/mpicxx` is broken on this host (its `x86_64-conda-linux-gnu-c++` wrapper isn't on PATH) — use `/usr/bin/mpicxx` explicitly.

`searchMatrix` (no MPI):
```
cd searchMatrix/src
g++ -O2 -DSILENT searchMatrix.cpp -o ../build/searchmatrix
```
Builds clean against g++ 11 (warnings only — `%d` vs pointer/`size_t` format mismatches; cosmetic). The result is a 64-bit native binary at `searchMatrix/build/searchmatrix`.

`-DSILENT` redirects stdout to `/dev/null` at startup. The original source has 200+ leftover debug `cout`/`printf` calls inside the per-DB-entry match loop; on the 710k-entry F70 DB that's ~85 MB of stdout per query (measured), which crippled v1/v2 sweeps via NFS I/O contention before we moved logs to compute-node local `/tmp`. The flag eliminates the write entirely without touching the call sites. Errors still surface via `cerr`. Drop the flag if you want the original verbose debug output for a one-off invocation. Note: do NOT enable `-O2` together with this build — it surfaces a pre-existing glibc FORTIFY check inside `searchControl.h` (a `strcpy`/`sprintf` overflow that the original release tolerated under no-opt).

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
searchmatrix <query.txt> <path/to/metamatricesDB> <output_dir>
```
Writes one file per PDB hit into `output_dir`, each listing the matched SSEs (type, position, residue range, chain, length).

### generateMatrix
Three invocation modes:
```
# single PDB (PALSSE .ssd file)
generateMatrix -os <pdbid.ssd> <palsse_dir/> <output_file>

# directory of PALSSE files
generateMatrix -ds <palsse_dir/> <output_file>

# explicit file list
generateMatrix -fs <listfile> <palsse_dir/> <output_file>
```
Concatenate per-PDB outputs into the searchable DB: `cat *.out > metamatricesDB`.

Generation depends on **PALSSE** SSE definitions as input — see the companion working copy at `~/dev/palsse_cl/` (or upstream: http://prodata.swmed.edu/palsse/).

## Query format

Plain text. First line: SSE indices. Second: SSE types (`H` helix, `E` strand, `X` any). Then the upper-triangular interaction matrix, one row per SSE. Followed by optional constraint lines.

Matrix symbols (query side):
| sym | meaning |
|---|---|
| `*` | diagonal |
| `-` | no interaction |
| `X` | don't care |
| `x` | interaction present, angle unchecked |
| `C` | contact, angle < 85° |
| `T` | contact, angle ≥ 95° |
| `N` | contact, 85° ≤ angle < 95° |
| `c` | parallel β-strand pair with H-bonds |
| `t` | antiparallel β-strand pair with H-bonds |
| `u` | contact (no H-bonds), angle < 85° |
| `v` | contact (no H-bonds), angle ≥ 95° |

Constraint lines:
- `handedness <i> <j> <k> R|L` — chirality of an SSE triple
- `length <i> <type> <min> [<max>]` — element-length bounds
- `sheet S|D <i> <j> ...` — same-sheet (S) or not-all-same-sheet (D)
- `chain S|D <i> <j> ...` — same/different chain
- `parallel <i> <j>` / `antiparallel <i> <j>` — orientation of non-H-bonded strands in the same sheet

Example — β-grasp ([`example/query.txt`](./example/query.txt)):
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

## Post-processing (`scripts/`)

- **`fetchmatrix/`** — perl helpers to extract and pretty-print a single PDB's matrix from a hit list, useful when designing queries.
- **`generateinsightIIlog/`** — InsightII visualization log generator.
- **`generatemolscript/`** — MolScript input generator for hit structures.

## References

Shi S, Zhong Y, Majumdar I, Sri Krishna S, Grishin NV. *Searching for three-dimensional secondary structural patterns in proteins with ProSMoS.* Bioinformatics 23(11):1331–8 (2007).

Majumdar I, Krishna SS, Grishin NV. *PALSSE: A program to delineate linear secondary structural elements from protein structures.* BMC Bioinformatics 6:202 (2005).

Contact (upstream): shuoyong.shi@UTsouthwestern.edu
