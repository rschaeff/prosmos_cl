# Session log: AFDB negspace sweep + searchmatrix patch (2026-06)

Working log for the multi-week session that produced the
`negspace_persistence` result and the searchmatrix buffer-overflow fix.
Organized as a timeline; each phase links to its commits and the
artifacts it produced.

## Goal

Test the v4-identified S5 negspace claim (25 of 198 S5 skeletons return
zero hits against 19k ECOD manual reps) at broader scale: against the
full experimental PDB (500k) and against AFDB (4.9M predicted). The
question driving everything: are these motifs absent from observable
protein structure, or only from the curated subset?

## Phase 1: AFDB infrastructure (start of session — d7b788a region)

- **Handoff bundle** for Gray at `/data/ECOD0/html/distributions/prosmos_afdb_negspace/`
  containing prosmos_cl + palsse_cl tarballs + RUNBOOK + the 800 negspace
  queries — anticipated AFDB sweep would happen on HGD.
- Direction changed mid-session: Gray moved the AFDB non-singleton
  subset (~4.9M structures) to leda at `~grey/afdb.200m/non_singleton_4p9m_structures/`,
  so the sweep ran here.

## Phase 2: AFDB metamatricesDB build (commits f09e5e5, eae5e57)

- Built `scripts/afdb_db_build/{submit.sh,process_chunk.sh,array.sbatch}`
  adapted from `ecod_db_build/` to handle AFDB's tarball-shard layout
  (`new/<shard>/<XX>.tar.gz` etc.) — extract-then-PALSSE-then-generateMatrix
  per chunk.
- First build (job 556195, 210 chunks × 506 tarballs): 134 of 210 chunks
  TIMED OUT at the 2h limit (PALSSE hangs accumulating). DB built but
  only 22% complete (~1.26M of expected ~5.7M entries).
- Recovery (job 557199): re-chunked at 100 work units each (679 chunks)
  with longer timeout headroom. 0 timeouts this run. Combined with the
  76 surviving original fragments → 4.92M-entry DB.
- Filter dropped 5 malformed entries (`>4-digit-coord` pattern from the
  earlier v3 build, same awk filter) → metamatricesDB.clean.

## Phase 3: First AFDB sweep "result" — false zero (commits c167c26)

- Auto-fired the 800-query negspace sweep against the 4.9M AFDB DB.
- Sanity checks passed: FHB returned 12,183 hits, RLM 1,193 hits.
- All 800 negspace queries returned 0 hits — initial framing in
  `negspace_persistence.md` claimed the AFDB-only quadrant was empty.
- This claim was **wrong**. Revisited later (see Phase 6).

## Phase 4: Source-type filter discovery (commit c167c26)

- User asked: are these PDB hits actually experimental crystal
  structures? Discovered that the F70 v3 manifest filter
  `derived_files.domain_source_type='pdb'` selects on **file format**
  (any `.pdb`-suffixed file), not experimental provenance. AFDB-predicted
  structures stored as `.pdb` files passed through.
- Actual F70 v3 composition: **74% AFDB-predicted, 26% experimental PDB**.
  The "PDB scale-up" framing of sweep A in earlier writeups was
  meaningless — both DBs being compared were AFDB-heavy.
- Manual-reps DB is unaffected (`ecod_rep.domain WHERE manual_rep` is
  100% experimental by curation).
- Fix: pull the correct filter via `ecod_commons.domain_summary.source_type='pdb'`
  → 500,313 truly experimental PDB-derived entries.

## Phase 5: Experimental-PDB DB build + crash discovery (commit d3e1e83)

- Built `ecod_db_pdb_exp/metamatricesDB.clean` (496,359 entries after
  filtering 439 malformed).
- Sanity check passed: FHB 22,627 hits, RLM 21,263 hits — confirmed the
  DB is well-populated with real Rossmann/four-helix structures (far
  more per-entry than AFDB).
- Negspace sweep: **all 800 queries crashed**. Initial reading was
  rc=139 (SIGSEGV) but bisection revealed the actual exit codes split
  rc=134 (SIGABRT, `std::out_of_range`) and rc=139 (SIGSEGV) depending
  on which malformed entry surfaced first.
- Binary-search bad-entry hunt (jobs 561358, 561408, 561446): bisected
  the DB into chunks, identified that 61–68% of sub-chunks crashed at
  every level. Bad entries were **dense**, not localized.

## Phase 6: Root cause + searchmatrix patch (commit d3e1e83)

Diagnosis: bisecting one chunk pinpointed entry 000154400 as the
crash trigger, but that entry **alone** worked fine. The crash was
state-accumulation: searchmatrix's `intMnumofele` and `getInterActionM`
mishandled DB transitions when an orphan block (matrix line without
preceding header) appeared.

Two bugs in `searchControl.h`:

1. **`intMnumofele` line ~1203**: `strcpy(pid, temp.c_str())` writes
   31 bytes into a 20-byte `pid` buffer. For normal headers the 12
   extra bytes are padding spaces (benign). For orphan matrix lines
   (no spaces in first 31 chars), 12 bytes of matrix data overflow
   adjacent stack — AND `while(pid[u]!=' ') u++` has no upper bound,
   walking past the buffer until finding a stray space byte
   somewhere in memory.
2. **`getInterActionM` line ~1135**: matrix-fill loop writes
   `intM[i][j]` without bounds-checking `i` or `j`. A matrix line
   longer than the header's declared `mrow` SSEs (which happens
   during misclassified-line cascades) overflows the heap allocation.

Root structural cause: the DB has malformed blocks (`HEADER-no-matrix`
and orphan matrix lines, ~57 such blocks scattered through 496k
entries) from rare `generateMatrix` output truncations. The original
`searchControl.h` line-type tracking is parity-based
(`readLincon % 2 == 0` → header, else matrix); orphan blocks throw
the parity off, leading to header-parsing code being run on matrix
lines.

Fix (commit d3e1e83, ~37 lines):
- Validate header in `intMnumofele` (must contain `.ssd`, start with
  digit). Return `-1` on parse failure.
- `strncpy` with explicit bound + bounded trim loop.
- Bounds-check `i, j` in `getInterActionM`.
- Caller in `oneprocess`: on `intMnumofele` returning -1, reset judges
  and `continue` without advancing `readLincon` — recovers parity.

Verified: 800-query sweep against 496k DB completes cleanly in ~200s
per query, 1 orphan-block diagnostic emitted, 0 crashes.

## Phase 7: Re-sweep, definitive result (commit f44a17e)

After the binary fix:

| DB | Entries | Source | 800 queries: rc=0 | Lit-up |
|---|---|---|---|---|
| ecod_rep manreps (v4) | 18,982 | experimental PDB (curated) | 800/800 | 0/800 |
| ecod_db_pdb_exp | 496,359 | experimental PDB (all) | 800/800 | 0/800 |
| afdb_db (v3 sweep) | 4,921,931 | AFDB v4 non-singleton | 800/800 | 0/800 |

**5.4M structures searched across experimental + predicted protein
space; 800/800 negspace queries return zero hits in every database.**

`truly_absent = 800, AFDB_only = 0, AFDB_loss = 0, common = 0`.

Earlier "AFDB-only = 0" finding in `negspace_persistence.md`
(c167c26) was based on 799/800 crashed queries, so superseded by
this clean result. Writeup amended in f44a17e to reflect the post-fix
data and document the bug-fix journey.

## Phase 8: Paper section + figures (commits 1357ccc, 655bf62)

- `enum/docs/paper_section.md` — Methods + Results draft (4 paragraphs
  each) for inclusion in the 200M paper.
- `enum/docs/figures/s5_hit_grid.png` — heatmap of all 6,336 S5
  queries × hit count, sorted by total hits, 25-row negspace block
  marked. Generated by `enum/scripts/plot_s5_hit_grid.py`.
- `enum/docs/figures/s5_skeleton_schematic.png` — 2 × 3 panel showing
  one hit-rich skeleton (s5-0098, RLM-canonical) and one zero-hit
  skeleton (s5-0001) side-by-side: lattice topology → typing applied →
  ProSMoS interaction matrix. Generated by
  `enum/scripts/plot_skeleton_schematic.py`.

## Phase 9 (in flight): three-DB comparable hit grids + spotlight

Planning doc embedded in session conversation; key elements:

- **Three-panel hit grid figure**: same 6,336 S5 queries against
  manreps (19k) | PDB-exp (496k) | AFDB (4.9M). All three sorted by
  manreps total hits so cells are positionally comparable. Reveals
  how the hit distribution shifts across DB scale.
- **Low-hit skeleton spotlight**: 5–10 skeletons that sit just above
  the negspace boundary (1–3 typings hit out of 32). For each, list
  the actual ECOD domains that lit each typing, render those
  structures in PyMOL to show what the rare topology realization
  looks like in nature.

Sweeps for the three-panel figure:
- PDB-exp full S5 (6,336 queries): jobs 563583–563590, ~5.6h wall
- AFDB full S5: 563591–563593 + background recovery launcher submitting
  the remaining 4 chunks + merge as SLURM rate-limit allows.
  ~13h wall total.

## Commits this session (chronological)

```
d7b788a  slurm_search: watch_failures.sh — flag nonzero-rc parts
80700c1  searchMatrix: -DSILENT mutes 200+ debug prints inside the match loop
821db38  slurm_search: chunked arrays, local-/tmp logs, recursive query manifest
8ab795e  enum: assignment.py — Skeleton -> ProSMoS query pipeline, S3-S5 corpus
1654def  results: v4 ECOD manual-reps sweep summary (7048 queries)
534ed29  enum: coverage_gaps.md — S3-S5 enum vs ECOD manual reps writeup
2cb77c6  enum: positive_controls.md — FHB + RLM recovery validation
f09e5e5  slurm_search: recover_s6.sh — throttled recovery for large sweeps
eae5e57  handoff: AFDB negative-space sweep — runbook + 800 queries
56117c8  enum: negspace_persistence.md — PDB+AFDB cross-DB findings + open Qs
ccf3e42  cleanup: consolidate ~/work/ paths under prosmos_2026/
c167c26  enum: amend negspace_persistence.md — correct source-type framing
d3e1e83  searchMatrix: fix pid buffer overflow and matrix-fill OOB write
f44a17e  negspace: definitive post-fix result + new quadrant analyzer
1357ccc  enum: paper_section.md + s5_hit_grid figure
655bf62  enum: skeleton schematic — lattice -> typing -> ProSMoS query
```

## Durable findings worth carrying forward

- The searchmatrix binary on this branch (d3e1e83) is the canonical
  one; the previous binary silently crashes on any DB containing
  malformed `generateMatrix` outputs (which means almost any newly-
  built DB).
- `derived_files.domain_source_type` ≠ experimental provenance. To
  filter for actual crystal structures, join against
  `ecod_commons.domain_summary.source_type='pdb'`.
- The 800-query negspace set is now confirmed against 5.4M
  structures. Any further validation should target the geometric
  reachability question (constructive 3D-backbone modeling), not
  another scale-up sweep.
- Phase 9's three-panel figure will be the final headline visual
  for the paper; the spotlight subset will be the structural
  ground-truth examples.
