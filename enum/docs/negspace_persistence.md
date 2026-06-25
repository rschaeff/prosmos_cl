# Persistence of S5 negative space across PDB and AFDB

## Result (definitive, post-fix)

The 800 enumerated S5 negative-space queries (25 zero-hit S5 skeletons
in the v4 manual-rep sweep × 32 H/E typings, staged at
`enum/negspace_queries/`) were swept against:

| DB | Entries | Source | rc=0 | Lit-up queries |
|---|---|---|---|---|
| ecod_rep manual_reps | 18,982 | experimental PDB (curated) | 800/800 | 0/800 |
| ecod_db_pdb_exp | 496,359 | experimental PDB (all) | 800/800 | 0/800 |
| afdb_db | 4,921,931 | AFDB v4 non-singleton clusters | 800/800 | 0/800 |
| **Combined** | **~5.4M** | experimental + predicted | **2,400/2,400** | **0/2,400** |

**All 800 motifs are zero-hit at every scale tested. truly_absent = 800.
AFDB_only = 0. AFDB_loss = 0.** Sweep summary:
`~/work/prosmos_2026/afdb_negspace_final_v2/summary.txt`.

This is the strongest possible empirical absence claim from this
methodology: no motif appears in 5+ million real and predicted protein
structures.

## Methodology corrections (2026-06-24)

Earlier revisions of this document reported partial findings that two
critical fixes invalidated and that this revision supersedes.

### 1. Source-type filter (2026-06-19)

The original F70 v3 build used `derived_files.domain_source_type='pdb'`
as the "PDB" filter, but that column tracks the on-disk file format
(`.pdb`), not experimental provenance. AFDB-predicted domains stored as
`.pdb` files passed through. F70 v3 is actually 74% AFDB-predicted /
26% experimental PDB. The new `pdb_exp` DB uses
`ecod_commons.domain_summary.source_type='pdb'` as the join condition,
which is genuine experimental provenance.

### 2. searchmatrix buffer-overflow fix (2026-06-24, commit d3e1e83)

Both the PDB-experimental and the AFDB sweeps in earlier revisions
crashed with rc=139 (SIGSEGV) or rc=134 (std::out_of_range) and were
incorrectly reported as "zero-hit." The crashes were triggered by
malformed entries in the DB (orphan matrix lines with no preceding
header — `generateMatrix` truncations) that caused:

  (a) A 31-byte strcpy into a 20-byte `pid` buffer in
      `intMnumofele` → 12-byte overflow of adjacent stack memory.
  (b) An unbounded `while(pid[u]!=' ')` walking past the buffer end.
  (c) An unbounded matrix-fill loop in `getInterActionM` writing past
      the heap allocation when a matrix line is longer than the
      header-declared SSE count.

The fix (commit d3e1e83) adds:
  - Header validation in `intMnumofele` (must contain `.ssd`, start
    with digit); returns -1 on parse failure
  - `strncpy` with explicit bound + bounded trim loop
  - Bounds checks on `i` and `j` in the matrix-fill loop
  - Caller-side state recovery: when `intMnumofele` returns -1, reset
    judges and continue without advancing the parity counter

Verified on the 496k-entry experimental DB: 0 of 800 queries crash
post-fix; 1 orphan-block diagnostic emitted; full scan completes
cleanly in ~200s/query (mean).

## Background

`coverage_gaps.md` documented two mechanisms behind the 37% S5 hit rate
of the v4 sweep (7,048 enumerated S3-S5 queries against 19,015 ECOD
`manual_rep=true` experimental-PDB-derived domains):

1. A **composition gradient** — each additional β-strand in the H/E
   typing roughly halves the per-typing coverage (S5 5H/0E = 73%
   → 0H/5E = 7%).
2. **25 of 198 S5 skeletons** are zero-hit across all 32 H/E typings —
   "structurally-realized topologies absent from this set."

The natural follow-up question was whether the 25-skeleton gap is a
sampling artifact (ECOD manual reps is only ~19k of ~500k
experimental-PDB-derived domains, and only a vanishing fraction of all
sequenced protein space) or a real absence. The above result resolves
this: the gap holds at all tested scales, in both real and predicted
protein structure.

## What the result implies

Two robust claims:

1. **The 800 S5 negative-space motifs are absent from observable
   protein structure space.** 5.4M structures spanning curated +
   exhaustive experimental PDB + the AlphaFold-predicted clustered
   subset find zero matches.

2. **AFDB does not extend protein structure space beyond PDB for these
   queries.** Going from 500k experimental + 4.9M predicted produces
   zero new motif realizations. AFDB confirms PDB's negative space
   rather than filling it.

A note on interpretation (a): AFDB v4 was trained on PDB, so AFDB
inherits PDB-like structural priors. The absence in AFDB does not
independently prove non-realizability — it may reflect the predictor's
training-set distribution. The two findings are best read together:
*experimental coverage doesn't have these motifs, and the predictor
trained on that experimental coverage also doesn't predict them.*

## Caveats

- **AFDB singleton clusters are not searched.** The 4.9M subset
  excludes sequences that AFDB clustered as singletons (~210M of
  AFDB's ~214M total). In principle some singletons could realize
  motifs absent from non-singletons. The "absent from AFDB" claim is
  specifically about non-singleton clusters.
- **ProSMoS substructure semantics.** Hits represent embedded SSE-motif
  subgraphs within a domain, not whole-fold matches. A "zero hit"
  means the motif doesn't appear as a substructure of any domain in
  the DB. This is the same semantic as the original
  `coverage_gaps.md` analysis and the same semantic ProSMoS was
  designed for (Medvedev et al. 2021).

## The geometric reachability question

The original framing is unchanged:

> Are these motifs **geometrically possible** (sterically reachable)
> but **not realized** by evolution, or are they **geometrically
> impossible** under the constraints of protein chemistry?

The enum produces them by satisfying:

- 2D hex-lattice planarity (SCC-2 constraint, per Chitturi 2016)
- Compactness filters
- Per-triple handedness signatures (chirality-resolved)

What the enum does not check:

- Steric overlap between SSEs (lattice nodes have zero volume; real
  helices and strands don't)
- Loop-length feasibility for connecting chain segments
- Backbone hydrogen-bonding compatibility for β-sheet motifs (specific
  register requirements smoothed over by the lattice abstraction)
- Statistical preferences from Ramachandran / SSE-packing distributions

Any of these could render an enumerated motif geometrically
unrealizable in 3D, explaining its absence from both PDB and AFDB
without needing to invoke an evolutionary filter.

Two paths to resolve:

1. **Constructive geometry test**: for each of the 800 motifs, attempt
   to build a 3D backbone model satisfying the topology + handedness,
   with realistic SSE geometry + loop lengths. Use RFdiffusion or
   ProteinMPNN-style backbone-design tooling, or a constraint-solver
   approach.
2. **Negative geometry test**: derive analytic geometric constraints
   from the lattice → 3D map. A tighter enum that requires each
   adjacency to correspond to a realizable SSE-SSE contact distance +
   angular range, and each non-adjacency to a non-clashing separation,
   would pre-filter the 198-skeleton set to a geometrically-valid
   subset.

If most of the 800 turn out geometrically valid by approach 1, the
result is a strong claim: **evolution has not produced a substantial
fraction of designable protein topologies**. If most are
geometrically impossible, the lattice enum is over-generating relative
to 3D space.

## Status

| Component | State |
|---|---|
| v4 manual-rep sweep | done; `results/v4_ecod_manual_reps_summary.tsv` |
| pdb_exp_negspace_sweep_v2 (496k DB, post-fix) | done; 800/800 zero-hit |
| afdb_negspace_sweep_v3 (4.9M DB, post-fix) | done; 800/800 zero-hit |
| Final cross-reference | done; `~/work/prosmos_2026/afdb_negspace_final_v2/` |
| AFDB singleton sweep | not run — would close the singleton caveat |
| Constructive geometry test | not run — the next research direction |
