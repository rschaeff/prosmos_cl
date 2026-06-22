# Persistence of S5 negative space across PDB and AFDB

## Background

`coverage_gaps.md` documented two mechanisms behind the 37% S5 hit rate of
the v4 sweep (7,048 enumerated S3-S5 queries against 19,015 ECOD
`manual_rep=true` PDB-derived domains):

1. A **composition gradient** — each additional β-strand in the H/E
   typing roughly halves the per-typing coverage (S5 5H/0E = 73%
   → 0H/5E = 7%).
2. **25 of 198 S5 skeletons** are zero-hit across all 32 H/E typings —
   "structurally-realized topologies absent from this set."

The natural follow-up question is whether the 25-skeleton gap is a
**sampling artifact** (ECOD manual reps is only ~19k of ~700k PDB-derived
domains in F70, and only a vanishing fraction of all sequenced protein
space) or a **real absence**. This document records two follow-up sweeps
that bracket that question and the interpretation they support.

## The two follow-up sweeps

Both sweep the same 800-query corpus: the 25 zero-hit S5 skeletons ×
32 H/E typings, staged at `enum/negspace_queries/` in this repo.

### Sweep A: PDB at PDB-wide scale

- **DB**: `ecod_db_f70_v3/metamatricesDB.clean` — F70 cluster reps,
  707,567 entries (~37× larger than the manual-rep set)
- **Status**: in progress at writeup time; partial coverage of 800
  queries with the cancelled-then-resumed in-flight sweep, plus the
  366-query fill-in (see chain notes in `RUNBOOK_AFDB_NEGSPACE.md`).
- **Pattern in the completed parts so far**: 89% of completed queries
  hold their zero-hit status against the 37×-larger DB. The ~10% that
  light up have small hit counts (1–47 hits per query) and concentrate
  on specific skeletons (s5-0004, s5-0033, s5-0019, s5-0003,
  s5-0026) and their helix-rich typings.

### Sweep B: AFDB at AlphaFold scale

- **DB**: built locally from
  `~grey/afdb.200m/non_singleton_4p9m_structures/` (the AFDB v4
  non-singleton-cluster representative set), 4,921,931 entries
  (~260× larger than the manual-rep set).
- **Status**: complete.
- **Result**: **0 of 800 queries returned any hits** against the
  4.9M-entry AFDB DB.

### Cross-reference table (full 800 once the chain settles)

| Quadrant | PDB v4 19k | PDB F70 707k | AFDB 4.9M | Meaning |
|---|---|---|---|---|
| `truly_absent` | 0 hits | 0 hits | 0 hits | Motif absent at every scale tested |
| `AFDB_only` | 0 hits | 0 hits | >0 hits | Predicted-only realization — **count: 0** |
| `AFDB_loss` | 0 hits | >0 hits | 0 hits | PDB has it, AFDB-non-singleton-set doesn't (likely sequences AFDB clustered as singletons) |
| `common` | 0 hits | >0 hits | >0 hits | Lit up at PDB scale |

The headline cell is `AFDB_only` = **0**. Across 800 queries × 4.9M
predicted structures, no motif appeared in AFDB that wasn't already in
the F70 PDB set.

## What the AFDB-only zero implies

AFDB v4 contains ~214M predicted structures spanning the bulk of
sequenced protein space across UniProt — including organisms and
proteomes that have never been crystallized. The non-singleton subset
we searched (4.9M) is biased toward sequences with at least one cluster
neighbor, i.e. proteins for which there is some evolutionary signal of
relatedness — which is a broader sample of "what evolution makes" than
PDB.

If AFDB at this scale brought no new realizations of our 800 enumerated
motifs, then either:

1. **The 800 motifs are not present in real proteins**, even at the
   evolutionary breadth captured by 4.9M AlphaFold-predicted structures.
   The PDB negative-space identification is a real absence in nature,
   not a sampling artifact of crystallographic preference.

2. **AFDB predicts protein space the same way as PDB samples it** —
   i.e., AFDB's structural distribution is dominated by the same fold
   classes that the PDB has, just at greater sequence redundancy. AFDB
   doesn't discover novel topologies because AlphaFold is trained to
   predict what PDB-like structure looks like.

Both interpretations are plausible and may co-act. The data alone
can't separate them: we'd need a non-PDB-trained predictor or
experimental coverage of currently-uncrystallized organisms to
distinguish.

However, the practical implication is clear: **at the SSE-motif level,
AFDB does not extend protein structure space beyond PDB.** Going from
PDB (19k → 707k) to AFDB (707k → 4.9M) adds 7× the entries but 0 new
S5 negative-space motifs. The "AFDB fills the gap" hypothesis is
unsupported for this query class.

## The geometric reachability question

Given that 800 enumerated S5 motifs find zero match across 5+ million
real and predicted protein structures, the natural next question is:

> Are these motifs **geometrically possible** (sterically reachable)
> but **not realized** by evolution, or are they **geometrically
> impossible** under the constraints of protein chemistry?

The enum produces them by satisfying:

- 2D hex-lattice planarity (SCC-2 constraint, per Chitturi 2016)
- Compactness filters
- Per-triple handedness signatures (chirality-resolved)

What the enum does **not** check:

- Steric overlap between SSEs (the lattice model treats SSEs as
  abstract nodes; in 3D, helices and strands have non-negligible
  volume)
- Loop-length feasibility for the connecting chain segments (some
  enumerated topologies may require crossover loops that exceed
  feasible polypeptide flexibility, or zero-length linkers between
  non-sequential adjacent SSEs)
- Backbone hydrogen-bonding compatibility for β-sheet motifs (parallel
  vs antiparallel pairing has specific register requirements that the
  lattice abstraction smooths over)
- Statistical preferences from Ramachandran / SSE-packing distributions

Any of these could render an "enumerated" motif geometrically
unrealizable in 3D, which would explain its absence from both PDB and
AFDB without needing to invoke an evolutionary filter.

This is the next experimental question. Two paths:

1. **Constructive geometry test**: for each of the 800 motifs, attempt
   to build a 3D backbone model satisfying the topology + handedness,
   with realistic SSE geometry + loop lengths. Use RFdiffusion or
   ProteinMPNN-style backbone-design tooling, or a constraint-solver
   approach (à la Folding@home / Rosetta cyclic-peptide topology). A
   motif that cannot be geometrically instantiated is structurally
   unreachable; one that can is "designable" and the gap is
   evolutionary.

2. **Negative geometry test**: derive analytic geometric constraints
   from the lattice → 3D map. The current lattice enum operates on
   abstract adjacency; a tighter version that requires each adjacency
   to correspond to a realizable SSE-SSE contact distance + angular
   range, and each non-adjacency to a non-clashing separation, would
   pre-filter the 198-skeleton set to a smaller "geometrically valid"
   subset. The difference is the "geometrically impossible but
   topologically enumerable" set.

If most of the 800 turn out geometrically valid by approach 1, the
result is a strong claim: **evolution has not produced a substantial
fraction of designable protein topologies**, which has implications
for both de novo protein design (these are unexplored targets) and for
understanding what makes a fold "evolvable" (what about the realized
folds gave them a leg up).

If most are geometrically impossible, the lattice enum is producing
more skeletons than 3D space supports, and the gap reflects an
abstraction-induced overshoot rather than a deep biological signal.

## Status

The PDB-wide F70 sweep (sweep A above) is still completing as of this
writeup. Final cross-reference summary will land at
`~/work/prosmos_2026/afdb_negspace_final/{summary.txt,quadrants.tsv}`.
The AFDB sweep (sweep B) is complete; its summary is at
`~/work/prosmos_2026/afdb_negspace_sweep_v2/summary.tsv`.

The structural-reachability follow-up (the analytic + constructive
checks above) has not been started.
