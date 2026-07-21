# Session handoff 2026-07-20/21 (c) — Ruczinski "forbidden" 4-strand topologies

**STATUS: COMPLETE.** Sweeps finished, census computed. Earlier work this session
(design scaffold, instrument agreement, n=3/4/5 curve, cluster-unit significance,
prosmos_inspect significant-browse) is in `session_2026-07-20b_handoff.md` +
memory `project_{design_scaffold,instrument_agreement,why_s5_curve}`.

## Origin

Nick observed that **A0A085ZJC4_D1 (an AlphaFold model) realizes Ruczinski
topology 23** — a four-strand β "pretzel" (spatial strand order 1-4-2-3) that
Ruczinski et al. (Proteins 2002 48:85-97; PDF + `fig1.png` in `prosmos_cl/doc/`)
Fig 1 catalogued as one of **48 four-strand motifs that NEVER occurred** in the
2002 Dunbrack database. Question: do those 48 occur in the modern PDB and/or AFDB?

## RESULTS (final, both filters applied)

Of Ruczinski's 48 "never observed" four-strand topologies:

| | of 48 |
|---|---|
| occur in **TRUE experimental** structures | **26** |
| occur in **AFDB** | **47** |
| **predicted-only** (0 experimental) | **21** |
| **absent everywhere** | **1** |

*(Computed over the 48 ACTUAL Fig-1 panels, transcribed 2026-07-21. The earlier
31/47/16/1 was computed over the wrong 48 — see "Forbidden set" below.)*

**The forbidden set is largely not forbidden.** Their 2002 absence was a statement
about database size, not a structural prohibition. Per-motif detail in
`ruczinski/sweep_occurrence.json`.

### The key negative — S5 vs S4 (RS's diagnosis, confirmed)

| object | experimental instances |
|---|---|
| **S5 `sk171 ty07`** (2 helices + 3 strands, HHEEE) | **0** (of 500,313 searched) |
| **S4 topology 23** (`motif_15`, the pretzel sheet alone) | **21** (+224 pred, +321 AFDB) |

Encoding "topology 23" as an S4 **silently swapped a rare object for a common
one**. The β-sheet that caught the eye is ordinary; the PDB-dark object is the
full 5-SSE motif *with the helices*. Any future framing must keep these separate.

## THREE CORRECTIONS that changed conclusions (read before trusting any rerun)

1. **Sheet completeness.** A query saying "these 4 strands share a sheet" also
   matches a 4-strand SUBSET of a larger sheet. Ruczinski's motifs are COMPLETE
   four-stranded sheets. Verified: PDB record `000011053` matched 4 SSEs of an
   ELEVEN-strand sheet. Filter: matched SSEs must equal the members of a
   4-strand sheet (`sheets_{pdb,afdb}.pkl`, built by `sheet_index.py`).
2. **Provenance.** `derived_files.domain_source_type='pdb'` is **FILE FORMAT, not
   provenance**. The 1.77M "full redundant PDB" build is only **28.3%
   experimental** (68.3% afdb). Hits must be split via
   `domain_summary.source_type` (`provenance.pkl`). *(Third time hitting this trap
   in one session — see `feedback_source_type_filter`.)*
3. **Combined filter impact: 555,699 raw hits → 12,050 experimental + 36,099
   predicted (~91% were artifacts).** Without both filters the headline reads
   "all 48 forbidden topologies occur tens of thousands of times" — entirely false.

## Caveat limiting the "16 predicted-only"

Those 16 are **not** safely "forbidden yet predicted": the experimental side
searched ~500k domains = only **33.5% of ECOD's experimental domains**; the other
991,648 have no structure file at all. So "0 experimental" can mean "absent from
the searchable third." Documented in
**`/data/ecod/database_versions/EXPERIMENTAL_STRUCTURE_GAP.md`** (with the correct
SQL and the file-format-vs-provenance trap).

Corollary: **Nick's ask #1 ("search the full redundant experimental PDB") has no
answer beyond ~500k** — that IS the ceiling. The original 496k sweep was already
within 1% of it. Deepening requires *generating* ~991k structure files.

## Artifacts — `~/work/prosmos_2026/ruczinski/`

- `queries/motif_01..96.query` + `motifs.tsv` (motif_id, spatial_seq, orientation,
  code_matrix, jumps, crossings, `ruczinski_forbidden`).
- `sweep_pdb/`, `sweep_afdb/` — completed sweeps (500 + 800 chunks, no failures).
- `analyze_sweep.py` — the census, **both filters baked in**; rerun = one command.
- `sheet_index.py`, `sheets_{pdb,afdb}.pkl`, `provenance.pkl` — the filter inputs.
- `validate/` — topology-23 encoder validation (tight+loose, negative control).
- `../s5_fullpdb_test/` — the S5 `sk171`/`sk088` vs full-DB test (job 630129).
- `doc/fig1_panels/` (+`_big`, `fig1_rows/`) — Fig 1 split into 48 panels.
- `ecod_db_pdb_full/metamatricesDB{,.clean}` — the 1.77M build (mixed provenance;
  1,749,222 clean records after removing 23,650 malformed).

## The encoder (VALIDATED — reusable)

Motif = (spatial strand order, up/down dirs) → ProSMoS query: 4 `E`s; code matrix
`c`(parallel adj)/`t`(antiparallel adj)/`-`(non-adj); `sheetS 1 2 3 4`; lengths.
**Handedness** (tight variant) via the convention decoded from
`searchMatrix/src/searchControl.h::chirality`: align each strand's direction to
strand 1, average → **N**; **M** = (m₂−m₁)×(m₃−m₁) from midpoints; **N·M>0 → R,
<0 → L**. Python reimpl matches the binary (topology 23 → R,R,L,L; self-hits;
wrong handedness → 0 hits).

**Loose vs tight:** the code matrix is the COMPLETE invariant (384 perm×dir → 96
distinct = Ruczinski's ¼·4!·2⁴), so the **loose query is the faithful motif
encoding**; handedness is a per-structure chirality refinement, not part of a
Ruczinski motif. **The sweep used loose.**

## Forbidden set — TRANSCRIBED, not derived (the ≥2-jump inference was WRONG)

The 48 are now read directly off Fig 1 (`ruczinski/panels_user.tsv`, all 48 panels
read by RS; `verify_panels.py` validates). **48/48 give distinct code matrices.**

**The earlier derivation "Fig 1's 48 = the ≥2-jump motifs" is REFUTED.** Fig 1
contains one 0-jump motif (panel 1) and thirteen 1-jump motifs; the jump histogram
over the true 48 is **1/13/26/8** for 0/1/2/3 jumps (the full 3-jump class is
present, consistent with "3 jumps never occur"). The ≥2-jump set overlapped the
true set by only **34/48** — 14 wrong, 14 missing — and inflated experimental hits
**24×** (12,050 vs 501). The 40+8=48 count match was coincidence.

### NOTATION — both columns are "as drawn" (three failed readings before this)

- `strand_order` = spatial order of sequence-strands **left-to-right as drawn**.
  "Topology 23 = pretzel = 1-4-2-3" is this quantity. It is **not** the inverse
  permutation (position of sequence-strand *i*); reading it that way makes
  physically distinct panels collide under the 180° rotation symmetry.
- `ori_drawn` = up/down **read left-to-right as drawn**, NOT indexed by sequence
  strand. Sequence-indexing sends panel 23 to motif 4 and breaks the anchor.

**Anchor (the discriminator): panel 23 → code matrix `--ttc-` → motif_15**, which
is both the structurally-validated A0A085ZJC4 encoding and a match to the name
"topology 23." Only the as-drawn/as-drawn reading satisfies it.

There is **no canonicalization** in the figure (panels 33/34/37/39/43/44/47/48
have the leftmost strand pointing down), and the panel ORDER is not reconstructible
from (jumps, strand_order, #down) — block order is not numeric (1423 precedes 1342)
and #down is not the within-block tiebreak. Ordering is advisory; **distinctness is
the hard gate**.

## Operational notes

- Always `db_validate.py <db> --clean <out>` before searching — it exits **1 when
  it finds** malformed records (a signal, not a failure); do not let `set -e` kill
  the caller (this broke the first launcher run).
- Build file lists with `find`, never `ls` (aliased to long format — it killed the
  first S5 test submission; the `array_inverted.sbatch` guard caught it instantly,
  zero compute wasted).
- Use the post-d3e1e83 `searchMatrix/build/searchmatrix`; `unset OMP_PROC_BIND`;
  strip SLURM/PMI vars for generateMatrix.

## Open / possible next steps

- The 16 predicted-only motifs need the experimental structure gap closed (or a
  targeted per-motif check) before any "forbidden yet predicted" claim.
- The S5 object (`sk171 ty07`, 0/500,313 experimental) is the real rare finding —
  it is one of the 5 depth-significant AFDB-only cells; see
  `project_ruczinski_forbidden` and the design-scaffold thread.
