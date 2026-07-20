# Session handoff 2026-07-20 (b): design deliverable, instrument agreement, n-curve

Three new results built on the corrected full S5 sweeps (see
`session_2026-07-20_handoff.md` for the sweep state). Analysis lives in
`~/work/prosmos_2026/`; scripts mirrored into `enum/scripts/{design_scaffold,
curve,instr_agree}/`. Three companion decks (same visual system) published as
Claude artifacts (URLs below).

## 1. The 200M reframed for design (RosettaCon)  — `enum/scripts/design_scaffold/`

The negative discovery result (no novel local topology) is a positive design
result: the 200M gives sequence-diverse, high-confidence scaffolds over the
shared, experimentally-validated topology vocabulary. The depth argument that
kills novelty makes the scaffold claim.

- **Matched sequence-diversity multiplier = 4.6× median** per shared topology
  (13 PDB vs 58 AFDB sequence families), ~10× total pool, 95.4% of topologies
  gain. Both DBs clustered identically (MMseqs2 50% id / 90% cov). The naive
  ECOD F-group count inflates this to 8× because F-group is coarser than a 50%
  cluster; the correction falls entirely on the PDB side (AFDB grey-reps do not
  merge at 50% — 100% retained).
- **31 pLDDT-validated under-templated cells** (PDB ≤5, AFDB ≥50 confident
  families); 27 survive the strict ≥80 gate. Shared cells → every one is
  PDB-confirmed, no novelty risk.
- pLDDT: `afdb_200m.domain_plddt` is EMPTY; per-protein proxy exact for only 37%.
  Gold standard = mean CA B-factor of the DPAM-chopped structure. The proxy was
  conservative (median 80.8 vs exact 85.2). Pipeline:
  headline → recount_matched → plddt_filter → extract_plddt → validate_perdomain.
- Artifact: https://claude.ai/code/artifact/464273d9-3d09-430e-82bf-b7d00c5dca36

## 2. Instrument agreement: PALSSE across provenance  — `enum/scripts/instr_agree/`

Same-molecule control: 1,529 accessions with one exp + one AF2 ECOD domain in the
same T-group, same store, run through the identical pipeline. Provenance is the
only variable.

- **SSE inventory 100% identical** across provenance.
- Cell agreement median Jaccard 0.75–0.80.
- **Predicted models light +14% more cells** (17.98 vs 15.74/domain; sign test
  z=+8.62) — idealized geometry satisfies more ProSMoS constraints.
- The AFDB-only falsification is **underpowered by construction** (AFDB-only cells
  are PDB-dark; pairs require a PDB domain → only 3 of 619 touched).
- Read: idealization is a second confound on top of depth, same direction — a
  footnote for discovery, harmless for design (multiplier lives on shared cells).
- Method note: shard records across an array (32 tasks); one serial searchmatrix
  over all records is ~6h/side and times out (same chunk_0059 lesson).

## 3. Why S5: the n=3/4/5 curve  — `enum/scripts/curve/`

Computed on identical record sets (491,963 AFDB / 49,640 PDB via
`s5_grid/{afdb,pdb}_sub_records.pkl`).

|  n | cells | AFDB sat | PDB sat | AFDB share (mean/med) | Jaccard |
|---:|------:|---------:|--------:|----------------------:|--------:|
|  3 |    40 |   100%   |  97.5%  |        950 / 1034      |  0.975  |
|  4 |   672 |    89%   |  87%    |        159 / 88        |  0.968  |
|  5 |  6336 |    45%   |  36%    |         29 / 4         |  0.674  |

- Two boundaries: n≤4 saturated (no signal), n≥6 fingerprint (sharing ÷6/step →
  median ≈1 at n=6). **n=5 is the unique unsaturated-yet-still-shared level.**
- **The n=5 divergence is DEPTH, not alphabet.** Decisive control = within-database
  null at matched hitter count: cross-DB Jaccard 0.645 ≈ within-AFDB null 0.664
  (within-PDB 0.791). The DBs differ no more than AFDB differs from itself.
- ⇒ the S5 alphabet is common to experimental and predicted space, scale-robust
  through n=5. **n=6 (24× cost) would measure undersampling, not new vocabulary.**
- S3 is a clean boundary — no enumerator pathology.
- Artifact: https://claude.ai/code/artifact/2bcdd280-3abd-4e53-a2b1-9cd262fd5488

## Also this session
- `enum/docs/ted_vs_prosmos_novel_fold_2026-07-20.md` — TED's novel-fold detection
  vs the S5 grid; SymD is ~46% of TED's novelty output (arrangement, not local
  geometry) — the axis the grid discards. Added as a deck slide.
- Memory: `project_design_scaffold`, `project_instrument_agreement`,
  `project_why_s5_curve` written.

## Still open
- prosmos_inspect branch `corrected-sweep-reexport` UNMERGED.
- 8Å cutoff job (626145) timed out — needs longer limit or split.
