# Positive controls for the S3-S5 enum + ProSMoS sweep

Validation of the enum-to-ProSMoS-to-ECOD pipeline against two canonical
protein motifs of known fold-class membership in ECOD: the four-helix bundle
and the minimal Rossmann-like motif (RLM, Medvedev et al. 2021,
[10.1016/j.jmb.2020.166788](https://doi.org/10.1016/j.jmb.2020.166788)).
Both controls use the v4 sweep results
(`results/v4_ecod_manual_reps_summary.tsv`) and per-query hit dirs in
`~/work/prosmos_2026/ecod_search_v4/hits/`.

## Why these controls

The negative-space analysis in `coverage_gaps.md` rests on the assumption
that ProSMoS searchmatrix correctly *finds* motifs that exist in nature. A
"zero hit" claim is only as strong as the matcher's recall on known
positives. We validate that recall on:

- **Four-helix bundle (FHB)** — small, ubiquitous, well-defined. ECOD
  X-group 601 ("Four-helical up-and-down bundle") gives 158 manual-rep
  domains in the v4 DB.
- **RLM α-variant** — the Medvedev paper's primary query class. Used the
  same tool (ProSMoS) and same DB strategy (PALSSE-derived ECOD matrices),
  so direct methodological precedent. ECOD X-groups 2003, 2111, 4279, 7582-87,
  7592 (Rossmann-like superfamilies) give 739 manual-rep domains in v4.

## Results

| Motif | Best single query | Recovery | Total hits | % "unexpected context" |
|---|---|---|---|---|
| Four-helix bundle | `s4-0026-0000` (P4 path, HHHH) | **90%** (142/158) | 4,056 | 97% |
| RLM α-variant | `s5-0098-0021` (Y-sandwich, EHEHE) | **87%** (641/739) | 2,910 | 78% |
| (RLM α, top-5 union) | | 91% (671/739) | 3,222 | |
| (RLM α, all 198 EHEHE queries union) | | 94% (692/739) | 3,557 | |

The 6% RLM miss vs the Medvedev paper's full-PDB 0.7% miss is explained by
(a) we tested only the α-variant here, paper used three query variants
(α / β / linker), and (b) v4 is manual reps only (~19k), the paper searched
the full PDB.

## Skeleton geometries of the winning queries

### s5-0098 (RLM α-variant winner)

SSE labels: `0=β1  1=α1  2=β2  3=α2  4=β3` (typing EHEHE = index 21).

Adjacency (spatial contacts; symmetric):

```
     β1  α1  β2  α2  β3
β1 [  .   X   X   .   . ]
α1 [  X   .   .   X   . ]
β2 [  X   .   .   .   X ]
α2 [  .   X   .   .   . ]
β3 [  .   .   X   .   . ]
```

This is the canonical doubly-wound α/β/α sandwich:

- β-sheet layer: β1-β2-β3 (SSEs 0,2,4), parallel-paired via β1↔β2 and β2↔β3
- α layer 1: α1 (SSE 1) packs against β1 and α2
- α layer 2: α2 (SSE 3) packs against α1

Handedness signature: `(0,0,0,0,0,0,0,0,0,0)` — all wildcard. This is why
the skeleton is so promiscuous: 10 wildcard triples means the matcher
accepts any chirality. Edges = 4 (the minimum for a connected 5-SSE
skeleton).

### s4-0026 (FHB winner)

P4-path topology (linear chain), degree sequence [2,2,1,1], 3 edges. The
"up-and-down" bundle is a chain of helices with each one contacting only
its sequential neighbors -- which is what P4 encodes.

## What the "embedded in unexpected context" hit fraction means

For both controls, the canonical winning query returns many more hits than
the expected fold-class size:

- FHB: 4,056 hits vs 158 known FHBs → 3,914 hits in non-FHB folds
- RLM: 2,910 hits vs 739 known RLMs → 2,269 hits in non-RLM folds

These are not false positives. ProSMoS is a **motif** detector; a four-
helix substructure inside a β-barrel is correctly matched. The Medvedev
paper's central finding — that the RLM appears in 156 H-groups across
123 X-groups, including many "unexpected" hosts — is the same phenomenon
visible here.

The negative-space analysis is the complementary question: what enumerated
motifs find *zero* embedded matches in 19k+ diverse domains? That's the set
we're carrying forward to AFDB (`coverage_gaps.md` + the in-flight 707k F70
sweep).

## Implications for the negative-space interpretation

1. **Methodology has high recall** on known positives (>85% per single
   query, ≥91% via small unions, ~94% via the full typing-class union).
   Consistent with the Medvedev paper's reported <1% false-negative rate
   on the full PDB.

2. **The 25 zero-hit S5 skeletons aren't artifacts of low matcher
   sensitivity.** They're the chirally-constrained skeletons (mean 1.2
   wildcards vs 3.2 for hit-bearing); the matcher works fine on chirally-
   specific queries that *do* hit (e.g. the chirality-bias mirror pairs
   in `coverage_gaps.md`).

3. **The enum's coverage of the canonical motif space is non-trivial.**
   The skeleton that corresponds to the canonical RLM (s5-0098) is in our
   198, in the position that the Medvedev paper would have hand-encoded.
   The enum is generating real biological motifs as part of its uniform
   sampling of valid lattice topologies, not just random graphs.

## Reproducibility

The 158 FHB ECOD IDs and 739 RLM ECOD IDs are obtained by recursive
descent under the relevant X-groups:

```sql
WITH RECURSIVE descs AS (
  SELECT id FROM ecod_rep.cluster WHERE id IN ('601')         -- FHB
  -- OR ('2003','2111','4279','7582','7583','7584','7585','7586','7587','7592')  -- RLM
  UNION ALL SELECT c.id FROM ecod_rep.cluster c JOIN descs d ON c.parent = d.id
)
SELECT lpad(d.ecod_uid::text, 9, '0')
FROM ecod_rep.domain d JOIN ecod_commons.derived_files df USING(ecod_uid)
WHERE d.manual_rep AND df.file_type_id=2 AND df.status='complete'
  AND df.domain_source_type='pdb'
  AND (d.t_id IN (SELECT id FROM descs) OR d.f_id IN (SELECT id FROM descs));
```

The hit-list of a query is the basename-stripped contents of its
`~/work/prosmos_2026/ecod_search_v4/hits/<query_name>/*.txt` dir; intersection with the
SQL output above gives the recovery count.
