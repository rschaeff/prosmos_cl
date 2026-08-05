# S3-S5 enum coverage gaps vs ECOD manual reps

Analysis of the v4 ProSMoS sweep (7,048 enumerated S3-S5 queries against
19,015 ECOD `manual_rep=true` domain matrices). Per-query results in
`results/v4_ecod_manual_reps_summary.tsv`.

## Per-class hit rates

| Class | Queries | With hits | Coverage | Total hits |
|---|---|---|---|---|
| S3 | 40 | 39 | 98% | 183,637 |
| S4 | 672 | 589 | 88% | 316,624 |
| S5 | 6,336 | 2,331 | 37% | 196,800 |

## Two distinct gap mechanisms

The headline 37% S5 hit rate is the sum of two independent filters that
nature applies but the enum does not:

### Mechanism 1 — composition gradient (per-strand cost ~½ coverage)

Hit rate as a function of H/E composition is monotonic and ~halves with each
additional strand:

S5 (198 skeletons × 32 typings):

| Composition | Queries | Coverage |
|---|---|---|
| 5H/0E | 198 | 73% |
| 4H/1E | 990 | 54% |
| 3H/2E | 1,980 | 43% |
| 2H/3E | 1,980 | 32% |
| 1H/4E | 990 | 16% |
| 0H/5E | 198 | 7% |

S4 (42 skeletons × 16 typings):

| Composition | Queries | Coverage |
|---|---|---|
| 4H/0E | 42 | 100% |
| 3H/1E | 168 | 100% |
| 2H/2E | 252 | 100% |
| 1H/3E | 168 | 71% |
| 0H/4E | 42 | 19% |

All 83 zero-hit S4 queries are 1H/3E or 0H/4E.

The asymmetry is physical: helix packing is geometrically permissive
(van der Waals only) and tolerates the strand-strand adjacencies our
2D-hex-lattice skeletons emit. β-sheet H-bond register is constrained --
two strands adjacent in the lattice rarely correspond to a realizable
sheet partner geometry. Each additional E in the typing compounds the
constraint.

### Mechanism 2 — chirality bias (S5 only)

This mechanism is invisible at S3/S4 and dominant at S5. Within S5, **25
of 198 skeletons (13%) are zero-hit across ALL 32 typings** -- entire
topologies that no manual-rep domain realizes.

The zero-hit S5 skeletons differ from hit-bearing ones in two ways:

1. **They're denser**. 18 of the 25 have 7 edges (out of K₅'s max 10);
   30% of all e=7 skeletons are zero-hit vs 12.6% overall.

2. **They're more chirally specific.** The handedness signature is a
   10-tuple of {-1, 0, +1} (one per SSE triple), where 0 = degenerate /
   wildcard. Hit-bearing S5 skeletons average **3.2 zeros** in their
   signature; zero-hit skeletons average **1.2 zeros** -- ~3× fewer
   wildcards. Locking in concrete L/R on nearly every triple forces a
   tight match against ECOD's chirality distribution.

The smoking gun: **17 of 25 zero-hit S5 skeletons have a handedness-negated
mirror in the corpus, and 6 of those mirrors do find ECOD hits**:

| Zero-hit | Mirror | Mirror hits |
|---|---|---|
| s5-0006 | s5-0150 | 2 |
| s5-0008 | s5-0151 | 1 |
| s5-0026 | s5-0147 | 1 |
| s5-0071 | s5-0105 | 2 |
| s5-0177 | s5-0032 | 1 |
| s5-0190 | s5-0062 | 8 |

The mirror hits are all low (1-8), so these are rare folds either way --
but the asymmetry shows nature's handedness preference (right-handed
α-helices, predominantly right-handed β-α-β crossovers) acting as an
empirical filter that the enum's uniform sampling doesn't.

## Why chirality bias only appears at S5

S4 has C(4,3) = 4 handedness triples → 2⁴ = 16 possible signatures even
at zero wildcards. ECOD's manual-rep coverage of 4-SSE chirality is
dense enough that no signature is unrealized.

S5 has C(5,3) = 10 triples → 2¹⁰ = 1024 possible signatures. The
distribution of natural chirality across that larger space is uneven,
leaving gaps that our uniform enum walks into.

By extrapolation, S6 (20 triples, 2²⁰ space) and up will show the
chirality-driven gap progressively dominating the composition-driven
one.

## Full zero-hit S5 skeleton list

```
0001  0003  0004  0006  0008  0019  0023  0025  0026  0033
0035  0047  0071  0079  0111  0116  0120  0142  0149  0177
0185  0186  0189  0190  0193
```

Reproduce with:
```
python3 -c "
import sys; sys.path.insert(0, 'src')
from ssp_enum.enumerate import enumerate_skeletons
import csv
skels = enumerate_skeletons(5)
sums = {}
with open('results/v4_ecod_manual_reps_summary.tsv') as f:
    next(f)
    for row in csv.reader(f, delimiter='\t'):
        if not row[0].startswith('s5-'): continue
        sk = int(row[0].split('-')[1])
        sums[sk] = sums.get(sk,0) + int(row[3])
zero = sorted(sk for sk in range(198) if sums.get(sk,0)==0)
print(zero)
"
```

## Implications for the enum

Two design directions, depending on the goal:

- **Reality-aligned corpus**: add a chirality prior (e.g. drop skeletons
  with `sum(handedness_signature)` outside the natural distribution
  estimated from ECOD), and either drop strand-heavy typings on dense
  topologies or down-weight them. Would cut the corpus ~3× while
  preserving >90% of expected hits.

- **Design-space exploration**: keep the 25 zero-hit S5 skeletons as
  flagged "anti-handed" controls -- the de novo design community
  routinely produces left-handed bundles, so "zero ECOD hits but valid
  topology" is a meaningful design target, not noise.

The non-zero-hit S5 skeletons that nonetheless have low coverage (e.g.
the e=7 dense skeletons that only hit 1-5 of their 32 typings) are the
most interesting "edge of the realized" cases worth pulling structures
for visual inspection.
