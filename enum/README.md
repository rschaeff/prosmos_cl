# SSP enumeration

Fresh Python implementation of the super-secondary structure pattern (SSP)
enumeration model from Chitturi B, Shi S, Kinch LN, Grishin NV.
*Compact Structure Patterns in Proteins.* J Mol Biol 428(21):4392–4412 (2016).
DOI [10.1016/j.jmb.2016.07.022](https://doi.org/10.1016/j.jmb.2016.07.022).

## Why this exists

Chalam Chitturi's 2012-vintage C# enumerator (`CompactGenerator2.exe`, recovered
from his archived home directory) produces a 1472-SSP set at dimension 5;
the published paper reports 1239. Running a decade-old `.exe` from a former
lab member's archive as the input to new publishable work is not defensible —
this package re-implements the enumeration in modern, readable, paper-traceable
code so every decision can be tied back to a figure or paragraph in the paper.

The 2012 binary is *not* a dependency. It is used only as an **output oracle**:
its IA.txt files (in `reference/`) define the SSP set we must match, line by
line, before our implementation is considered correct on the same rules.

## What ships here

```
enum/
├── README.md               this file
├── pyproject.toml          package metadata
├── src/ssp_enum/           library code
│   ├── lattice.py          hexagonal lattice, LatticePoint, neighbors
│   ├── node.py             CGNode (directed SSE node on the lattice)
│   ├── motif.py            CGMotif (an SSP); isMotifCompact, surrCompact
│   ├── matrix.py           CGMatrix (2D interaction matrix, IA.txt I/O)
│   ├── enumerate.py        the dimension-by-dimension growth loop
│   ├── compactness.py      PCC / SCC criteria from Appendix Fig. S2
│   ├── prosmos.py          ProSMoS query.txt writer (length constraints, etc.)
│   └── oracle.py           IA.txt parser for validation against CG-2012
├── reference/              CG-2012 outputs used as validation oracle
│   └── IA-*.txt            symlinks or copies from ~/chalam/CG-2012/Sn/IA.txt
└── tests/                  pytest suite
    ├── test_lattice.py     neighbor enumeration matches paper Fig. 1a
    ├── test_compactness.py PCC/SCC against Appendix Fig. S2 cases
    └── test_oracle.py      our enum output matches reference/ counts
```

## Validation strategy

Each dimension N has a known SSP count from the CG-2012 reference:

| Dim | CG-2012 count | Paper count | Status |
|---:|---:|---:|---|
| 2 | 1 (skeleton) / 5 (typed) | 5 | binary counts skeletons, paper counts typed assignments |
| 3 | 11 | 23 | binary smaller — rule difference |
| 4 | 203 (IA.txt blocks) / 41 unique skeletons | 221 | close |
| 5 | 1472 | 1239 | binary larger — paper added filtering |

Our implementation targets matching the **binary's** counts initially
(reproducing the 2012 algorithm faithfully), then extends to match the
**paper's** counts by applying the additional rules described in paper
Methods that the 2012 binary apparently lacked.

This two-step lets us catch implementation bugs vs. methodology bugs:
divergence at step 1 means we got the lattice/compactness wrong;
divergence at step 2 means we have the algorithm right but missed a
paper-specified post-filter or refinement.

## Paper section ↔ code module map

| Paper | Code | What's implemented |
|---|---|---|
| Methods, "Generating SSPs with a lattice" | `lattice.py` | Hexagonal lattice over Z = {−1, 0, 1}, neighbor enumeration |
| Methods, "compactness was determined by three criteria" | `compactness.py` | PCC (perimeter), SCC criterion-1 (collinearity / ≥2-adjacency), specific eliminations |
| Appendix Fig. S2c | `compactness.py` | The trapezoid-vs-rhombus elimination |
| Appendix Fig. S2d | `compactness.py` | Specific S5 eliminations |
| Appendix Fig. S3 | `lattice.py` | Allowed grids at N=4 and N=5 |
| Methods, "Assignment of SSE type" | `node.py` | H/E assignment to each lattice node |
| Methods, "specification of interactions" | `matrix.py` | Pairwise interaction type per SSE pair |
| Methods, "fold growth through stepwise addition" | `enumerate.py` | Combine S_p + S_q → S_{p+q}, dedupe via mirror symmetry |
| (out of scope here — search step) | `prosmos.py` | Write final SSPs as ProSMoS query.txt with length constraints |

The `oracle.py` module reads the abbreviated IA.txt format
(`sS`/`sD`/`hand` keywords, `[5-skel-X][sub-Y]` headers, upper-triangular
matrix) and yields canonical SSP records for comparison with our enumeration's
output. This is the only place where the legacy format is consumed; everywhere
else we use ProSMoS-format directly.

## What we are NOT reimplementing

- **ProSMoS search itself.** Lives in the parent repo (`searchMatrix/`).
- **The post-ProSMoS hit filter** (paper Appendix Fig. S4). Separate concern.
- **PDB → SCOP/ECOD mapping.** Handled by parent repo `scripts/`.

This package outputs *queries*. Downstream pipelines consume them.

## Setting up `reference/`

The validation oracles in `reference/` are recovered from Chalam Chitturi's
archived home directory at `~/chalam/CG-2012/` and are not committed to this
repo. Populate locally:

```
cd reference/
ln -s ~/chalam/CG-2012/S3/IA.txt IA-S3.txt   # 0 bytes — S3 only logged in Info.txt
ln -s ~/chalam/CG-2012/S4/IA.txt IA-S4.txt   # 203 SSP blocks
ln -s ~/chalam/CG-2012/S5/IA.txt IA-S5.txt   # 2807 SSP blocks (1472 unique)
```

If `~/chalam/` is not on your machine, contact the lab — the CG-2012 tree
was restored from archive in May 2026 and lives on the leda fileserver.

## Status

S3 full enumeration + SCC-2 whitelist + combine-pairs Phase A/B + Phase C1 land green (31 tests pass):

```
tests/test_oracle.py     ✓  4 passed
tests/test_enumerate.py  ✓  7 passed
tests/test_grids.py      ✓  7 passed
tests/test_combine.py    ✓ 13 passed   (single+multi node combine, canonical-key
                                        rotation-only quotient with mirror-pair
                                        distinguishing, S3=5 S4=14 S5=70 regression
                                        pins, S2+S2 / S3+S2 contribution checks)
```

Implemented:
- `oracle.py` — IA.txt parser, yields SSPRecord stream
- `lattice.py` — hex axial coords + neighbor enumeration (in-plane + cross-layer)
- `skeleton.py` — ordered LatticePoint tuple + adjacency matrix + chirality label
- `compactness.py` — PCC, SCC-1, and SCC-2 (Appendix Fig. S3 whitelist)
- `grids.py` — unlabeled-adjacency-graph canonical form + per-dimension whitelists
  (S3: 2 grids; S4: 3 grids; S5: 4 unique signatures, see g≡h note below)
- `enumerate.py` — dimension dispatcher; `enumerate_s3_planar` (4 base shapes);
  `enumerate_s3` (11 SSPs = planar + chiral); `enumerate_skeletons(n)` via combine,
  iterating all (p, q) splits with p+q=n, p≥q
- `combine.py` — Chitturi 2016 Appendix §1.1 + §1.2 combine-pairs growth:
  `valid_extension_points`, `combine_with_single_node` (|s2|=1),
  `combine_two_skeletons` (|s2|≥2: anchor × ext_pt × 60°·k × {front,end} join),
  `canonical_key` (Phase C1: rotation-only quotient — translation × 6 hex
  rotations; reflections and z-flip flip handedness and so stay outside the
  symmetry group, leaving mirror pairs as distinct skeletons)

S3 planar produces the 4 base spatial-sequence patterns from CG-2012's
S3/Stru.txt. Crossing each base shape with chirality (per `enumerate_s3`)
gives the full 11 SSPs:

| Adjacency (1-2, 1-3, 2-3) | Name | CG-2012 panels (chirality) |
|---|---|---|
| (T, F, T) | Linear   | 3-0 (None), 3-5 (L), 3-6 (R) |
| (T, T, T) | Triangle | 3-1 (L), 3-2 (R)             |
| (F, T, T) | Bent-A   | 3-3 (None), 3-7 (L), 3-8 (R) |
| (T, T, F) | Bent-B   | 3-4 (None), 3-9 (R), 3-10 (L)|

Rule: acyclic shapes (Linear, Bent-A, Bent-B) admit an unhanded variant
in addition to L/R; the cyclic Triangle has only L and R variants — the
closed 3-cycle is intrinsically chiral via the orientation of the
sequence walk around the loop. Chirality is recorded as a `'L' | 'R' |
None` label on `Skeleton`, mapping directly to CG-2012's `hand i j k L|R`
and ProSMoS's `handedness i j k L|R` lines. The full Z-aware geometric
realization (assigning concrete z ∈ {-1, 0, 1} per node) is deferred to
SSE-type/direction assignment downstream.

**SCC-2 / S5-grid note:** Appendix Fig. S3 lists 5 induced grids at S5
(d–h), but the Appendix's per-grid handedness table says Grid 5 (S3-h)
is "Same as Grid 4 (S3-g)" — g and h share an unlabeled adjacency
graph but differ in lattice embedding (mirror pair). Our canonical
form quotients to the adjacency graph, so we observe 4 distinct S5
signatures, not 5. Every oracle SSP record's grid is in the whitelist,
so SCC-2 accept/reject is identical to the paper's. The lattice-
embedding chirality is captured separately by the chirality label
introduced in `enumerate_s3()`.

**Combine-pairs Phase C1 status (rotation-only canonical_key):**

| Dim | Phase B | Phase C1 | oracle (distinct skel_id) | paper |
|---:|---:|---:|---:|---:|
| 3 | 4   | 5   | 11   | 23   |
| 4 | 10  | 14  | 41   | 221  |
| 5 | 41  | 70  | 648  | 1239 |

Narrowing the canonical-form quotient from 12-hex × z-flip to
6-rotation only separates mirror pairs the previous canonical
collapsed. S3 Triangle splits into CW/CCW chirality variants
(4 → 5). S4 +4, S5 +29 — same mechanism. The remaining S4 gap to
oracle 41 is mainly skeleton-level chirality labeling (some
skeletons admit L and R distinct, others are self-mirror and stay
1). S5's 648 - 70 = 578 gap is dominated by z-displacement
variants the paper allows ("up to three layers") which our
enumeration doesn't yet produce — seeds + valid(s) extensions all
stay in the z=0 plane.

Next (remaining Phase C):
1. SSE orientation tracking (up/down per node, alternating along
   the sequence per paper §1.1). Phase C2.
2. Handedness-based equivalence: two skeletons equivalent iff for
   every triple of node labels the (p × q) · r handedness sign
   matches. Phase C3 — requires (1).
3. RCC tie-breaking on equivalent skeletons (paper Appendix §1.2:
   pair-wise distance sum, then unique-y count). Phase C4 —
   requires (2).
4. Cross-layer / z-displaced extensions: revisit `valid_extension_points`
   and the seeds to allow placement at z ∈ {-1, +1}. Closes the
   bulk of the remaining S5 gap. Phase C5.
5. SSE-type (H/E) and interaction-type assignment to grow
   skeletons into full SSPs. Post-Phase-C.
6. ProSMoS-format writer (`prosmos.py`) — depends on (5).
