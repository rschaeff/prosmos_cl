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
  symmetry group, leaving mirror pairs as distinct skeletons),
  `handedness_signature` (paper §1.1.2 scalar-triple-product per-triple sign),
  `dist_sum` + `layers` + `rcc_dedup` (RCC port from decoded CG-2012 IL —
  see "RCC port" section below)
- `combine_b.py` — experimental Model B (combinatorial adjacency)
  enumeration; reaches K1,4 / tripod+pendant declared graphs but
  over-counts and misses C5. See "Model A vs Model B" below.
- `prosmos.py` — ProSMoS query.txt writer from `oracle.SSPRecord`. Used
  to auto-generate design-target queries (`scripts/generate_design_target_queries.py`).

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

**Combine-pairs Phase C status:**

| Dim | Phase B | Phase C1 | Phase C-S4 | Phase C-S5 | Phase C2 | **Phase C2-revert** | oracle (distinct skel_id) | paper |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 3 | 4   | 5   | 5   | 5   | 10  | 5   | 11   | 23   |
| 4 | 10  | 14  | 42  | 42  | 84  | **42** | 41   | 221  |
| 5 | 41  | 70  | 72  | 198 | 396 | **198** | 648  | 1239 |

Phase C2 (added `start_up` to `canonical_key`) was reverted. Pre-C2 S4=42
matched oracle 41 almost exactly; post-C2 we overshot to 84. The paper
says UP-start and DOWN-start give distinct linker-plane assignments
(so distinct SSPs at the SSP level), but CG-2012's oracle counts them
as one skeleton at the skeleton-enumeration level. `Skeleton.start_up`
and `Skeleton.orientations` are retained for downstream handedness /
chirality work (and `handedness_signature` from Phase C4 uses them);
the field is just no longer part of the canonical equivalence.

S4 now essentially matches the oracle (42 vs 41; +1 likely an RCC
tie-break case). S5 remains 3.3× short.

**Model A vs Model B for adjacency** (resolution of the S5 gap):

There are two paper-supportable readings of what counts as a skeleton edge:

- **Model A** (geometric, default): edges = pairs of hex-adjacent lattice
  points. The full geometric adjacency graph of any placement is the
  skeleton's edge set automatically. Implemented by `combine.py`;
  cannot reach K1,4 / C5 (not 2D-hex-realizable as induced subgraphs).
- **Model B** (combinatorial, experimental): edges are declared during
  combine. Each step adds one new node with one explicitly-declared
  join edge; geometric coincidences are not declared (would be `X` in
  the ProSMoS query). Implemented in `combine_b.py`; reaches the K1,4
  and tripod+pendant declared graphs Model A can't (but doesn't reach
  C5 yet — cycles need declaring 2+ edges per add).

`Skeleton.edges` (`frozenset[tuple[int,int]] | None`) is the data-structure
piece: when set, `adjacency_matrix()` uses the declared edges; when None,
falls back to geometric. Both modes coexist; existing Model A code paths
are unchanged.

**Phase 2 status (Model B enumeration):** experimental.

  | Dim | Model A | Model B | oracle |
  |---:|---:|---:|---:|
  | 3 |   5 |    19 |    11 |
  | 4 |  42 |   358 |    41 |
  | 5 | 198 |  4842 |   648 |

Model B's enumeration over-counts (3858 P5 variants at S5 alone) because
`canonical_key` quotients by lattice rotation but the declared graphs
have many lattice realizations the rotation quotient doesn't merge. The
right dedup is graph-isomorphism on declared edges; future work.

For the **14 design-target work**: `prosmos.py` writes a ProSMoS
query.txt directly from an `SSPRecord` (oracle data — already has
adjacency, types, sheets, handedness). The script
`scripts/generate_design_target_queries.py` produces queries for the 9
of 14 design targets present in `~/chalam/CG-2012/S5/IA.txt`; output
lands in `../example/ssp_design_targets/queries_enum/`. The 5 missing
targets (`5-283-1-2`, `5-307-1-2`, `5-243-1-2`, `5-265-7-7`,
`5-234-7-7`) are from a post-2012 website enumeration not in
CG-2012 — hand-construction via `Skeleton.edges` + the writer is the
path for those.

**S5 oracle gap analysis (post-investigation, pre-Model-B):** the gap is structural.
Direct per-grid comparison of combine vs oracle labeled-adjacency
patterns shows the two produce **different sets of unlabeled grids**
— only P5 overlaps:

| Unlabeled grid | Combine | Oracle |
|---|--:|--:|
| K1,4 star (4e)        | 0  | 5  |
| tripod+pendant (4e)   | 0  | 58 |
| P5 (4e)               | 12 | 45 |
| triangle+2pendants (5e)| 24 | 0  |
| C5 cycle (5e)         | 0  | 12 |
| K1,4 sparse (6e)      | 9  | 0  |
| K1,4 dense (7e)       | 30 | 0  |

We verified exhaustively (box [-3,3]²) that pure C5 is **not
realizable as a hex-induced graph** — closing a 5-cycle on hex always
introduces a shortcut edge. The same applies to pure K1,4 (4 leaves
of any vertex have at least 2 leaf-pair adjacencies on hex). Combine
respects strict hex realizability; oracle CG-2012 evidently treats
skeletons as abstract adjacency graphs without enforcing it, so
records "C5 lattice" or "K1,4 lattice" exist there but their pairs
with `X` markings would have to be lattice-adjacent — contradicting
paper §1.1.1 ("Non adjacent SSEs can either interact optionally (X)
or not interact at all (-)").

The 198 vs 648 gap therefore isn't fixable by tweaking `canonical_key`
or combine logic — it reflects a different definition of "skeleton
lattice". Our 198 are paper-Methods-correct hex-induced arrangements;
matching CG-2012's 648 would require accepting abstract non-realizable
skeletons. See `project_s5_oracle_gap.md` in session memory for the
full record.

**S4 essentially matches (42 vs 41); S5 closes a major gap (72 → 198).**
WHITELIST_S4 and WHITELIST_S5 are now derived directly from paper Fig. S3
+ hex lattice geometry, replacing the earlier oracle-X-included extractions
(which incorrectly included K4 / K5 / K5-1e — none of which are realizable
on the 2D hex lattice).

S4 lattice whitelist:
  - **(a)** P4 path:               `((0,1),(0,2),(1,3))`         — 3 edges
  - **(b)** triangle + pendant:    `((0,1),(0,2),(0,3),(1,2))`   — 4 edges
  - **(c)** C4 + short diagonal:   `((0,1),(0,2),(0,3),(1,2),(1,3))` — 5 edges

S5 lattice whitelist:
  - **(d)** P5 path:                `((0,1),(0,2),(1,3),(2,4))`           — 4 edges
  - **(e)** triangle+2pendants:     `((0,1),(0,2),(0,3),(1,2),(1,4))`     — 5 edges
  - **(g)/(h)** sparse K1,4:        `((0,1),(0,2),(0,3),(0,4),(1,2),(3,4))` — 6 edges
  - **(f)** / dense K1,4:           `((0,1),(0,2),(0,3),(0,4),(1,2),(1,3),(2,4))` — 7 edges

Combine S5 distribution across the 4 grids: 12 (d) + 98 (e) + 28 (g sparse) + 60 (f/g dense) = 198.

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

**Phase C3 attempt: reverted.** Sequence-reversal + start_up-flip
collapses Bent-A/Bent-B (sequence reversals of each other) and CW/CCW
triangle mirrors — both kept distinct by CG-2012 and the paper.
Per paper Appendix §1.1.2, handedness is the scalar triple product
`(p × q) · r`, which is sign-sensitive to argument order; sequence
reversal flips every triple's handedness sign, so it's not a valid
quotient under the paper's label-by-label equivalence definition.

**Phase C4 investigation: handedness_signature utility added; no
canonical_key change.** Paper §1.1.2 defines equivalence as
label-by-label matching handedness signatures. Empirically this
over-collapses against the oracle:

| Dim | enum_skeletons (rotation + start_up) | pure handedness eq | oracle |
|----:|---:|---:|---:|
| 3 | 10 | 3 | 11 |
| 4 | 84 | 23 | 41 |
| 5 | 396 | 193 | 648 |

CG-2012's oracle counts per-labeling without applying the paper's
handedness equivalence rule. Our current canonical (rotation + start_up
+ chirality) sits closer to oracle behavior than pure handedness would.
`handedness_signature(skel)` is retained as a utility — it'll be
load-bearing once SSE-type assignment + ProSMoS query writing land.

**RCC port (Chitturi 2016 Appendix §1.2)**: `combine.dist_sum(skel)`,
`combine.layers(skel)`, `combine.rcc_dedup(skels)`, and
`enumerate.enumerate_skeletons_rcc(n)` implement the paper's Relative
Compactness Criterion. The port was derived by decompiling CG-2012's
binary; full reference at [docs/cg2012_internals.md](docs/cg2012_internals.md).

Production `enumerate_skeletons` does *not* apply RCC. Empirical
comparison:

| Dim | enumerate_skeletons | enumerate_skeletons_rcc | oracle |
|---:|---:|---:|---:|
| 3 | 5   | 3   | 11  |
| 4 | 42  | 23  | 41  |
| 5 | 198 | 164 | 648 |

CG-2012 applies RCC in a **two-tier** way: locally inside
`extend()` for each parent's extension batch, globally only on the
"carry-forward" set (motifs retained as building blocks for larger
dimensions). The motifs *emitted as output* at each dim get only
local-tier dedup — which is why oracle records keep some
handedness-equivalent siblings from different parents. Closing the
last +1 at S4 would require modeling the output/carry-forward split
end-to-end; out of scope for now. See the docs for details.

Next:
1. **Decide whether to revert Phase C2's start_up doubling**. Pre-C2,
   S4=42 matched oracle 41 almost exactly; post-C2 we overshoot to 84.
   The paper says UP-start and DOWN-start give distinct linker planes
   (so distinct SSPs), but oracle treats them as one. Reverting C2
   would tighten S4 match at the cost of throwing away orientation info
   that's needed downstream.
2. **Investigate the S5 oracle gap directly** — sample specific oracle
   S5 records, compare against our enumeration, find what we're
   missing structurally (lattice variants? labeling conventions?).
3. **RCC tie-breaking** (paper Appendix §1.2: pair-wise distance sum,
   then unique-y). Picks canonical representative among equivalent.
   Useful once we settle on the right equivalence.

## Skeleton → ProSMoS query pipeline (S3-S5)

`assignment.skeletons_to_records` closes the gap from "we have N
hex-induced skeletons at dimension k" to "we have a corpus of ProSMoS
query.txt files we can search against an ECOD metamatricesDB". For
each skeleton it yields one `SSPRecord` per Cartesian H/E type
assignment, with:

- **Interaction matrix** derived from lattice adjacency, types, and
  the sequence-direction alternation in `Skeleton.orientations`.
  Codes follow the paper conventions plus the lowercase `u`/`v`
  helix-helix variants seen in CG-2012 IA.txt (see module docstring).
- **Sheet partition** as connected components of E-E lattice
  adjacency: one `sheetS` directive per sheet of ≥ 2 strands, one
  `sheetD` directive per pair of E SSEs in different sheets.
- **Per-triple handedness** from `handedness_signature(skel)` — one
  directive per (i, j, k) triple with non-zero signed triple product.
  Planar acyclic shapes whose triples are all coplanar emit no
  handedness directives, matching the looser-search intent.
- **Length constraints** via `prosmos.write_query` (per-SSE, defaults
  E ≥ 5, H ≥ 8, max 1000 — see [[feedback-prosmos-query-length]]).

The driver `scripts/generate_enumerated_queries.py` writes the full
S3-S5 corpus to `../example/ssp_enumerated/queries_typed/sN/` (gitignored
because the 7,048 files are reproducible from this script):

| Dim | Skeletons | × 2^N typings | = queries |
|---:|---:|---:|---:|
| 3 | 5    | 8   |    40 |
| 4 | 42   | 16  |   672 |
| 5 | 198  | 32  | 6,336 |
| total |  |  | **7,048** |

This is the S5-bounded scope. Extending to higher N is the same code
path — `enumerate_skeletons(6)` already returns 2,372 skeletons in
under a second; the only added cost is the per-skeleton Cartesian
explosion (S6 → 2,372 × 64 = 151,808 queries).
