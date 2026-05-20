# CG-2012 binary internals (decompiled reference)

Decoded from `CompactGenerator2.exe` (299KB .NET Framework 4 console app,
2012-01-16 build, recovered from `~/chalam/CG-2012/S5/ia.zip`) using
`monodis` from the `mono` conda env. Source not preserved; this doc
captures the model derived from IL reading.

**Status of CG-2012 as a reference**: per
[[feedback-reproducibility-over-legacy]] memory, the binary is a
non-authoritative artifact. This doc is for *understanding* — to inform
paper-traceable decisions in our Python re-implementation under
`enum/src/ssp_enum/`. Implementing against this doc (bit-matching
CG-2012 behavior) is discouraged; cite the paper instead.

## How to reproduce / extend

```bash
mkdir -p /tmp/cg2012 && cd /tmp/cg2012
unzip -j ~/chalam/CG-2012/S5/ia.zip Debug/CompactGenerator2.exe Debug/CompactGenerator2.pdb
PATH=~/.conda/envs/mono/bin:$PATH monodis --output=dump.il CompactGenerator2.exe
# dump.il is ~137K lines. Useful greps:
PATH=~/.conda/envs/mono/bin:$PATH monodis --typedef CompactGenerator2.exe   # type table
PATH=~/.conda/envs/mono/bin:$PATH monodis --method  CompactGenerator2.exe   # method table
```

If you want C# decompilation rather than IL, install `dotnet` (conda-forge
has it) and `ilspycmd` — that gives readable source with debug symbols.
We stayed at the IL level because the questions we needed answers to
were tractable from method signatures + targeted IL reading.

## Top-level class structure

From the typedef table (`monodis --typedef`):

```
CGElementType  — enum (probably H/E/X)
CGHand         — enum (L/R/None)
CGDirection    — enum (UP/DOWN per node orientation)
CGType         — enum
CGMotif        — main skeleton class
Program        — Main + helpers
CGMList        — list-container of motif-related objects (per dim)
CGNode         — single SSE in the skeleton
LatticePoint   — hex lattice position (xy floats)
CGMatrix       — adjacency/interaction matrix
```

`Program::Main` is the orchestration entry point.

## LatticePoint

```
fields:
  _x : float64       Euclidean x
  _y : float64       Euclidean y

static constants:
  ZERO_THRESHOLD     = 0.01     epsilon for equality
  PRECISION_ROUNDING = 3        decimal places for coord rounding

static methods:
  distance(p1, p2)   = (p1.x - p2.x)^2 + (p1.y - p2.y)^2     # squared Euclidean
  op_Equality(p1, p2) = |p1.x - p2.x| < 0.01 AND |p1.y - p2.y| < 0.01
  op_Inequality       = !op_Equality
```

**Adjacency model = Model A (geometric)**: CG-2012 stores positions in
Euclidean float coords (not axial hex integers), uses epsilon equality
for lookup, and computes neighbors via `(x + cos(i·π/3), y + sin(i·π/3))`
for i=0..5. Two SSEs are adjacent iff their xy positions are at unit
Euclidean distance, within the 0.01 epsilon.

**Speculation on K1,4 / C5 oracle records**: cos/sin-based neighbor
computation can accumulate floating-point error, and `op_Equality`'s
0.01 tolerance is generous in some configurations but tight in others.
The "K1,4" oracle records may correspond to skeletons where geometric
leaf-leaf adjacency *should* exist on a clean hex but the float math
puts them just outside the 0.01 ball, so they end up X-marked instead
of `u`/`v`. Not confirmed by reading the matrix-construction code in
detail.

## CGNode (one SSE in a skeleton)

```
fields (selection):
  Location  : LatticePoint    xy position on the lattice
  PrevNode  : CGNode          predecessor in sequence
  NextNode  : CGNode          successor in sequence
  Direction : CGDirection     UP / DOWN orientation
  Leftmost  : bool            optimization flag for neighbor enumeration

methods (selection):
  enumerateNeighborPositions(bool flag) -> LatticePoint[]
    Returns 6 hex neighbors of this.Location, computed in Euclidean coords:
       for i in 0..5: result[i] = (x + cos(i·π/3), y + sin(i·π/3))
    The bool flag (= Leftmost) switches between full-6 and a left-half-only
    variant (`enumerateNeighborPositionsLeft`), an optimization to avoid
    emitting symmetric extensions twice.

  distance(node1, node2) -> float64
    Wrapper around LatticePoint.distance for two nodes.

  clone() -> CGNode
```

## CGMotif (one skeleton)

The big one. Fields tell you what's tracked:

```
linked-list view:
  Start, End      : CGNode    head/tail pointers
  Count           : int32     number of nodes

lattice view:
  Lattice         : Dictionary<LatticePoint, CGNode>
                              position → node lookup. THE source of truth
                              for adjacency. Adjacency check = is the
                              neighbor position a key in this dict?

family / tree view:
  parent, child   : int32     parent/child motif IDs (combine history)
  Id              : int32     sequential id within Sn
  dominatingId    : int16     RCC canonical-representative ID
  StringId        : string

RCC scoring (matches paper Appendix §1.2):
  distSum         : float64   sum of pairwise Euclidean distances
  layers          : int16     count of unique-y coordinates (eps 0.01)
  MinX            : float64   (canonical-form bounding hint)

handedness (for equivalence check):
  handSum1        : uint64    hash of handedness across all triples (?)
  handSum2        : uint64    extra hash bits for Count == 7
  handSum3        : uint64    extra hash bits for Count == 8

filters / flags:
  scc             : bool      passes SCC-2 (induced-grid whitelist)
  forHigherDim    : bool      retained for use in larger dim's combine
  helix           : bool      SSE-type flag
  helix6          : uint8     packed helix info
  threec          : bool      three-collinear-points present
  fourc           : bool      four-collinear
  threenum        : uint8     count of 3-collinear triples
  Mark4Delete     : bool      sweep-deletion flag (used by removeDupMotifs)
  surroundNum     : int16

caches (arrays, fixed size 8 — the model's max dim is 5 so 8 is slack):
  _distances      : float64[8,8]  pairwise distance matrix
  _addIA          : int16[8,8]    additional interactions
  _elements       : CGElementType[8]  SSE types
  _surroundPoints : float64[2,34]  surrounding non-skeleton points
  _surrCompact    : byte[34]       PCC scratchpad

canonical position:
  canonicalStart  : LatticePoint  (default (10, 10) constant)
```

## Compactness criteria

### PCC ("primary compactness")

Not a separate field — checked inline during combine. Predicate:

```
for each lattice point p not in skeleton.Lattice:
    if p has 4 or more neighbors that ARE in skeleton.Lattice:
        FAIL
```

i.e., no "hole" in the skeleton interior has too many neighbors. Walks
the `_surroundPoints` cache built during combine.

### SCC-1 ("secondary compactness, criterion 1")

When adding a new node `p` to skeleton `s`:

```
candidates_in_s_adj_to_p = [n for n in s.Lattice.values() if Location.is_neighbor(p, n)]
if len(candidates_in_s_adj_to_p) >= 2:
    PASS
elif len(candidates_in_s_adj_to_p) == 1:
    # Must extend a collinear line in s
    anchor = candidates_in_s_adj_to_p[0]
    for other in s nodes adjacent to anchor:
        if collinear(other, anchor, p):
            PASS
    FAIL
else:
    FAIL  (disconnected)
```

The `enumerateNeighborPositions` + `getNodesAtPosition` machinery
finds the candidate positions and the per-position adjacency count.

### SCC-2 ("secondary compactness, criterion 2")

CG-2012 does **not** store an explicit grid whitelist (Fig. S3 in the
paper). Instead, the `scc` boolean on CGMotif is set based on inline
predicates (`threec`, `fourc`, `threenum`, `helix`, etc.). A motif
passes SCC if its specific structural flags match an allowed pattern.

This is more rigid than our `WHITELIST_S4` / `WHITELIST_S5` approach,
which checks against canonical unlabeled grids derived from Fig. S3.
Both yield the same accept/reject outcome for paper-allowed grids,
but the representations differ.

## Combine driver

The full pipeline at each dimension *n*:

```python
def driver(n):
    sn_motifs = []
    for parent in s_{n-1}_motifs_passed_to_this_dim:
        extensions = parent.extend()  # local RCC dedup happens inside
        sn_motifs.extend(extensions)

    sn_motifs = removeDupMotifs(sn_motifs)  # sweep Mark4Delete-flagged

    # Partition by (scc=true AND forHigherDim=false) — these get emitted
    output_for_this_dim = [m for m in sn_motifs if m.scc and not m.forHigherDim]
    carry_forward     = [m for m in sn_motifs if (not m.scc) or m.forHigherDim]

    # Output: no further RCC. Goes to IA.txt.
    write_ia_txt(output_for_this_dim)

    # Carry-forward: pairwise global RCC dedup before passing to next dim
    carry_forward = global_rcc_dedup(carry_forward)
    pass_to_next_dim(carry_forward)
```

### `CGMotif::extend()`

```python
def extend(self) -> List[CGMotif]:
    if self.Count > 0:
        self.setLeftmost()
    ext_start = self.extendMotifStart()   # front-join extensions
    ext_end   = self.extendMotifEnd()     # end-join extensions
    result = list(ext_start)
    for e in ext_end:
        result.append(e)
        if len(result) > 1:
            eliminateEquivalentMotifs1(result)  # local RCC dedup
    return result
```

### `CGMotif::checkEquivalence5Gr(m1, m2) -> int32`

The RCC equivalence + selection function. Returns:
- `0` if not equivalent (sizes mismatch or handedness hash mismatch)
- `1` if m1 dominates (lower `distSum`, or tied with lower-or-equal `layers`)
- `11` if m2 dominates

Decoded logic:

```python
def checkEquivalence5Gr(m1, m2):
    if m1.Count != m2.Count: return 0
    if m1.Count < 3:         return 1   # trivially equivalent, m1 wins

    result = 0
    if abs(m1.distSum - m2.distSum) > 0.01:
        if m1.distSum >= m2.distSum:  # m2 has lower distSum
            result = 10
        # else m1 wins (result stays 0)
    else:
        # distSums tied → layers tie-break
        if m2.layers < m1.layers:
            result = 10
        # else m1 wins or true tie (result stays 0)

    if not chiralSameMay12(m1, m2):
        return 0  # not equivalent regardless of distSum/layers

    return result + 1   # 1 if m1 wins, 11 if m2 wins
```

The "first skeleton wins by default" rule from paper Appendix §1.2
corresponds to the case where distSums and layers both tie:
`result = 0`, chirality matches, return `0 + 1 = 1` → m1 wins (m1 was
the iteration's first).

### `chiralSameMay12(m1, m2) -> bool` (equivalence predicate)

```python
def chiralSameMay12(m1, m2):
    assert m1.Count == m2.Count   # else throw ApplicationException
    n = m1.Count
    if n < 7:   return m1.handSum1 == m2.handSum1
    if n == 7:  return m1.handSum1 == m2.handSum1 and m1.handSum2 == m2.handSum2
    if n == 8:  return (m1.handSum1 == m2.handSum1
                        and m1.handSum2 == m2.handSum2
                        and m1.handSum3 == m2.handSum3)
    return True   # n>8 never happens in Chitturi enumeration
```

For S1..S5 (the only sizes the paper enumerates), equivalence reduces
to **a single 64-bit hash comparison on `handSum1`**. The hash packs
per-triple handedness signs. Our Python port uses the full
`handedness_signature(skel)` tuple in place of the uint64 hash — finer
than CG-2012's hash but conservative (no false positives from hash
collisions).

## Two-tier RCC application (the key finding)

CG-2012 applies RCC **in two distinct contexts**:

| Context | Where | Effect |
|---|---|---|
| Local (per-parent batch) | `extend()` via `eliminateEquivalentMotifs1` | Dedupe handedness-equivalent siblings before emitting |
| Global (carry-forward only) | `Program::Main` after partitioning | Dedupe the carry-forward set before it feeds the next dim's combine |

The **output stream** (IA.txt records) gets only local-tier dedup.
Siblings of *different* parents that happen to be handedness-equivalent
end up as distinct records in the output. This is why the oracle
records 41 distinct skel_ids at S4 while a global RCC pass would
collapse them to ~23.

Our Python port (`enumerate_skeletons`) doesn't model the
output/carry-forward split. To exactly match CG-2012's per-dimension
counts we'd need to:

1. Track `scc` and `forHigherDim` per emitted Skeleton.
2. Apply local-RCC inside `combine_with_single_node` and
   `combine_two_skeletons` (per-parent's batch of emissions).
3. Apply global RCC only on the (`scc=false` OR `forHigherDim=true`)
   subset that feeds the next dim.

We've decided this isn't worth the refactor for a +1 at S4 — the
paper-Methods-traceable count is what we want to defend in publications,
and exact CG-2012 reproduction risks encoding implementation quirks
(see [[feedback-reproducibility-over-legacy]]).

## Cross-references in our codebase

| CG-2012 concept | Our port | Notes |
|---|---|---|
| `LatticePoint::distance` | `lattice.py` `is_adjacent` + `combine.py` `_axial_to_euclidean` | We use integer axial coords, no float drift |
| `CGNode::enumerateNeighborPositions` | `lattice.py` `LatticePoint.neighbors` | 6 hex neighbors |
| `CGMotif.Lattice` (dict) | `Skeleton.points` (tuple) | We have a dict-lookup equivalent via `is_adjacent` |
| `CGMotif.distSum` | `combine.dist_sum(skel)` | Ported |
| `CGMotif.layers` | `combine.layers(skel)` | Ported |
| `CGMotif.handSum1` | `combine.handedness_signature(skel)` | Full sig instead of uint64 hash |
| `checkEquivalence5Gr` | `combine.rcc_dedup(skels)` | Group + select |
| `eliminateEquivalentMotifs1` | (not separately ported) | Local-RCC isn't applied; see two-tier discussion |
| `CGMotif.scc` | (implicit in `passes_scc_2`) | We check WHITELIST_Sn |
| `CGMotif.forHigherDim` | (not modeled) | Would need to refactor enumerate_skeletons |
| SCC-2 grid whitelist | `grids.py` `WHITELIST_S3/S4/S5` | We store explicit lists; CG-2012 uses inline predicates |

## Open puzzles we did *not* fully trace

1. **How `handSum1` is computed.** We can see it's a uint64 accumulated during construction. Exact hash function (which triples, ordering, packing) wasn't traced; our `handedness_signature` is a conservative substitute.
2. **K1,4 / C5 oracle records.** We have a strong working hypothesis (float-precision artifacts in cos/sin neighbor computation), but didn't trace matrix construction to confirm.
3. **`forHigherDim` setter.** What predicate flips a motif into "retained for higher dim"? Likely related to SCC-2 failure but with a "useful as building block" override.
4. **Exact `eliminateEquivalentMotifs1` ordering.** Whether it scans pairwise or uses a smarter data structure. Affects whether the local RCC is deterministic across runs.

These can be resolved by deeper IL reading. None are blockers for our
current scope (paper-traceable enumeration for ProSMoS queries).
