# Methods and Results — S5 negative space in PDB and AFDB

Draft for inclusion in the 200M paper. Methods + Results sections only;
omits Introduction and Discussion.

---

## Methods

### Enumeration of small SSE arrangements

We re-implemented the lattice enumeration of secondary-structure
patterns (SSPs) described by Chitturi et al. (2016), which represents
each SSP as an induced subgraph of the 2D hexagonal lattice annotated
with per-vertex up/down orientation, subject to a strict size-2
self-complementarity constraint (SCC-2) and a compactness filter. For
SSP size N ∈ {3, 4, 5} the enumeration produces 5, 42, and 198
distinct skeletons respectively. Each skeleton was crossed with the
full 2ᴺ Cartesian product of helix (H) and strand (E) assignments,
yielding 40, 672, and 6,336 typed queries at S3, S4, and S5
respectively (7,048 total). For each typed query we computed a
handedness signature — the sign of the scalar triple product (p×q)·r
for every ordered SSE triple — and emitted a ProSMoS-format query
matrix (Shi et al., 2007) encoding SSE types, adjacency, and
per-triple chirality directives.

### Construction of search databases

Three metamatricesDBs were built using the ProSMoS `generateMatrix`
pipeline applied to PALSSE secondary-structure assignments (Majumdar
et al., 2005). For PALSSE we used a Python 3 re-implementation
(palsse_cl) verified byte-for-byte against the original Python 2
oracle on a 30-domain validation set.

- **Manual-rep DB (v4)** — 19,015 entries, all experimental crystal
  structures, derived from `ecod_rep.domain WHERE manual_rep = true`
  in ECOD develop281.
- **Experimental-PDB DB** — 496,359 entries, built from all ECOD
  domains with `domain_summary.source_type = 'pdb'` (i.e., backed by
  a PDB crystallographic or NMR structure, excluding AFDB-predicted
  domains stored as `.pdb` files).
- **AFDB-non-singleton DB** — 4,921,931 entries, built from the
  non-singleton-cluster representative subset of AFDB v4 (~4.9M
  domains derived from sequences with at least one cluster neighbor
  at the ECOD-curated clustering threshold).

Each build produced a metamatricesDB file via per-chunk PALSSE +
generateMatrix processing under SLURM, followed by post-filter to
remove malformed entries (those containing coordinate fields wider
than the fixed-width format permits). Final entry counts above are
post-filter.

### Identification of negative-space queries

The 7,048-query corpus was first swept against the manual-rep DB
using `searchmatrix` (Shi et al., 2007). Of the 198 S5 skeletons,
**25 returned zero hits across all 32 H/E typings**, defining a
candidate set of 800 motifs absent from the curated experimental
core of protein structure space. These 800 queries — the "negspace
set" — were then re-swept against the larger experimental-PDB DB and
the AFDB-non-singleton DB to test whether the absence held at
broader sampling.

### Bug fix to searchmatrix

During the broader sweeps we identified two buffer-overflow bugs in
the released ProSMoS `searchmatrix` source that produced silent
SIGSEGV and `std::out_of_range` crashes when the input
metamatricesDB contained malformed blocks (orphan matrix lines with
no preceding header, an artifact of `generateMatrix` output
truncation on rare inputs): a 31-byte `strcpy` into a 20-byte `pid`
buffer in `intMnumofele`, and an unbounded matrix-fill loop in
`getInterActionM`. We patched both — bounds-checking the buffers and
adding header validation — and verified that the fixed binary
completes a full 496k-entry DB scan in ~200 s per query without
crashes. All sweeps reported below use the patched binary; source
patches are released alongside this paper.

---

## Results

### S5 motif space contains topologies absent from curated PDB

The full 7,048-query S3–S5 sweep against the v4 manual-rep DB
returned hits for 98%, 88%, and 37% of S3, S4, and S5 queries
respectively. The dropping S5 hit rate reflects two distinct
mechanisms: a composition gradient in which each additional
β-strand in the H/E typing roughly halves coverage (S5 all-helix
73% → S5 all-strand 7%), and a skeleton-level gap in which 25 of
198 S5 skeletons (12.6%) return zero hits across **all** 32 H/E
typings, indicating topology classes entirely absent from the
curated set. The 25 skeletons span the range of S5 graph
densities and handedness profiles; they are not concentrated at
any single graph-theoretic feature.

### The absence persists in the broader experimental PDB

To test whether the manual-rep gap was a curation artifact, we
re-swept the 800 negative-space queries (25 zero-hit skeletons ×
32 typings) against the 496k-entry experimental-PDB DB —
including every PDB-derived domain in ECOD develop281, not only
the manually selected representatives. **All 800 queries returned
zero hits.** Mean per-query runtime against this DB was 192 s
(full linear scan of all 496,359 entries). The 26×-larger
experimental sample produces no new matches.

### The absence persists in AFDB-predicted protein space

The same 800 queries were swept against the AFDB-non-singleton DB
(4,921,931 entries, AlphaFold predictions for the clustered
fraction of UniProt). **All 800 queries returned zero hits.** Mean
per-query runtime 462 s. AFDB extends the sample by a further
~10× (real + predicted = 5.4M structures total) and produces no
matches. Combined per-quadrant tally vs experimental PDB is
`truly_absent = 800`, `AFDB_only = 0`, `AFDB_loss = 0`,
`common = 0`.

### Interpretation

The 800 S5 motifs identified by the lattice enumeration are
**absent from 5.4 million observed or predicted protein structures**
spanning experimental crystallography (curated and exhaustive) and
the largest extant predicted-structure database. Two readings are
compatible with this result. (i) These topologies are unrealized in
nature — either because they are sterically inaccessible under
protein chemistry, or because the evolutionary process has not
sampled them. (ii) AFDB's predictor inherits PDB's structural
distribution from its training set, so the agreement between PDB
and AFDB absence is partly a training-prior artifact and the
underlying claim is "absent from PDB and from PDB-like
predictions." Distinguishing these two readings requires
constructive geometric tests — attempting to build 3D backbones
satisfying the topology and handedness constraints with realistic
SSE geometry — and falls outside the scope of this paper.
