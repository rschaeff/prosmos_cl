# archive/ — exploratory work, not part of the published census

Everything under this directory is retained for provenance and reproducibility of the
project's history. **None of it backs any claim in the four-strand β-sheet census paper**
(Schaeffer, Guo, Cong & Grishin). It is kept because deleting work that informed a decision
makes the decision harder to audit later, not because it supports a published result.

If you arrived here from the paper's code citation, the material you want is in the
repository root (`searchMatrix/`, `generateMatrix/`) and in the Zenodo deposition — not here.

## `s5_negspace/` — the five-element lattice enumeration

An enumeration of five-SSE arrangements on a hexagonal lattice, following the model of
Chitturi et al. (2016), and a survey of which enumerated cells go unobserved in the PDB and
the AlphaFold database ("negative space").

**Cut from the paper on 2026-07-30.** The census that was published is confined to the 96
four-strand β-sheet topologies enumerated by Ruczinski et al. (2002) — a small, exactly
specified state space whose complete invariant is strand order and per-strand orientation.
The five-element lattice enumeration is a different and much larger claim about protein
architecture in general, it was not carried through the geometric-validity work that the
four-strand census rests on, and its numbers were superseded. It should not be read as a
result, and in particular:

- `docs/paper_section.md` is a **draft that was never published**. It describes methods and
  results for the S5 negative-space analysis and is retained only as a record of that
  direction. Nothing in it was included in the four-strand census paper.
- The 800 `queries/s5-*.query` files are the typed query set for that enumeration.
- The `scc2_*` and `cg2012_*` documents concern the lattice model's self-complementarity
  constraint and its relationship to the Chitturi enumeration.

The published paper is explicit that ProSMoS is used here as an implementation of an
*already-defined* state space, and that it makes no claim about ProSMoS's adequacy as a
general representation of protein fold space. This directory is precisely the material that
claim is scoped away from.
