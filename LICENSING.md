# Licensing — decided

**Status as of 2026-08-05: decided and in force.** The selection is made by the Grishin
Laboratory as the originating group; the reasoning and its limits are recorded in
[`NOTICE.md`](./NOTICE.md) and should be read as part of this decision. This file is the summary
and the residual action list.

## Decision

| Artifact | License | Status |
|---|---|---|
| This repository (`searchMatrix`, `generateMatrix`) | **PolyForm Noncommercial 1.0.0** | in force |
| `palsse_cl` | **PolyForm Noncommercial 1.0.0** | in force |
| Zenodo analysis deposition (10.5281/zenodo.21814885) | **CC BY 4.0** | in force |

The split is deliberate. The engines are derivative works and inherit a restriction. The
deposition contains no upstream code — it is the *output* of running these tools, not a
modification of them — so it inherits nothing and is released openly. Reasoning in `NOTICE.md`.

Rejected: a hand-written "academic use" license (do not freelance license text when a
professionally drafted one exists); the **Academic Free License 3.0**, which despite its name is
permissive and carries no academic restriction at all; and Prosperity 3.0.0, which has no SPDX
identifier and so cannot be expressed to GitHub, Zenodo or package tooling.

## Why no upstream sign-off was sought

Deliberately, not by oversight. The original ProSMoS author left the laboratory roughly fifteen
years ago; the 2010 release is the final statement of terms in existence, and no clarification of
intent is obtainable. The PALSSE terms are not merely unclear but absent — no copy survives. There
is no historical reach here that would produce firmer ground than the laboratory's own position.

That position is stronger than a proxy's, for the reason set out in `NOTICE.md`: copyright in
software written at UT Southwestern in the course of employment vests in the institution, not in
the individual author, so the laboratory is acting inside the entity that holds the right rather
than standing in for someone outside it.

The selection is conservative by design. PolyForm Noncommercial permits the same users for the
same purposes as "free for academic use only" and adds only the redistribution mechanics required
to release the software at all. A permissive license would have granted rights the original
expressly withheld. Acting without an obtainable permission carries little risk when the action
does not expand anyone's rights — and that, rather than anyone's consent, is what makes this
defensible.

## Residual items

1. **Confirm the copyright holder line.** `LICENSE.md` names UT Southwestern Medical Center on the
   assumption above. If the lab has a different arrangement, correct the `Required Notice:` lines.
   If certainty is ever wanted, this is a single question to UTSW's technology office — not an
   archaeology problem.
2. **Clean `scripts/map2scop/mapscop.pl`** — the placeholder credential line trips secret scanners
   on public repositories.
3. **Retain `readme.original` permanently** — it is the only surviving record of the original
   terms, and the `Required Notice:` attribution depends on it.

## Standing caution

Do not describe this repository as open source. A noncommercial license is source-available by
definition and can never be OSI-approved, because restricting commercial use fails clause 6 of the
Open Source Definition. The Zenodo deposition is a different matter and *is* openly licensed.
