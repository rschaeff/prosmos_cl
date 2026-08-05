# Licensing — unresolved, needs sign-off before public release

**Status: this repository has no LICENSE file, and one should not be added without a decision
from the original authors.** This note records the situation so the decision can be made
deliberately rather than by default.

## Why it is not straightforward

`searchMatrix/` and `generateMatrix/` are **derivative works**, not original code. They descend
from the ProSMoS release by Shuoyong Shi (UT Southwestern, ~2010), whose own terms are stated
in [`readme.original`](./readme.original) as *"free for academic use"* — a phrase that is not a
license. It grants no explicit rights to redistribute, modify, or sublicense, says nothing
about commercial use, and carries no warranty disclaimer. Journals and Zenodo both expect a
named license, and "free for academic use" will not satisfy either.

This tree adds substantial modifications (defect fixes, a build system, a hardened database
reader, loop inversion) but those are modifications *of* that code, so the downstream license
cannot be chosen unilaterally.

## What needs to happen

1. **Confirm the upstream position with Shuoyong Shi and Nick Grishin.** They are the parties
   who can say what the original terms were intended to be and whether the lab is willing to
   relicense.
2. **Pick a named license** for the combined work. If the goal is maximum reuse and the
   upstream authors agree, MIT or BSD-3-Clause are the obvious candidates and are what most
   structural-bioinformatics tooling uses. If the intent is genuinely to restrict commercial
   use — which "free for academic use" gestures at — then a custom academic license is needed,
   and it should be written explicitly rather than implied; note that non-commercial terms are
   not OSI-approved and some funders and repositories treat them as non-open.
3. **Add `LICENSE` at the repository root**, and state the upstream provenance in it so the
   derivation is not lost.
4. **Mirror the choice in the Zenodo deposition**, which asks for a license at upload time.

## Third-party code to check while you are at it

- `scripts/map2scop/mapscop.pl` ships upstream with a placeholder credential line
  (`my $password="yourpassword"`). Harmless, but it trips automated secret scanners on public
  repositories; consider replacing it with a `~/.my.cnf` read or removing the script if it is
  no longer used.
- `readme.original` should be retained regardless of the outcome — it is the record of the
  original terms.

## Until this is resolved

Do not present the repository as open source, and do not select a license on the Zenodo
deposition that is inconsistent with whatever the upstream authors agree to.
