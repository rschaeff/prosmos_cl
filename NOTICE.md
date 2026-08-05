# Provenance and attribution

**SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0**

This repository is a **derivative work**. It is not original software, and the notices below
travel with every copy of it under the `Required Notice:` lines in [`LICENSE.md`](./LICENSE.md).

## Upstream

`searchMatrix/` and `generateMatrix/` descend from the ProSMoS release by Shuoyong Shi
(UT Southwestern, ~2010), preserved verbatim in [`readme.original`](./readme.original).

> Shi S, Zhong Y, Majumdar I, Sri Krishna S, Grishin NV. Searching for three-dimensional
> secondary structural patterns in proteins with ProSMoS. *Bioinformatics*
> 2007;23(11):1331–1338.
>
> Shi S, Chitturi B, Grishin NV. ProSMoS server. *Nucleic Acids Res*
> 2009;37(Web Server):W526–W531.

The original terms, quoted in full from `readme.original`:

> ProSMoS is free for academic use only. For non-academic use, please contact Shuoyong Shi.

`readme.original` is retained permanently as the record of those terms and must not be removed.

## What this tree adds

Defect fixes, a build system, a hardened meta-matrix database reader, and loop inversion in the
search path. These are modifications *of* the upstream code, not a replacement for it. The
copyright in the underlying work is unchanged; this repository claims only the modifications.

## Why PolyForm Noncommercial

"Free for academic use only" is a statement of intent, not a license. It grants no explicit
right to redistribute, modify or sublicense, says nothing about warranty, and cannot be
expressed as a machine-readable identifier — journals, repositories and package tooling all
require a named license.

PolyForm Noncommercial 1.0.0 is the closest faithful rendering of that intent as an actual
license. Its *Noncommercial Organizations* clause permits use by "any charitable organization,
educational institution, public research organization, public safety or health organization,
environmental protection organization, or government institution ... **regardless of the source
of funding or obligations resulting from the funding**" — which resolves the ambiguity that
makes a bare "academic use" restriction hard to rely on for a university lab holding industry
funding. It grants the redistribution and modification rights the original text omits, and it
carries an SPDX identifier (`PolyForm-Noncommercial-1.0.0`).

It is deliberately **not** an OSI-approved license. No license restricting commercial use can be:
that restriction fails clause 6 of the Open Source Definition. This repository is
*source-available*, not open source, and should not be described as the latter.

## Basis for this license selection

Recorded 2026-08-05, so that the reasoning is documented rather than reconstructed later.

The original terms grant *use* and are silent on *redistribution*. No further clarification of
authorial intent is obtainable: the original author left the laboratory approximately fifteen
years ago, and the 2010 release is the last statement of terms that exists. The selection below
is therefore made by the Grishin Laboratory as the originating group, on three grounds.

**1. The copyright almost certainly never vested in the author personally.** ProSMoS was written
at UT Southwestern Medical Center in the course of employment. Under the ordinary rule for
institutionally-authored software, the copyright is the institution's. If that is so, the
laboratory is not acting as anyone's proxy — it is acting within the same institution that holds
the right, and an individual author's blessing was never the operative permission. The one party
who could speak with more authority than the laboratory is UT Southwestern's technology
office, not the original author.

**2. The selection grants no more than the original offered.** This is the load-bearing point.
"Free for academic use only" and PolyForm Noncommercial permit the same class of user for the
same class of purpose; PolyForm adds only the redistribution and modification mechanics without
which the software cannot lawfully be released at all, plus a warranty disclaimer that protects
the original authors rather than exposing them. A permissive license — MIT, BSD — would have
extended rights the original expressly withheld, and was rejected for that reason. Acting on an
unobtainable permission is low-risk precisely because the action does not expand anyone's rights.

**3. Attribution is preserved and strengthened.** The original terms are retained verbatim in
`readme.original`; the authorship is carried in `Required Notice:` lines that PolyForm obliges
every downstream recipient to pass on. The 2010 release carried no such obligation.

**What this is not.** It is not a claim that the original author consented, and it should not be
described as one. It is a good-faith selection, on the reasoning above, by the group in the best
position to make it.

## Third-party code

`scripts/map2scop/mapscop.pl` ships upstream with a placeholder credential line
(`my $password="yourpassword"`). Harmless, but it trips automated secret scanners; replace it
with a `~/.my.cnf` read or drop the script before any public release.

## Related deposition

The *analysis* built on this engine — the census, controls and figures for the four-strand
β-sheet topology census — is deposited separately and contains none of this code. It is
released under **CC BY 4.0** and does not inherit the restriction above.
