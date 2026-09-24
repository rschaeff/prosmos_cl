# acceptance: the post-ProSMoS acceptance rule

`acceptance.py` turns `searchmatrix` hits into accepted occurrences. It applies the four-clause
rule of the four-strand topology census (Methods, "Candidate matches and acceptance rule"):

| stage | clause | rule |
|---|---|---|
| `hits` | none | read hits, resolve record names, attach residue ranges and the smallest containing sheet |
| `sieve` | 1. disjointness | the four residue ranges share no residue |
| `sieve` | 2. axial overlap | on every interior strand, both flanks overlap it by > 3 Å along its PALSSE axis |
| `pairing` | 3. local pairing | on every interior strand, ≥ 2 consecutive residues whose nearest flank Cα atoms are within 6.0 Å of Cα(i) and on opposite sides of the plane Cα(i−1), Cα(i), Cα(i+1) |
| `bridge` | 4. DSSP bridges | every pair the query marks `c`/`t` carries ≥ 1 DSSP bridge (BP1/BP2) |
| `census` | none | per query: accepted records at sheet size ≤ 4 (complete), ≤ 6 (core) and any size (raw) |

Every stage writes every row back out. The `status` column says what happened to each row:

- `ok` means the row is still in.
- `reject:<clause>` means the rule was applied and said no.
- `fail:<reason>` means the row could not be judged, and it is **not** accepted.

`census` prints the count for every status. A non-zero exit means an input problem, such as a
missing file, a hits path with no hits, a structure path that matches nothing, or a malformed
query. The message on stderr says which one.

`RUNBOOK_GREY.md` in the paper's working tree has the worked example.

## Failure reasons

| status | cause |
|---|---|
| `fail:segment_count` | The hit has fewer segments than the query has elements. `archive_hits.py` dropped segments whose range carries an insertion code. Read the hit tree instead of an archived table to recover them. |
| `fail:record_not_in_db`, `fail:range_mismatch`, `fail:ambiguous_record` | The hit's record name is not in the DB, even after resolving searchmatrix's 19-character truncation against the element ranges. |
| `fail:malformed_record` | The metamatrix record has overflowed fixed-width fields. |
| `fail:slot_mismatch` | The record's codes disagree with the query on a required pair. The hit's element order is not the query's slot order. Never seen; this is a guard. |
| `fail:no_sheet` | No PALSSE sheet holds all four elements. |
| `fail:no_structure`, `fail:no_ca_atoms` | The structure file is missing or has no Cα atoms. |
| `fail:dssp_ca_only` | DSSP cannot assign a Cα-only model. |
| `fail:dssp_no_carbonyl` | The model has only N, Cα and C, so there is no H-bond acceptor. |
| `fail:dssp_atom_names` | The atom names are left-justified in columns 13–16, and DSSP rejects every residue. |
| `fail:dssp_no_output` | DSSP assigned nothing, for another reason. |

## Agreement with the paper's census code

Implementations compared: the paper's intermediate files, produced by `ruczinski/bridge_prep.py`
with the sign clause off, `nvg_worker.py`, `bridge_worker.py`/`bridge_rules.py` and
`census_plane.py`.

Test: 1,400 experimental domains (1,200 random, plus 200 that DSSP cannot read), 5,166
candidates, compared instance by instance (2026-09-23).

- Clauses 1–2: all 966 rejects and 4,180 passes agree.
- Clause 3: 180 of 181 rejects agree. The exception is a domain with Kabat insertion codes
  (`52A`, `100A`–`100D`). The old code keyed Cα atoms by residue number only and collapsed
  those residues together. This version keys on insertion codes too.
- Clause 4: all 73 rejects and 3,391 passes agree.

There are three deliberate differences from the paper's code:

1. **No DSSP verdict now means not accepted.** `census_plane.py` accepted such instances without
   checking bridges. In the experimental archive that is 11,239 instances over 4,089 domains.
   The census occupancy (73 / 86 / 93) is unchanged. Counts fall by about 0.5%.
2. **Insertion codes are handled.** Residues are keyed (chain, number, insertion code).
   Hit-tree ranges like `82A` are read, not dropped.
3. **Nothing is dropped silently.** The old path skipped truncated rows, malformed records and
   degenerate axes without a trace. They are now `fail:*` rows.

Residue order within a range is the order of (number, insertion code). A structure that
numbers insertions out of alphabetical order (`82B` before `82A`) can therefore be misjudged
by disjointness. This is rare, and the stage does not read the structure.
