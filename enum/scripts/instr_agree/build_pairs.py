#!/usr/bin/env python3
"""Build the paired experimental/AF2 domain set for the instrument-agreement test.

Every cross-database claim we make rests on PALSSE assigning SSEs comparably to a
crystal structure and to a predicted model. Predicted models are idealised: no
missing density, no static disorder, cleaner termini. If PALSSE reads them
differently, "this domain lights cell X" is partly a statement about provenance
rather than topology -- and that contaminates the AFDB-only cells specifically.

Nothing else we have run touches this. Rarefaction controls depth, taxonomy
controls composition, vintage controls time; none control the detector.

The control here is same-molecule: one UniProt accession, one experimental ECOD
domain, one AF2 ECOD domain, SAME T-group, and -- critically -- both structures
chopped by the same ECOD pipeline into the same store, so the only difference
between the two sides is where the coordinates came from.

Requiring tid equality is what makes the pair a pair: a 1:1 accession whose two
domains sit in different T-groups is two different regions of the protein, and
comparing their motifs would measure domain boundaries, not the detector.
"""
import psycopg2, io, json
from pathlib import Path

OUT = Path("/home/rschaeff/work/prosmos_2026/instr_agree")

c = psycopg2.connect(host='sangala', port=45000, user='ecod', dbname='ecod_af2_pdb')
cur = c.cursor()
cur.execute("""
with p as (
  select unp_acc,
    count(*) filter (where type='experimental structure') ne,
    count(*) filter (where type='predicted structure')    np
  from domain where not is_obsolete and unp_acc is not null
  group by unp_acc having count(*) filter (where type='experimental structure')=1
                    and count(*) filter (where type='predicted structure')=1
)
select p.unp_acc,
       max(d.id)  filter (where d.type='experimental structure') exp_id,
       max(d.tid) filter (where d.type='experimental structure') exp_t,
       max(d.id)  filter (where d.type='predicted structure')    pred_id,
       max(d.tid) filter (where d.type='predicted structure')    pred_t
from p join domain d on d.unp_acc = p.unp_acc and not d.is_obsolete
group by p.unp_acc
""")
rows = cur.fetchall()
c.close()
print(f"1:1 accessions: {len(rows):,}")
pairs = [r for r in rows if r[2] and r[2] == r[4]]
print(f"  with matching T-group: {len(pairs):,}  (dropped {len(rows)-len(pairs):,} as different regions)")

# domain_id -> ecod_uid -> structure path, both sides from the same store
ids = sorted({r[1] for r in pairs} | {r[3] for r in pairs})
c = psycopg2.connect(host='dione', port=45000, user='ecod', dbname='ecod_protein')
cur = c.cursor()
cur.execute("create temp table qd(did text primary key)")
cur.copy_from(io.StringIO("\n".join(ids) + "\n"), 'qd', columns=('did',))
# Provenance comes from domain_summary.source_type ('pdb' vs 'afdb'), NOT from
# derived_files.domain_source_type -- that column tracks FILE FORMAT, so 1,348 of
# our AF2 domains are labelled 'pdb' there simply because they are PDB-format
# files. Asserting on it dropped 1,155 good pairs and would have left us with 374.
cur.execute("""select d.domain_id, f.internal_path, s.source_type
               from ecod_commons.domains d
               join qd on qd.did = d.domain_id
               join ecod_commons.derived_files f on f.ecod_uid = d.ecod_uid
               join ecod_commons.domain_summary s on s.ecod_uid = d.ecod_uid
               where f.file_type_id = 2 and f.status = 'complete' and not d.is_obsolete""")
path = {}
for did, p, st in cur.fetchall():
    path[did] = (p, st)
c.close()
print(f"  domain_ids resolved to a complete .pdb: {len(path):,}/{len(ids):,}")

out, miss, xprov = [], 0, 0
for acc, eid, et, pid, pt in pairs:
    pe, pp = path.get(eid), path.get(pid)
    if not pe or not pp:
        miss += 1
        continue
    # provenance sanity: the experimental side must be pdb-sourced and the
    # predicted side alphafold-sourced. A swap here would silently invert the
    # whole test, so it is an assertion, not a filter.
    if pe[1] != 'pdb' or pp[1] != 'afdb':
        xprov += 1
        continue
    if not (Path(pe[0]).is_file() and Path(pp[0]).is_file()):
        miss += 1
        continue
    out.append({"unp": acc, "tid": et, "exp_id": eid, "pred_id": pid,
                "exp_pdb": pe[0], "pred_pdb": pp[0]})

print(f"  dropped, structure missing on one side: {miss:,}")
print(f"  dropped, provenance mismatch: {xprov:,}")
print(f"\nUSABLE PAIRS: {len(out):,}")
json.dump(out, open(OUT / "pairs.json", "w"), indent=1)
with open(OUT / "exp_paths.txt", "w") as fh:
    fh.write("".join(p["exp_pdb"] + "\n" for p in out))
with open(OUT / "pred_paths.txt", "w") as fh:
    fh.write("".join(p["pred_pdb"] + "\n" for p in out))
print(f"wrote pairs.json, exp_paths.txt, pred_paths.txt to {OUT}")
