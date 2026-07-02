#!/usr/bin/env python3
"""Select a curated set of cells (unitypical exemplars + promiscuous montages)
and gather everything needed to render them: for each exemplar domain, its
pdb/chain/domain-range (ecod_rep) and the 5 matched SSE ranges (hit file).

Emits render_spec.json consumed by render_pymol.py.
"""
import json, re, os
from pathlib import Path
import psycopg2

SP = Path("/tmp/claude-1219/-home-rschaeff-dev-prosmos-cl/02f20625-920b-46f1-a2bf-bc06d84727af/scratchpad")
HITROOT = Path.home() / "work/prosmos_2026/ecod_search_v4/hits"
HITROOT2 = Path.home() / "work/prosmos_2026/ecod_search_v4_retry/hits"
D = json.load((SP / "s5_promiscuity.json").open())
row = D["rowSkeleton"]


def tystr(ty):
    return "".join("E" if (ty >> b) & 1 else "H" for b in range(4, -1, -1))


def read_segments(sk, ty, euid):
    """Return list of (type,start,end,chain) from first MOTIF block of the hit file."""
    fn = f"pdb{euid:09d}.txt"
    for rootdir in (HITROOT, HITROOT2):
        p = rootdir / f"s5-{sk:04d}-{ty:04d}" / fn
        if p.exists():
            segs = []
            for line in p.read_text().splitlines():
                m = re.match(r"segment-Type:\s*(\S+)\s+Position:\s*\d+\s+Range:\s*(\d+)\s*--\s*(\d+)\s+(\S+)", line)
                if m:
                    segs.append((m.group(1), int(m.group(2)), int(m.group(3)), m.group(4)))
                elif line.strip() == "END" and segs:
                    break
            return segs[:5]
    return None


# ---- choose curated cells --------------------------------------------------
cells = {}
for k, v in D["cellGroups"].items():
    r, ty = map(int, k.split(","))
    cells[k] = (r, ty, len(v["T"]), v["nhit"], v)

# unitypical: nT==1, decent depth, diverse typings (by strand count)
uni = [c for c in cells.values() if c[2] == 1 and c[3] >= 4]
uni.sort(key=lambda c: -c[3])
picked_uni, seen_ty = [], set()
for c in uni:
    scount = bin(c[1]).count("1")
    if scount in seen_ty:
        continue
    # need an experimental exemplar
    g = c[4]["T"][0]
    if not g[2].startswith("e"):
        continue
    picked_uni.append(c); seen_ty.add(scount)
    if len(picked_uni) == 5:
        break

# promiscuous: large nT, one per DISTINCT skeleton and distinct typing pattern
pro = [c for c in cells.values() if c[2] >= 40]
pro.sort(key=lambda c: -c[2])
picked_pro, seen_sk, seen_ty = [], set(), set()
for c in pro:
    sk = row[c[0]]
    scount = bin(c[1]).count("1")
    if sk in seen_sk or scount in seen_ty:
        continue
    picked_pro.append(c); seen_sk.add(sk); seen_ty.add(scount)
    if len(picked_pro) == 3:
        break

# ---- gather domain info from ecod_rep --------------------------------------
need_dids = set()
for c in picked_uni:
    need_dids.add(c[4]["T"][0][2])
for c in picked_pro:
    for g in c[4]["T"][:8]:
        if g[2].startswith("e"):
            need_dids.add(g[2])

conn = psycopg2.connect(host="dione", port=45000, dbname="ecod_protein", user="ecod")
cur = conn.cursor()
cur.execute("select ecod_domain_id,chain_id,pdb_range,seqid_range from ecod_rep.domain "
            "where ecod_domain_id = any(%s)", (list(need_dids),))
dinfo = {}
for did, chain, pr, sr in cur.fetchall():
    dinfo[did] = {"chain": chain, "pdb_range": pr, "seqid_range": sr}

def dom_spec(sk, ty, group):
    gid, cnt, did, euid = group
    pdb = did[1:5]
    info = dinfo.get(did, {})
    segs = read_segments(sk, ty, euid)
    return {"did": did, "pdb": pdb, "group": gid, "count": cnt,
            "chain": info.get("chain"), "pdb_range": info.get("pdb_range"),
            "segments": segs}

jobs = {"unitypical": [], "promiscuous": []}
for c in picked_uni:
    r, ty, nt, nhit, v = c
    sk = row[r]
    jobs["unitypical"].append({
        "sk": sk, "ty": ty, "typing": tystr(ty), "nT": nt, "nhit": nhit,
        "exemplar": dom_spec(sk, ty, v["T"][0]),
    })
for c in picked_pro:
    r, ty, nt, nhit, v = c
    sk = row[r]
    ex = []
    for g in v["T"]:
        if len(ex) == 6:
            break
        if not g[2].startswith("e"):
            continue
        d = dom_spec(sk, ty, g)
        if d["segments"] and d["pdb_range"]:
            ex.append(d)
    jobs["promiscuous"].append({
        "sk": sk, "ty": ty, "typing": tystr(ty), "nT": nt, "nhit": nhit,
        "exemplars": ex,
    })

(SP / "render_spec.json").write_text(json.dumps(jobs, indent=1))
print("unitypical picks:")
for j in jobs["unitypical"]:
    e = j["exemplar"]
    print(f"  sk{j['sk']:04d} ty{j['ty']:02d} {j['typing']}  T={e['group']}  {e['did']}  "
          f"segs={'ok' if e['segments'] else 'MISSING'} range={e['pdb_range']}")
print("promiscuous picks:")
for j in jobs["promiscuous"]:
    print(f"  sk{j['sk']:04d} ty{j['ty']:02d} {j['typing']}  nT={j['nT']}  "
          f"montage of {len(j['exemplars'])} exemplars: "
          + ", ".join(e['did'] for e in j['exemplars']))
