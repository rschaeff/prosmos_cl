"""Cross-DB Ig topology comparison WITHOUT template fitting.

Question: does a covering set of PDB Ig still leave AFDB Ig domains dark?

A ProSMoS motif matches ANY 5 SSEs in N->C order, so the natural unit is the
5-SSE subset, not the domain. For every Ig domain in each DB we enumerate every
connected, two-sheet 5-SSE subset and reduce it to a topology SIGNATURE:

    (canonical sheet-partition, 10-bit adjacency of the 5 chosen SSEs)

An AFDB Ig domain is "dark to the PDB vocabulary" iff NONE of its subset
signatures occurs in ANY PDB Ig domain -- which is exactly the condition under
which no PDB-derived template could ever hit it.

This avoids three defects of the template-fitting design:
  * no circularity (nothing is fitted, so "derived from AFDB => covers AFDB" cannot arise)
  * uses the whole population, not merged_tmpl.py's `n_pass == 5` rule, which keeps
    5.5% of AFDB Ig but only 0.8% of PDB Ig
  * lets us condition on n_pass, so "PDB Ig domains are bigger" (median 9 vs 7)
    cannot masquerade as "AFDB Ig is novel"

NOTE the SSE regex takes ANY alnum chain. PDB records carry chains A/B/C/D/X/1/H;
merged_tmpl.py's `([HE])A` silently drops every non-A SSE (harmless on AFDB, which
is chain-A only; fatal on PDB).
"""
import re, glob, pickle, sys
from collections import Counter, defaultdict, deque
from itertools import combinations

NAME = re.compile(r'^(\S+)\.ssd\s')
SSEC = re.compile(r'([HE])([A-Za-z0-9])\s*(-?\d+[A-Za-z]?)\s*--\s*(-?\d+[A-Za-z]?)\s+(\d+)\s+-?\d+\.\d')
SHEET = re.compile(r'^\s*sheet\s+(\d+)\s+(\d+)\s+(.*)')
MIN = {'H': 8, 'E': 5}
CAP = 14                      # C(14,5)=2002; larger domains are excluded and counted


def adjacency(m, n):
    segs = m.split('*')[1:]
    adj = [set() for _ in range(n)]
    for i in range(n):
        if i >= len(segs):
            break
        for p, ch in enumerate(segs[i]):
            if ch != '-':
                j = i + 1 + p
                if j < n:
                    adj[i].add(j); adj[j].add(i)
    return adj


def connected(sub, adj):
    s = set(sub); seen = {sub[0]}; q = deque([sub[0]])
    while q:
        u = q.popleft()
        for v in adj[u]:
            if v in s and v not in seen:
                seen.add(v); q.append(v)
    return len(seen) == 5


def sigs_of(chunkglob, keep, tag):
    """-> (sig -> n_domains), (domain -> set(sig)), n_pass histogram, skipped"""
    sig_dom = defaultdict(set)
    dom_sig = {}
    npass = Counter(); toobig = 0; ndom = 0
    lines = []

    def flush(lines):
        nonlocal toobig, ndom
        if not lines:
            return
        m = NAME.match(lines[0])
        if not m or m.group(1) not in keep:
            return
        nm = m.group(1)
        toks = SSEC.findall(lines[0])
        sses = [(t[0], int(t[4])) for t in toks]
        pi = [i for i, (ty, l) in enumerate(sses) if l >= MIN[ty]]
        npass[len(pi)] += 1
        if len(pi) < 5:
            return
        ndom += 1
        if len(pi) > CAP:
            toobig += 1
            return
        mat = None; sm = {}
        for ln in lines[1:]:
            if ln.startswith('*'):
                mat = ln.strip()
            else:
                s = SHEET.match(ln)
                if s:
                    sid = int(s.group(1))
                    for x in s.group(3).split():
                        sm[int(x)] = sid
        if mat is None:
            return
        adj = adjacency(mat, len(sses))
        out = set()
        for sub in combinations(pi, 5):
            sids = [sm.get(i + 1) for i in sub]
            if len({x for x in sids if x}) != 2:      # two-sheet subsets only
                continue
            if not connected(sub, adj):
                continue
            order = {}; part = []
            for s in sids:
                if s not in order:
                    order[s] = len(order)
                part.append(order[s])
            adjb = tuple(1 if sub[b] in adj[sub[a]] else 0
                         for a in range(5) for b in range(a + 1, 5))
            out.add((tuple(part), adjb))
        if out:
            dom_sig[nm] = out
            for s in out:
                sig_dom[s].add(nm)

    for ch in sorted(glob.glob(chunkglob)):
        for line in open(ch, errors='replace'):
            if NAME.match(line):
                flush(lines); lines = [line.rstrip('\n')]
            else:
                lines.append(line.rstrip('\n'))
        flush(lines); lines = []
    print(f"{tag}: {ndom:,} Ig domains with >=5 passing SSEs | "
          f"{len(dom_sig):,} have >=1 two-sheet 5-subset | "
          f"{len(sig_dom):,} distinct signatures | {toobig:,} skipped (n_pass>{CAP})", flush=True)
    return sig_dom, dom_sig, npass


if __name__ == '__main__':
    import psycopg2
    c = psycopg2.connect(host='lotta', port=45000, user='ecod', dbname='ecod_protein')
    cur = c.cursor()
    cur.execute("SELECT cluster_name FROM afdb_200m.ecod_cluster_summary WHERE t_group_id LIKE '11.%';")
    afdb_ig = {r[0] for r in cur.fetchall()}
    c.close()
    pdb_ig = {f"{int(l.split(chr(9))[0]):09d}" for l in
              open('/home/rschaeff/work/prosmos_2026/pdb_exp_build/uid_class.tsv')
              if l.split(chr(9))[2].startswith('11.')}
    A = sigs_of('/home/rschaeff/work/prosmos_2026/s5_inv/chunks/chunk_*.db', afdb_ig, 'AFDB')
    P = sigs_of('/home/rschaeff/work/prosmos_2026/s5_pdb_inv/chunks/chunk_*.db', pdb_ig, 'PDB ')
    pickle.dump((A, P), open('/home/rschaeff/work/prosmos_2026/s5_grid/igcross/sigs.pkl', 'wb'))
    print('wrote igcross/sigs.pkl')
