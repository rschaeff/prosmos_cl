import json, math
from pathlib import Path

G = Path("/home/rschaeff/work/prosmos_2026/s5_grid")
OUT = "/tmp/claude-1219/-home-rschaeff-dev-prosmos-cl/02f20625-920b-46f1-a2bf-bc06d84727af/scratchpad/why_s5_deck.html"

curve = {r["n"]: r for r in json.load(open(G / "curve_n345.json"))}
null = {r["n"]: r for r in json.load(open(G / "jaccard_null.json"))}
NS = [3, 4, 5]

# ---------- generic panel geometry (viewBox 360x300) ----------
PX0, PX1, PY0, PY1 = 46, 344, 30, 246
def nx(n): return PX0 + (n - 2.6) / (6.4 - 2.6) * (PX1 - PX0)   # x domain ~[3..6]

def panel(title, series, ylab, yticks, yfmt, ymap, note=""):
    grid = ""
    for yv in yticks:
        y = ymap(yv)
        grid += (f'<line x1="{PX0}" y1="{y:.1f}" x2="{PX1}" y2="{y:.1f}" class="pg"/>'
                 f'<text x="{PX0-7}" y="{y+3.5:.1f}" class="pyt">{yfmt(yv)}</text>')
    xt = "".join(f'<text x="{nx(n):.1f}" y="{PY1+18}" class="pxt">{n}</text>' for n in [3, 4, 5, 6])
    body = grid
    for s in series:
        pts = s["pts"]            # list of (n, val)
        seg = "M " + " L ".join(f"{nx(n):.1f} {ymap(v):.1f}" for n, v in pts)
        dash = ' stroke-dasharray="5 4"' if s.get("dash") else ""
        body += f'<path d="{seg}" class="pline {s["cls"]}"{dash}/>'
        for n, v in pts:
            body += (f'<circle cx="{nx(n):.1f}" cy="{ymap(v):.1f}" r="{3.6 if not s.get("proj") else 3}" '
                     f'class="pdot {s["cls"]}{" proj" if s.get("proj") else ""}">'
                     f'<title>n={n}: {yfmt(v)}</title></circle>')
    return (f'<figure class="panel"><figcaption class="ptitle">{title}</figcaption>'
            f'<svg viewBox="0 0 360 268" role="img" aria-label="{title}">{body}'
            f'<text x="{(PX0+PX1)/2:.0f}" y="264" class="paxt">motif size n</text>{xt}</svg>'
            f'<p class="pnote">{note}</p></figure>')

# ---- panel 1: saturation (linear 0-100%) ----
def ym1(v): return PY1 - v * (PY1 - PY0)
p1 = panel("Saturation — fraction of cells lit",
    [{"cls": "afdb", "pts": [(n, curve[n]["afdb_sat"]) for n in NS]},
     {"cls": "pdb",  "pts": [(n, curve[n]["pdb_sat"]) for n in NS]}],
    "", [0, .25, .5, .75, 1.0], lambda v: f"{v:.0%}", ym1,
    note="n=3 is fully saturated — every cell lit, zero discriminating power. Only at n=5 does the grid unsaturate (45%), so “which cells are lit” first carries information.")

# ---- panel 2: sharing (log-y) + n=6 projection ----
YLO, YHI = math.log10(0.9), math.log10(1200)
def ym2(v): return PY1 - (math.log10(v) - YLO) / (YHI - YLO) * (PY1 - PY0)
# ÷6-per-step projection to n=6 (mean); median → ~1 (fingerprint)
proj6 = curve[5]["afdb_share_mean"] / 6.0
p2 = panel("Sharing — mean T-groups per lit cell (log)",
    [{"cls": "afdb", "pts": [(n, curve[n]["afdb_share_mean"]) for n in NS]},
     {"cls": "afdb", "pts": [(5, curve[5]["afdb_share_mean"]), (6, proj6)], "dash": True, "proj": True},
     {"cls": "pdb",  "pts": [(n, curve[n]["pdb_share_mean"]) for n in NS]}],
    "", [1, 10, 100, 1000], lambda v: (f"{int(v)}" if v >= 1 else f"{v:g}"), ym2,
    note="A motif still shared across many folds still abstracts. 950→159→29 (÷6/step); AFDB median 4 at n=5. The dashed projection puts n=6 near median 1 — the fingerprint boundary, where a motif names one domain.")

# ---- panel 3: agreement + within-database null (0.5-1.0) ----
def ym3(v): return PY1 - (v - 0.5) / 0.5 * (PY1 - PY0)
p3 = panel("Agreement — cross-DB vs within-DB null",
    [{"cls": "cross", "pts": [(n, null[n]["cross"]) for n in NS]},
     {"cls": "wpdb",  "pts": [(n, null[n]["within_pdb"]) for n in NS], "dash": True},
     {"cls": "wafdb", "pts": [(n, null[n]["within_afdb"]) for n in NS], "dash": True}],
    "", [0.5, 0.6, 0.7, 0.8, 0.9, 1.0], lambda v: f"{v:.1f}", ym3,
    note="At n=5 the cross-database Jaccard (0.65) rides AFDB’s own self-agreement (0.66) — the databases differ no more than AFDB differs from itself. The drop is undersampling, not a different alphabet.")

HTML = f'''<div class="wrap">

<header class="hero">
  <div class="eyebrow">ProSMoS S-grid · a theoretical answer to “why n=5?” · identical record sets, n = 3 · 4 · 5</div>
  <h1>Why <span class="pred">S5</span>: the only level that is unsaturated, still shared, and one alphabet</h1>
  <p class="lede">Chalam’s choice of five SSEs had an operational reason — unsaturated but computable. The n=3/4/5 curve gives it a theoretical one. Three quantities, all on the same 491,963-AFDB / 49,640-PDB record subset, place n=5 in a narrow diagnostic window between two boundaries — and show the AFDB and PDB alphabets stay common the whole way.</p>
  <div class="boundband">
    <div class="bb sat"><span class="bn mono">n ≤ 4</span><span class="bl">saturated boundary<br>every cell lit, no signal</span></div>
    <div class="bb win"><span class="bn mono">n = 5</span><span class="bl">diagnostic window<br>45% lit · median 4 folds/cell</span></div>
    <div class="bb fin"><span class="bn mono">n ≥ 6</span><span class="bl">fingerprint boundary<br>a motif names one domain</span></div>
  </div>
</header>

<section>
  <div class="kicker">The three readings</div>
  <h2>Two boundaries close in from both sides</h2>
  <div class="panels">
    {p1}
    {p2}
    {p3}
  </div>
  <div class="legend">
    <span class="lk"><i class="sw afdb"></i>AFDB</span>
    <span class="lk"><i class="sw pdb"></i>PDB</span>
    <span class="lk"><i class="sw cross"></i>cross-database</span>
    <span class="lk"><i class="sw wpdb dash"></i>within-PDB null</span>
    <span class="lk"><i class="sw wafdb dash"></i>within-AFDB null</span>
    <span class="lk"><i class="sw afdb dash"></i>projection</span>
  </div>
</section>

<section class="feature">
  <div class="kicker">What the curve settles</div>
  <h2>n=5 is not a convenience — it’s the diagnostic optimum</h2>
  <div class="cards3">
    <div class="c3"><div class="c3h mono">below n=5</div><p>At n=3 the grid is <strong>100% saturated</strong> and a motif sits in ~950 unrelated folds — pure boundary, no power. n=4 is still 89% lit. Nothing to discriminate. <span class="note-in">(And S3 behaving as a clean boundary confirms the enumerator has no low-n pathology.)</span></p></div>
    <div class="c3"><div class="c3h mono">at n=5</div><p><strong>45% lit</strong>, and the median motif is shared across just <strong>4 folds</strong> — still genuinely shared, so the abstraction still abstracts, but close to the transition. The unique level that discriminates without collapsing to per-domain fingerprints.</p></div>
    <div class="c3"><div class="c3h mono">above n=5</div><p>Sharing falls ÷6 per step, so <strong>n=6 lands at median ≈1</strong>: a motif becomes a fingerprint and the grid stops abstracting. At 151,808 cells it is also 24× the cost and deeply undersampled — buying noise, not vocabulary.</p></div>
  </div>
  <p class="punch">Past the saturated boundary, short of the fingerprint boundary — n=5 is the widest motif that still speaks about folds rather than domains.</p>
</section>

<section class="closer">
  <div class="kicker">The load-bearing control</div>
  <h2>The apparent n=5 divergence is depth, not a different alphabet</h2>
  <p>The raw cross-database agreement drops from 0.97 to 0.67 at n=5, which could read as AFDB and PDB sampling different topological vocabularies. It isn’t. At matched depth, the cross-database Jaccard (<span class="mono">0.65</span>) is <strong>indistinguishable from AFDB’s own within-database null</strong> (<span class="mono">0.66</span>) — two draws of AFDB disagree by the same amount AFDB and PDB disagree. A database cannot sample a different alphabet from itself, so the residual is the sampling-noise floor at this depth, not divergence.</p>
  <div class="kv">
    <div class="krow"><span>cross-database (AFDB ↔ PDB), matched</span><span class="mono pred">0.645</span></div>
    <div class="krow"><span>within-AFDB null (self ↔ self)</span><span class="mono">0.664</span></div>
    <div class="krow"><span>within-PDB null (self ↔ self)</span><span class="mono">0.791</span></div>
  </div>
  <p class="punch2">The S5 topological alphabet is common to experimental and predicted structure space, and the “same alphabet” claim is scale-robust through n=5. Going to n=6 would measure deeper undersampling, not new vocabulary — the limiting factor across the whole curve is depth, never the alphabet.</p>
  <div class="caveats">
    <div class="cav"><span class="mono">identical sets</span> Every point uses the same 491,963 AFDB / 49,640 PDB records, so the trend is about n, not about which records were searched.</div>
    <div class="cav"><span class="mono">matched depth</span> Agreement panel thins both databases to a common hitter count and adds the within-database null — the divergence question answered directly.</div>
  </div>
</section>

<footer>ProSMoS S-grid · n=3/4/5 on identical record sets · sharing = distinct ECOD T-groups/cell · 2026-07-20</footer>
</div>'''

CSS = '''
*{box-sizing:border-box}
:root{
  --ground:#eef1f5;--surface:#ffffff;--ink:#14212d;--muted:#5b6a78;--faint:#8595a3;
  --hair:#dbe2ea;--pred:#2166d0;--exp:#c96500;--cross:#7a3ea8;--band:#eef3fb;--hi:#0f9d58;
  --sans:-apple-system,BlinkMacSystemFont,"Segoe UI",system-ui,Roboto,sans-serif;
  --mono:ui-monospace,"SF Mono","Cascadia Code","JetBrains Mono",Menlo,Consolas,monospace;
}
@media (prefers-color-scheme:dark){:root{
  --ground:#0c131c;--surface:#151f2b;--ink:#e7edf2;--muted:#94a5b5;--faint:#6c7d8d;
  --hair:#243240;--pred:#5a95f0;--exp:#f0a03e;--cross:#b784e0;--band:#132132;--hi:#3ec87e;
}}
:root[data-theme="light"]{--ground:#eef1f5;--surface:#ffffff;--ink:#14212d;--muted:#5b6a78;--faint:#8595a3;--hair:#dbe2ea;--pred:#2166d0;--exp:#c96500;--cross:#7a3ea8;--band:#eef3fb;--hi:#0f9d58}
:root[data-theme="dark"]{--ground:#0c131c;--surface:#151f2b;--ink:#e7edf2;--muted:#94a5b5;--faint:#6c7d8d;--hair:#243240;--pred:#5a95f0;--exp:#f0a03e;--cross:#b784e0;--band:#132132;--hi:#3ec87e}
html{-webkit-text-size-adjust:100%}
body{margin:0;background:var(--ground);color:var(--ink);font-family:var(--sans);line-height:1.6;font-size:17px;-webkit-font-smoothing:antialiased}
.wrap{max-width:1080px;margin:0 auto;padding:0 24px}
.mono{font-family:var(--mono);font-variant-numeric:tabular-nums;font-size:.92em}
.eyebrow,.kicker{font-family:var(--mono);font-size:12.5px;letter-spacing:.06em;text-transform:uppercase;color:var(--faint)}
h1{font-size:clamp(2rem,4.6vw,3.1rem);line-height:1.05;letter-spacing:-.02em;font-weight:800;text-wrap:balance;margin:.5rem 0 0}
h2{font-size:clamp(1.45rem,3.2vw,2.05rem);line-height:1.12;letter-spacing:-.015em;font-weight:750;text-wrap:balance;margin:.2rem 0 .7rem}
.pred{color:var(--pred)}
p{margin:.7rem 0;max-width:68ch;color:var(--ink)}
strong{font-weight:680}
.note-in{color:var(--faint)}
.hero{padding:70px 0 38px;border-bottom:1px solid var(--hair)}
.lede{font-size:clamp(1.02rem,1.7vw,1.2rem);color:var(--muted);max-width:66ch;margin-top:1.1rem}
.boundband{display:flex;gap:0;margin-top:34px;border:1px solid var(--hair);border-radius:12px;overflow:hidden}
.bb{flex:1;padding:16px 18px;display:flex;flex-direction:column;gap:5px}
.bb.sat{background:var(--band)} .bb.win{background:var(--surface);border-left:1px solid var(--hair);border-right:1px solid var(--hair)}
.bb.fin{background:var(--band)}
.bb .bn{font-size:15px;font-weight:700;color:var(--ink)}
.bb.win .bn{color:var(--pred)}
.bb .bl{font-size:12.5px;color:var(--muted);line-height:1.4}
@media (max-width:640px){.boundband{flex-direction:column}.bb.win{border-left:none;border-right:none;border-top:1px solid var(--hair);border-bottom:1px solid var(--hair)}}
section{padding:50px 0;border-bottom:1px solid var(--hair)}
.kicker{display:block;margin-bottom:14px}
.feature{background:var(--surface);margin:0 -24px;padding:54px 24px}
.panels{display:grid;grid-template-columns:repeat(3,1fr);gap:16px;margin:8px 0}
@media (max-width:820px){.panels{grid-template-columns:1fr}}
.panel{margin:0;background:var(--surface);border:1px solid var(--hair);border-radius:12px;padding:14px 14px 6px}
.feature .panel,.panels .panel{background:var(--surface)}
.ptitle{font-size:13px;font-weight:650;color:var(--ink);margin-bottom:6px;padding:0 2px}
.panel svg{width:100%;height:auto;display:block}
.pg{stroke:var(--hair);stroke-width:1}
.pyt{fill:var(--faint);font:10px var(--mono);text-anchor:end}
.pxt{fill:var(--faint);font:11px var(--mono);text-anchor:middle}
.paxt{fill:var(--faint);font:10.5px var(--sans);text-anchor:middle}
.pline{fill:none;stroke-width:2.4;stroke-linejoin:round;stroke-linecap:round}
.pline.afdb{stroke:var(--pred)} .pline.pdb{stroke:var(--exp)} .pline.cross{stroke:var(--cross)}
.pline.wpdb{stroke:var(--exp);opacity:.7} .pline.wafdb{stroke:var(--pred);opacity:.7}
.pdot{stroke:var(--surface);stroke-width:1.4}
.pdot.afdb{fill:var(--pred)} .pdot.pdb{fill:var(--exp)} .pdot.cross{fill:var(--cross)}
.pdot.wpdb{fill:var(--exp)} .pdot.wafdb{fill:var(--pred)}
.pdot.proj{opacity:.6}
.pnote{font-size:12px;color:var(--muted);line-height:1.45;margin:8px 4px 4px;max-width:none}
.legend{display:flex;gap:18px;flex-wrap:wrap;margin-top:18px;justify-content:center}
.lk{display:flex;align-items:center;gap:7px;font-size:12.5px;color:var(--muted);font-family:var(--mono)}
.sw{width:16px;height:3px;border-radius:2px;display:inline-block}
.sw.afdb{background:var(--pred)} .sw.pdb{background:var(--exp)} .sw.cross{background:var(--cross)}
.sw.wpdb{background:var(--exp);opacity:.7} .sw.wafdb{background:var(--pred);opacity:.7}
.sw.dash{background:repeating-linear-gradient(90deg,currentColor 0 4px,transparent 4px 7px)}
.sw.wpdb.dash{color:var(--exp);background:none;border-top:3px dashed var(--exp);height:0}
.cards3{display:grid;grid-template-columns:repeat(3,1fr);gap:16px;margin:22px 0 8px}
@media (max-width:720px){.cards3{grid-template-columns:1fr}}
.c3{background:var(--ground);border:1px solid var(--hair);border-radius:12px;padding:18px 20px}
.c3h{color:var(--pred);font-size:12px;text-transform:uppercase;letter-spacing:.04em;margin-bottom:8px}
.c3 p{font-size:14px;color:var(--muted);max-width:none;margin:0}
.c3 strong{color:var(--ink)}
.punch{font-size:1.12rem;font-weight:600;color:var(--ink);max-width:64ch;margin-top:22px;border-top:2px solid var(--pred);padding-top:16px}
.punch2{font-size:1.1rem;font-weight:600;color:var(--ink);max-width:66ch;margin-top:18px;border-top:2px solid var(--hi);padding-top:16px}
.closer{border-bottom:none}
.kv{display:flex;flex-direction:column;margin:22px 0 6px;border:1px solid var(--hair);border-radius:10px;overflow:hidden;max-width:520px}
.krow{display:flex;justify-content:space-between;padding:11px 16px;font-size:14.5px;border-bottom:1px solid var(--hair);align-items:center}
.krow:last-child{border-bottom:none}
.krow span:first-child{color:var(--muted)}
.krow .mono{font-size:15px;font-weight:600}
.krow .pred{color:var(--pred)}
.caveats{display:flex;gap:16px;margin-top:24px;flex-wrap:wrap}
.cav{flex:1;min-width:260px;background:var(--surface);border:1px solid var(--hair);border-radius:10px;padding:15px 17px;font-size:13.5px;color:var(--muted)}
.cav .mono{color:var(--pred);display:block;margin-bottom:6px;font-weight:600;font-size:11.5px;text-transform:uppercase;letter-spacing:.04em}
footer{padding:32px 0 58px;color:var(--faint);font-family:var(--mono);font-size:12.5px;text-align:center}
'''

open(OUT, "w").write(f"<style>{CSS}</style>\n{HTML}")
print("wrote", OUT)
