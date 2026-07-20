import json, math, base64, io
from pathlib import Path
from PIL import Image

D = Path("/home/rschaeff/work/prosmos_2026/design_scaffold")
OUT = "/tmp/claude-1219/-home-rschaeff-dev-prosmos-cl/02f20625-920b-46f1-a2bf-bc06d84727af/scratchpad/design_scaffold_deck.html"


def uri(fname, maxw=560):
    im = Image.open(D / "renders" / fname)      # keep alpha (transparent bg)
    if im.width > maxw:
        im = im.resize((maxw, round(im.height * maxw / im.width)), Image.LANCZOS)
    buf = io.BytesIO()
    im.save(buf, "PNG", optimize=True)
    return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode()


def sse_str(ty):
    # 5-bit typing, MSB = SSE1; bit 1 -> E, 0 -> H  (ty31=EEEEE, ty0=HHHHH)
    return "".join("E" if (ty >> (4 - i)) & 1 else "H" for i in range(5))


EX = json.load(open(D / "exemplars6.json"))
ex_order = ["43_0", "74_13", "40_0", "82_5", "27_0", "29_25"]
exemplar_cards = ""
for key in ex_order:
    e = EX[key]
    img = uri(f"ex_{key}.png")
    topo = sse_str(e["ty"])
    exemplar_cards += (
        f'<figure class="excard">'
        f'<div class="eximg"><img src="{img}" alt="AlphaFold scaffold for topology sk{e["sk"]:03d} ty{e["ty"]:02d}, coloured by pLDDT"/>'
        f'<span class="plbadge mono">pLDDT {e["plddt"]:.0f}</span></div>'
        f'<figcaption><span class="ectopo mono">{topo}</span>'
        f'<span class="eccell mono">sk{e["sk"]:03d}·ty{e["ty"]:02d}</span>'
        f'<span class="eccount"><b class="exp">{e["pdb"]}</b> PDB → '
        f'<b class="pred">{e["afdb70"]}</b> AFDB families</span></figcaption>'
        f'</figure>')

H = json.load(open(D / "mult_hist.json"))
V = json.load(open(D / "perdomain_validated.json"))
cells = sorted(V["cells"], key=lambda c: -c["afdb70_dom"])

# ---- multiplier distribution histogram (log-x) ----
edges, hist = H["edges"], H["hist"]
HX0, HX1, HY0, HY1 = 70, 712, 40, 300
lo, hi = math.log10(edges[0]), math.log10(edges[-1])
def hx(v): return HX0 + (math.log10(v) - lo) / (hi - lo) * (HX1 - HX0)
hmax = max(hist)
def hy(n): return HY1 - n / hmax * (HY1 - HY0)
bars = ""
for i, n in enumerate(hist):
    x0, x1 = hx(edges[i]), hx(edges[i + 1])
    cls = "hbar lose" if edges[i + 1] <= 1 else "hbar win"
    bars += (f'<rect x="{x0:.1f}" y="{hy(n):.1f}" width="{x1-x0-1.5:.1f}" height="{HY1-hy(n):.1f}" '
             f'class="{cls}"><title>{n} cells, {edges[i]}–{edges[i+1]}×</title></rect>')
xt = ""
for t in [0.5, 1, 2, 4, 8, 20, 120]:
    xt += (f'<text x="{hx(t):.1f}" y="{HY1+18}" class="hxtick">{("%g"%t)}×</text>')
yt = ""
for n in [0, 150, 300, 450]:
    yt += (f'<line x1="{HX0}" y1="{hy(n):.1f}" x2="{HX1}" y2="{hy(n):.1f}" class="hgrid"/>'
           f'<text x="{HX0-8}" y="{hy(n)+4:.1f}" class="hytick">{n}</text>')
medx = hx(H["median"])
onex = hx(1.0)
hist_svg = f'''<svg viewBox="0 0 740 340" class="chart" role="img" aria-label="Distribution of the per-topology sequence-diversity multiplier">
  {yt}
  <line x1="{onex:.1f}" y1="{HY0}" x2="{onex:.1f}" y2="{HY1}" class="hone"/>
  <text x="{onex-5:.1f}" y="{HY0+12}" class="hlabel" text-anchor="end">parity</text>
  {bars}
  <line x1="{medx:.1f}" y1="{HY0-6}" x2="{medx:.1f}" y2="{HY1}" class="hmed"/>
  <text x="{medx+6:.1f}" y="{HY0+4}" class="hmedlabel">median 4.6×</text>
  {xt}
  <text x="{(HX0+HX1)/2:.0f}" y="335" class="haxt">AFDB sequence families ÷ PDB sequence families, per shared topology (log scale)</text>
</svg>'''

# ---- the 31-cell table rows ----
def rowcls(c):
    return "grade80" if c["afdb80_dom"] >= 50 else "grade70"
trows = ""
for c in cells:
    mult80 = c["afdb80_dom"] / max(c["pdb"], 1)
    tag = '<span class="g80">≥80 ✓</span>' if c["afdb80_dom"] >= 50 else '<span class="g70">≥70</span>'
    trows += (f'<div class="crow {rowcls(c)}">'
              f'<span class="mono cid">sk{c["sk"]:03d}·ty{c["ty"]:02d}</span>'
              f'<span class="mono pdbn">{c["pdb"]}</span>'
              f'<span class="mono afn">{c["afdb70_dom"]}</span>'
              f'<span class="mono af8">{c["afdb80_dom"]}</span>'
              f'<span class="mono mx">{mult80:.0f}×</span>'
              f'<span class="gtag">{tag}</span></div>')

HTML = f'''<div class="wrap">

<header class="hero">
  <div class="eyebrow">ProSMoS S5 · AlphaFold DB vs experimental PDB · a scaffold-diversity readout for design</div>
  <h1>The 200M's gift to design isn't new folds —<br>it's <span class="pred">sequence-diverse scaffolds</span> over proven space</h1>
  <p class="lede">The same result that says predicted structures hold no new 5-SSE topology says something a design audience actually wants: over the shared, experimentally-validated vocabulary of local topologies, AlphaFold DB multiplies the pool of sequence-diverse, high-confidence backbones many-fold — and it concentrates that gift on exactly the topologies the PDB barely templates.</p>
  <div class="statrow">
    <div class="stat hi"><div class="n">4.6×</div><div class="l">median sequence-family<br>gain per topology</div></div>
    <div class="stat"><div class="n">~10×</div><div class="l">total scaffold pool<br>over shared space</div></div>
    <div class="stat"><div class="n">31</div><div class="l">pLDDT-validated<br>under-templated cells</div></div>
    <div class="stat"><div class="n">95<span class="of">%</span></div><div class="l">of topologies gain<br>diversity from AFDB</div></div>
  </div>
</header>

<section>
  <div class="kicker">The reframe</div>
  <h2>Every result that disappoints discovery <em>delivers</em> for design</h2>
  <p>We spent this project proving the 200M adds no novel local topology. Term for term, that negative is a design positive — because designers don't want exotic topology, they want many sequence-diverse, foldable realizations of <em>known-good</em> backbones.</p>
  <div class="flip2">
    <div class="fcol fneg"><div class="fh">What we proved (discovery)</div>
      <div class="fr">AFDB S5 space ≈ PDB S5 space topologically</div>
      <div class="fr">“AFDB-only cells” are mostly sampling depth</div>
      <div class="fr">A cell lit by many AFDB, few PDB domains</div>
      <div class="fr">Our AFDB set is non-singleton cluster reps</div>
      <div class="fr">52% of lit AFDB domains are “TED-novel”</div>
    </div>
    <div class="fcol fpos"><div class="fh">Same fact (design)</div>
      <div class="fr">The 200M covers the <strong>real, useful</strong> vocabulary — nothing exotic to validate</div>
      <div class="fr"><strong>10× the scaffold sampling</strong> over that vocabulary</div>
      <div class="fr">That topology is <strong>under-templated</strong> in the PDB, richly populated in AFDB</div>
      <div class="fr">A <strong>pre-deduplicated, sequence-diverse</strong> scaffold set by construction</div>
      <div class="fr">Half your scaffolds carry <strong>sequence signal the PDB never had</strong></div>
    </div>
  </div>
  <p class="note">The depth argument that <em>kills</em> the novelty claim is the argument that <em>makes</em> the scaffold claim. Nothing gets walked back — the axis is relabelled.</p>
</section>

<section class="feature">
  <div class="kicker">The measurement · matched 50%/90% clustering, both databases</div>
  <h2>4.6× more sequence families per topology — like-for-like</h2>
  <p>Naively, ECOD family counts put the multiplier at 8×. But ECOD F-groups are a coarser unit than AlphaFold’s clustering, so we clustered <em>both</em> databases’ hitting domains identically — MMseqs2 at 50% identity, 90% coverage — and counted distinct clusters per topology. The honest, matched number is <strong>4.6× median</strong>. The correction fell entirely on the PDB side; AFDB’s structurally-distinct reps barely merged (100.0% retained).</p>
  <figure class="chartfig">
    {hist_svg}
    <figcaption>3,086 shared satisfiable topologies · per-cell multiplier = AFDB 50%-clusters ÷ PDB 50%-clusters · median 13 PDB families vs 58 AFDB</figcaption>
  </figure>
  <div class="readout">
    <div class="ro"><span class="mono big pred">95.4%</span><span>of shared topologies have more sequence families in AFDB than the PDB</span></div>
    <div class="rovs">·</div>
    <div class="ro"><span class="mono big pred">867k → 8.56M</span><span>total non-redundant scaffold families over shared space (~10×)</span></div>
  </div>
  <p class="punch">Over the same, experimentally-validated topological vocabulary, predicted structure space offers roughly five times the sequence-diverse scaffolds per topology — and ten times the pool overall.</p>
</section>

<section>
  <div class="kicker">The design target list · confidence-gated</div>
  <h2>31 topologies the PDB barely templates — and AFDB richly does</h2>
  <p>Restricting to topologies where the PDB offers <strong>≤5</strong> sequence families and AFDB offers <strong>≥50</strong> gives the design-relevant targets. Crucially these are <em>shared</em> cells — every one is crystallographically confirmed, so there is no “is it real” risk. A designer currently has almost nothing to template here from the PDB; the 200M supplies dozens of diverse backbones.</p>
  <div class="funnel">
    <div class="fstep"><span class="fn mono">53</span><span class="fl">under-templated cells<br>(raw)</span></div>
    <div class="farrow">→</div>
    <div class="fstep hi"><span class="fn mono">31</span><span class="fl">hold at per-domain<br>pLDDT ≥ 70</span></div>
    <div class="farrow">→</div>
    <div class="fstep hi2"><span class="fn mono">27</span><span class="fl">hold at the strict<br>pLDDT ≥ 80</span></div>
  </div>
  <div class="ctable">
    <div class="crow chead"><span>topology</span><span>PDB fam</span><span>AFDB ≥70</span><span>AFDB ≥80</span><span>×(≥80)</span><span>grade</span></div>
    {trows}
  </div>
  <p class="note">pLDDT is the exact per-domain mean, read from the B-factor of each chopped AlphaFold domain — 2,716 scaffolds, 99% coverage. Green = survives the strict ≥80 gate.</p>
</section>

<section>
  <div class="kicker">What a designer actually receives</div>
  <h2>The highest-confidence scaffold from each of the top six cells</h2>
  <p>One representative AlphaFold backbone per top under-templated topology, coloured by per-residue pLDDT on the standard scheme — <span class="afk high">dark blue ≥90</span>, <span class="afk conf">light blue 70–90</span>. Each label reads the topology’s H/E string, then how the experimental record compares: a handful of PDB families against dozens of confident, sequence-diverse AlphaFold ones behind each of these.</p>
  <div class="exgrid">
    {exemplar_cards}
  </div>
  <p class="note">Backbones are AlphaFold v4 models chopped to the DPAM domain, oriented and rendered fresh. These are the single most-confident example per cell (all pLDDT 95–97); the design payload is the full set of ≥50 diverse families standing behind each.</p>
</section>

<section>
  <div class="kicker">Why the list is trustworthy · the confidence gate was validated, not assumed</div>
  <h2>The scaffolds are well-folded — and we checked it the hard way</h2>
  <p>A design audience will not template off low-confidence models, so the under-templated cells only count if their AFDB backbones are genuinely well-predicted. A fast per-protein pLDDT proxy said yes; we didn’t trust it, because it is exact for only 37% of these scaffolds (the single-domain ones). So we extracted the actual chopped domain structures and read per-residue pLDDT directly.</p>
  <div class="vgrid">
    <div class="vcard"><div class="vk mono">exact per-domain pLDDT</div>
      <div class="vstat"><span class="mono big pred">85.2</span><span>median — comfortably in AlphaFold’s “confident” band, above the 80.8 the proxy estimated</span></div>
      <p class="vfoot">The proxy was <strong>conservative</strong>: whole-protein means are dragged down by disordered regions elsewhere, but the domains that light these cells are well-folded on their own.</p>
    </div>
    <div class="vcard"><div class="vk mono">verdict under exact pLDDT</div>
      <div class="vstat"><span class="mono big hi">31/31</span><span>survive at ≥70; <strong>27/31</strong> survive at ≥80 — 15 more than the proxy allowed</span></div>
      <p class="vfoot">Every cell’s confident-scaffold count went <em>up</em> versus the proxy. The under-templated finding strengthens under scrutiny rather than eroding.</p>
    </div>
  </div>
  <p class="note">Cross-check on the detector itself: on 1,529 same-molecule pairs, PALSSE reads an <strong>identical</strong> SSE inventory from a crystal structure and its AlphaFold model (100%). A mild +14% idealization tilt exists on marginal cells, but the multiplier lives on shared cells where it barely registers — and idealized, constraint-satisfying backbones are precisely what design wants to template.</p>
</section>

<section class="closer">
  <div class="kicker">Takeaway for RosettaCon</div>
  <h2>The 200M is a sequence-diverse scaffold engine over validated topology space</h2>
  <p>Not a source of new folds — a source of <strong>many diverse, high-confidence realizations</strong> of the folds that already work, indexed by local topology, with a concrete list of 31 topologies where it most outstrips the experimental record. That is the contribution a template- or fragment-driven designer can use tomorrow, and it is defensible end to end because every claim is anchored to PDB-confirmed structure.</p>
  <div class="caveats">
    <div class="cav"><span class="mono">conservative</span> AFDB excludes 10.6M singleton clusters we never searched — real scaffolds, so the multiplier is a floor, not a ceiling.</div>
    <div class="cav"><span class="mono">matched</span> Both databases clustered identically at 50%/90%; the 4.6× is like-for-like, not a unit artifact.</div>
    <div class="cav"><span class="mono">validated</span> Per-domain pLDDT from chopped-structure B-factors, not a proxy — the 31-cell list is confidence-gated.</div>
  </div>
</section>

<footer>ProSMoS S5 · design-scaffold analysis · AFDB 4.9M cluster reps / PDB 496k domains · matched 50/90 MMseqs2 · 2026-07-20</footer>
</div>'''

CSS = '''
*{box-sizing:border-box}
:root{
  --ground:#eef1f5; --surface:#ffffff; --ink:#14212d; --muted:#5b6a78; --faint:#8595a3;
  --hair:#dbe2ea; --pred:#2166d0; --exp:#c96500; --pred-soft:#e7effb; --exp-soft:#fbeede;
  --band:#eef3fb; --hi:#0f9d58; --hi-soft:#e4f5ec;
  --sans:-apple-system,BlinkMacSystemFont,"Segoe UI",system-ui,Roboto,sans-serif;
  --mono:ui-monospace,"SF Mono","Cascadia Code","JetBrains Mono",Menlo,Consolas,monospace;
}
@media (prefers-color-scheme:dark){:root{
  --ground:#0c131c; --surface:#151f2b; --ink:#e7edf2; --muted:#94a5b5; --faint:#6c7d8d;
  --hair:#243240; --pred:#5a95f0; --exp:#f0a03e; --pred-soft:#16283f; --exp-soft:#33260f;
  --band:#132132; --hi:#3ec87e; --hi-soft:#122a1e;
}}
:root[data-theme="light"]{--ground:#eef1f5;--surface:#ffffff;--ink:#14212d;--muted:#5b6a78;--faint:#8595a3;--hair:#dbe2ea;--pred:#2166d0;--exp:#c96500;--pred-soft:#e7effb;--exp-soft:#fbeede;--band:#eef3fb;--hi:#0f9d58;--hi-soft:#e4f5ec}
:root[data-theme="dark"]{--ground:#0c131c;--surface:#151f2b;--ink:#e7edf2;--muted:#94a5b5;--faint:#6c7d8d;--hair:#243240;--pred:#5a95f0;--exp:#f0a03e;--pred-soft:#16283f;--exp-soft:#33260f;--band:#132132;--hi:#3ec87e;--hi-soft:#122a1e}
html{-webkit-text-size-adjust:100%}
body{margin:0;background:var(--ground);color:var(--ink);font-family:var(--sans);line-height:1.6;font-size:17px;-webkit-font-smoothing:antialiased}
.wrap{max-width:1080px;margin:0 auto;padding:0 24px}
.mono{font-family:var(--mono);font-variant-numeric:tabular-nums;font-size:.92em}
.eyebrow,.kicker{font-family:var(--mono);font-size:12.5px;letter-spacing:.06em;text-transform:uppercase;color:var(--faint)}
h1{font-size:clamp(2.1rem,5vw,3.4rem);line-height:1.04;letter-spacing:-.02em;font-weight:800;text-wrap:balance;margin:.5rem 0 0}
h2{font-size:clamp(1.5rem,3.4vw,2.15rem);line-height:1.12;letter-spacing:-.015em;font-weight:750;text-wrap:balance;margin:.2rem 0 .7rem}
.pred{color:var(--pred)} .exp{color:var(--exp)} .hi{color:var(--hi)}
em{font-style:normal;color:var(--pred);font-weight:650}
h2 em{border-bottom:2.5px solid var(--pred);padding-bottom:1px}
p{margin:.7rem 0;max-width:66ch;color:var(--ink)}
strong{font-weight:680}
.note{color:var(--muted);font-size:15.5px;border-left:2px solid var(--hair);padding-left:14px}
.hero{padding:72px 0 40px;border-bottom:1px solid var(--hair)}
.lede{font-size:clamp(1.02rem,1.8vw,1.22rem);color:var(--muted);max-width:64ch;margin-top:1.1rem}
.statrow{display:flex;align-items:stretch;gap:16px;margin-top:36px;flex-wrap:wrap}
.stat{background:var(--surface);border:1px solid var(--hair);border-radius:12px;padding:18px 22px;min-width:150px;flex:1}
.stat.hi{border-color:var(--pred);box-shadow:0 0 0 1px var(--pred)}
.stat .n{font-family:var(--mono);font-size:2.4rem;font-weight:700;line-height:1;color:var(--ink);font-variant-numeric:tabular-nums}
.stat.hi .n{color:var(--pred)}
.stat .n .of{font-size:1.2rem;color:var(--faint)}
.stat .l{color:var(--muted);font-size:13px;margin-top:8px;line-height:1.35}
section{padding:52px 0;border-bottom:1px solid var(--hair)}
.kicker{display:block;margin-bottom:14px}
.feature{background:var(--surface);border-radius:0;margin:0 -24px;padding:56px 24px}
.flip2{display:grid;grid-template-columns:1fr 1fr;gap:0;margin:24px 0 8px;border:1px solid var(--hair);border-radius:12px;overflow:hidden}
.fcol{padding:0}
.fh{font-family:var(--mono);font-size:12px;text-transform:uppercase;letter-spacing:.04em;padding:12px 18px;color:var(--faint);background:var(--ground)}
.fneg .fh{color:var(--exp)} .fpos .fh{color:var(--pred)}
.fr{padding:12px 18px;font-size:14px;border-top:1px solid var(--hair);min-height:64px;display:flex;align-items:center}
.fneg{border-right:1px solid var(--hair)} .fneg .fr{color:var(--muted)}
.fpos{background:var(--pred-soft)}
@media (max-width:640px){.flip2{grid-template-columns:1fr}.fneg{border-right:none}.fr{min-height:0}}
.chartfig{margin:22px 0 8px;background:var(--surface);border:1px solid var(--hair);border-radius:14px;padding:22px 18px 12px;overflow-x:auto}
.feature .chartfig{background:var(--ground)}
.chart{width:100%;min-width:600px;height:auto;display:block}
.hgrid{stroke:var(--hair);stroke-width:1}
.hytick,.hxtick{fill:var(--faint);font:11px var(--mono)}
.hytick{text-anchor:end}.hxtick{text-anchor:middle}
.haxt{fill:var(--muted);font:11.5px var(--sans);text-anchor:middle}
.hbar{}
.hbar.win{fill:var(--pred);opacity:.86}
.hbar.lose{fill:var(--exp);opacity:.7}
.hone{stroke:var(--faint);stroke-width:1;stroke-dasharray:3 3}
.hlabel{fill:var(--faint);font:11px var(--mono)}
.hmed{stroke:var(--ink);stroke-width:1.6}
.hmedlabel{fill:var(--ink);font:600 12px var(--mono)}
figcaption{color:var(--faint);font-size:13px;margin-top:10px;font-family:var(--mono);line-height:1.5}
.readout{display:flex;align-items:center;gap:22px;margin:26px 0 6px;flex-wrap:wrap}
.ro{display:flex;align-items:baseline;gap:12px}
.ro span:last-child{color:var(--muted);font-size:14.5px;max-width:26ch}
.big{font-size:2.4rem;font-weight:700;line-height:1}
.rovs{color:var(--faint);font-family:var(--mono)}
.punch{font-size:1.13rem;font-weight:600;color:var(--ink);max-width:64ch;margin-top:22px;border-top:2px solid var(--pred);padding-top:18px}
.funnel{display:flex;align-items:center;gap:16px;margin:26px 0;flex-wrap:wrap}
.fstep{background:var(--surface);border:1px solid var(--hair);border-radius:12px;padding:16px 24px;text-align:center;min-width:150px}
.fstep.hi{border-color:var(--pred);box-shadow:0 0 0 1px var(--pred)}
.fstep.hi2{border-color:var(--hi);box-shadow:0 0 0 1px var(--hi)}
.fstep .fn{font-size:2.2rem;font-weight:700;display:block;line-height:1}
.fstep.hi .fn{color:var(--pred)} .fstep.hi2 .fn{color:var(--hi)}
.fstep .fl{font-size:12.5px;color:var(--muted);margin-top:6px;display:block}
.farrow{color:var(--faint);font-size:1.5rem}
.ctable{margin:22px 0;border:1px solid var(--hair);border-radius:10px;overflow:hidden;background:var(--surface)}
.crow{display:grid;grid-template-columns:1.4fr .8fr .9fr .9fr .7fr 1fr;gap:8px;padding:8px 16px;align-items:center;font-size:14px;border-bottom:1px solid var(--hair)}
.crow:last-child{border-bottom:none}
.chead{background:var(--ground);color:var(--faint);font-family:var(--mono);font-size:11px;text-transform:uppercase;letter-spacing:.03em}
.crow span:not(.cid):not(:last-child){text-align:right;font-variant-numeric:tabular-nums}
.chead span:not(:first-child){text-align:right}
.cid{color:var(--ink);font-size:13px}
.pdbn{color:var(--exp);font-weight:600}
.afn{color:var(--pred);font-weight:600}
.af8{color:var(--muted)}
.mx{color:var(--ink)}
.grade80{background:linear-gradient(90deg,var(--hi-soft),transparent 60%)}
.gtag{text-align:right}
.g80{font-family:var(--mono);font-size:10.5px;color:var(--hi);border:1px solid var(--hi);border-radius:4px;padding:1px 6px}
.g70{font-family:var(--mono);font-size:10.5px;color:var(--faint);border:1px solid var(--hair);border-radius:4px;padding:1px 6px}
@media (max-width:640px){.crow{grid-template-columns:1.3fr .7fr .8fr .8fr;font-size:12.5px}.af8,.chead .h8,.gtag{display:none}}
.vgrid{display:flex;gap:18px;margin:24px 0 8px;flex-wrap:wrap}
.vcard{flex:1 1 320px;background:var(--surface);border:1px solid var(--hair);border-radius:12px;padding:20px 22px}
.vk{color:var(--pred);font-size:12px;letter-spacing:.03em;text-transform:uppercase;margin-bottom:12px}
.vstat{display:flex;align-items:baseline;gap:14px;margin:14px 0;padding:14px 0;border-top:1px solid var(--hair);border-bottom:1px solid var(--hair)}
.vstat span:last-child{font-size:14px;color:var(--muted)}
.vstat .hi{color:var(--hi)}
.vfoot{font-size:13.5px;color:var(--faint)!important;line-height:1.5}
.afk{font-family:var(--mono);font-size:.85em;padding:1px 6px;border-radius:4px;white-space:nowrap}
.afk.high{background:#0653D0;color:#fff}
.afk.conf{background:#65cbf3;color:#0a2a3a}
.exgrid{display:grid;grid-template-columns:repeat(3,1fr);gap:16px;margin:24px 0 8px}
@media (max-width:820px){.exgrid{grid-template-columns:repeat(2,1fr)}}
@media (max-width:520px){.exgrid{grid-template-columns:1fr}}
.excard{margin:0;background:var(--surface);border:1px solid var(--hair);border-radius:12px;overflow:hidden}
.eximg{position:relative;background:linear-gradient(160deg,#f3f7fc,#e7eef7);aspect-ratio:1/1;display:flex;align-items:center;justify-content:center}
.eximg img{width:100%;height:100%;object-fit:contain;display:block}
.plbadge{position:absolute;top:9px;right:9px;background:rgba(6,83,208,.92);color:#fff;font-size:11px;padding:2px 7px;border-radius:5px;letter-spacing:.02em}
.excard figcaption{padding:12px 14px 13px;display:flex;flex-direction:column;gap:3px;margin:0}
.ectopo{font-size:15px;font-weight:700;letter-spacing:.14em;color:var(--ink)}
.eccell{font-size:11px;color:var(--faint)}
.eccount{font-size:12.5px;color:var(--muted);margin-top:3px}
.eccount b{font-weight:700}
.closer{border-bottom:none}
.closer h2{max-width:22ch}
.caveats{display:flex;gap:16px;margin-top:26px;flex-wrap:wrap}
.cav{flex:1;min-width:250px;background:var(--surface);border:1px solid var(--hair);border-radius:10px;padding:15px 17px;font-size:14px;color:var(--muted)}
.cav .mono{color:var(--pred);display:block;margin-bottom:6px;font-weight:600;font-size:12px;text-transform:uppercase;letter-spacing:.04em}
footer{padding:34px 0 60px;color:var(--faint);font-family:var(--mono);font-size:12.5px;text-align:center}
@media (max-width:640px){body{font-size:16px}.statrow{gap:12px}.farrow{display:none}}
'''

doc = f"<style>{CSS}</style>\n{HTML}"
open(OUT, "w").write(doc)
print("wrote", OUT, round(len(doc) / 1024), "KB")
