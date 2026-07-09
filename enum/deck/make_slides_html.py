#!/usr/bin/env python3
"""Render the Marp markdown deck to a SINGLE standalone HTML: each '---'-separated slide becomes a
card, images are base64-embedded (no external files), Marp '![w:N](...)' sizing honored. No JS/CDN."""
import re, sys, os, base64, markdown
SRC = sys.argv[1] if len(sys.argv) > 1 else "PDB_ECOD_2026_overview.md"
OUT = os.path.splitext(SRC)[0] + ".html"
base = os.path.dirname(os.path.abspath(SRC))
text = open(SRC).read()
# strip YAML frontmatter, grab title
title = "PDB / ECOD overview"
if text.startswith("---"):
    m = re.match(r"^---\n(.*?)\n---\n", text, re.S)
    if m:
        for line in m.group(1).splitlines():
            if line.startswith("title:"): title = line.split(":",1)[1].strip().strip('"')
        text = text[m.end():]
slides = [s.strip() for s in re.split(r"(?m)^---\s*$", text) if s.strip()]

def embed(html):
    def repl(mt):
        attrs = mt.group(0); src = mt.group("src"); alt = mt.group("alt") or ""
        style = ""
        w = re.search(r"w:(\d+)", alt)
        if w: style = f' style="width:{w.group(1)}px;max-width:100%"'
        p = src if os.path.isabs(src) else os.path.join(base, src)
        if os.path.exists(p):
            b64 = base64.b64encode(open(p,"rb").read()).decode()
            return f'<img src="data:image/png;base64,{b64}"{style}>'
        return f'<img src="{src}"{style}><!-- MISSING -->'
    return re.sub(r'<img[^>]*?alt="(?P<alt>[^"]*)"[^>]*?src="(?P<src>[^"]+)"[^>]*?>|'
                  r'<img[^>]*?src="(?P<src2>[^"]+)"[^>]*?alt="(?P<alt2>[^"]*)"[^>]*?>',
                  lambda m: repl(type("M",(),{"group":lambda s,k:(m.group("src") or m.group("src2")) if k=="src" else (m.group("alt") or m.group("alt2"))})()),
                  html)

md = markdown.Markdown(extensions=["tables","fenced_code","sane_lists"])
cards = []
for i, s in enumerate(slides):
    html = embed(md.reset().convert(s))
    cls = "slide title" if i == 0 else "slide"
    cards.append(f'<section class="{cls}">{html}<div class="pg">{i+1}/{len(slides)}</div></section>')

CSS = """
:root{--fg:#1a2330;--mut:#5b6b7c;--ac:#2c6e9c}
*{box-sizing:border-box} body{margin:0;background:#e9edf1;font-family:-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;color:var(--fg)}
.slide{position:relative;width:1100px;min-height:620px;margin:22px auto;background:#fff;border:1px solid #d4dae0;border-radius:10px;
 box-shadow:0 2px 10px rgba(0,0,0,.08);padding:46px 56px}
.slide h1{font-size:34px;margin:.1em 0 .3em;color:var(--fg)} .slide.title h1{font-size:40px;margin-top:120px}
.slide h2{font-size:24px;color:var(--ac);margin:.2em 0 .5em;font-weight:600}
.slide p,.slide li{font-size:19px;line-height:1.5;color:#222} .slide li{margin:.25em 0}
.slide em{color:var(--mut)} .slide strong{color:#0b3d61}
.slide blockquote{border-left:4px solid var(--ac);margin:.6em 0;padding:.3em 1em;background:#f3f8fc;color:#234;font-size:18px}
.slide img{display:block;margin:14px auto;border-radius:4px}
.slide table{border-collapse:collapse;margin:12px auto;font-size:18px} .slide th,.slide td{border:1px solid #cdd5dd;padding:6px 12px;text-align:center}
.slide th{background:#f0f4f8}
.pg{position:absolute;bottom:14px;right:22px;color:#9aa7b3;font-size:13px}
@media print{body{background:#fff}.slide{margin:0;border:none;box-shadow:none;border-radius:0;page-break-after:always;width:100%;min-height:96vh}}
"""
open(OUT,"w").write(f"<!doctype html><html><head><meta charset='utf-8'><title>{title}</title><style>{CSS}</style></head><body>{''.join(cards)}</body></html>")
print(f"-> {OUT}  ({len(slides)} slides, {os.path.getsize(OUT)//1024} KB)")
