#!/usr/bin/env python3
"""Render one or more Marp markdown decks into a SINGLE multi-page PDF via matplotlib (no browser/LaTeX
needed). Each '---' slide -> one landscape page: title, subtitle, images (single centered, or a row from a
markdown image-table with header captions), and bullet text. Usage: make_slides_pdf.py out.pdf deck1.md [deck2.md ...]"""
import sys, os, re, textwrap
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

OUT=sys.argv[1]; DECKS=sys.argv[2:]
IMG=re.compile(r'!\[(?P<alt>[^\]]*)\]\((?P<src>[^)]+)\)')

def parse_slide(text):
    title=""; subs=[]; bullets=[]; paras=[]; imgs=[]; captions=[]
    for raw in text.splitlines():
        line=raw.rstrip()
        if not line.strip(): continue
        if line.startswith("# ") and not title: title=line[2:].strip(); continue
        if line.startswith(("## ","### ")): subs.append(re.sub(r'^#+ ','',line).strip()); continue
        # image table row
        if line.startswith("|") and IMG.search(line):
            for m in IMG.finditer(line):
                src=m.group("src"); w=re.search(r'w:(\d+)',m.group("alt")); imgs.append((src,int(w.group(1)) if w else 0))
            continue
        if line.startswith("|") and not IMG.search(line) and not re.match(r'^\|[\s:|-]+\|?$',line):
            captions=[c.strip() for c in line.strip("|").split("|") if c.strip()]; continue
        if re.match(r'^\|[\s:|-]+\|?$',line): continue
        m=IMG.search(line)
        if m:
            w=re.search(r'w:(\d+)',m.group("alt")); imgs.append((m.group("src"),int(w.group(1)) if w else 0)); continue
        if line.startswith(("- ","* ")): bullets.append(line[2:].strip()); continue
        if line.startswith(">"): paras.append(line.lstrip("> ").strip()); continue
        paras.append(line.strip())
    return title,subs,bullets,paras,imgs,captions

def clean(s):  # strip md emphasis/code/links for plain rendering
    s=s.replace('**','').replace('*','')           # drop all emphasis markers (handles nesting)
    s=re.sub(r'`([^`]+)`',r'\1',s); s=re.sub(r'\[([^\]]+)\]\([^)]+\)',r'\1',s)
    return s

def place_img(fig,path,box):  # box=(x,y,w,h) fig-fraction; fit preserving aspect
    if not os.path.exists(path): return
    im=plt.imread(path); ih,iw=im.shape[0],im.shape[1]; ar=iw/ih
    bx,by,bw,bh=box; bar=(bw*fig.get_figwidth())/(bh*fig.get_figheight())
    if ar>bar: w=bw; h=bw*fig.get_figwidth()/(ar*fig.get_figheight())
    else: h=bh; w=bh*fig.get_figheight()*ar/fig.get_figwidth()
    ax=fig.add_axes([bx+(bw-w)/2, by+(bh-h)/2, w, h]); ax.imshow(im); ax.axis("off")

def render(slides, pdf):
    for text in slides:
        title,subs,bullets,paras,imgs,caps=parse_slide(text)
        fig=plt.figure(figsize=(13.33,7.5)); fig.patch.set_facecolor("white")
        y=0.93
        if title: fig.text(0.06,y,clean(title),fontsize=26,fontweight="bold",color="#13233a",va="top"); y-=0.075
        for s in subs: fig.text(0.06,y,clean(s),fontsize=15,color="#2c6e9c",va="top"); y-=0.05
        top=y-0.01
        # place images; set ty = where body text starts (always render bullets/paras below images)
        if imgs and len(imgs)==1:
            place_img(fig,imgs[0][0],(0.10,0.40,0.80,top-0.42)); ty=0.36; bf=12.5; wrap=110
        elif imgs:
            gw=0.92/len(imgs)
            for i,(src,_) in enumerate(imgs):
                place_img(fig,src,(0.04+i*gw,top-0.40,gw-0.01,0.38))
                if i<len(caps): fig.text(0.04+i*gw+gw/2,top-0.42,caps[i],ha="center",fontsize=11,color="#444")
            ty=top-0.48; bf=11; wrap=112
        else:
            ty=top; bf=15; wrap=105
        for b in bullets:
            wrapped=textwrap.fill(clean(b),wrap)
            fig.text(0.06,ty,"• "+wrapped,fontsize=bf,va="top",color="#222")
            ty-=0.030*bf/12*(wrapped.count(chr(10))+1)+0.012
        for p in paras:
            if p==title or p in subs: continue
            wrapped=textwrap.fill(clean(p),wrap+10)
            fig.text(0.06,ty,wrapped,fontsize=max(bf-3,10),va="top",color="#555",style="italic")
            ty-=0.030*(max(bf-3,10))/12*(wrapped.count(chr(10))+1)+0.010
        pdf.savefig(fig); plt.close(fig)

with PdfPages(OUT) as pdf:
    for deck in DECKS:
        base=os.path.dirname(os.path.abspath(deck)); os.chdir(base)
        t=open(deck).read()
        t=re.sub(r'<!--.*?-->','',t,flags=re.S)   # drop speaker-note HTML comments (PDF renderer has no comment handling)
        if t.startswith("---"):
            m=re.match(r"^---\n.*?\n---\n",t,re.S); t=t[m.end():] if m else t
        slides=[s.strip() for s in re.split(r"(?m)^---\s*$",t) if s.strip()]
        render(slides,pdf)
print("wrote",OUT)
