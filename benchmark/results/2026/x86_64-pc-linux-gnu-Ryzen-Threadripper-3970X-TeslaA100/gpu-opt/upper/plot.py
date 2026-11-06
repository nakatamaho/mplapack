#!/usr/bin/env python3
"""Dependency-free plots of the paired item4b whole-entry timings."""
import csv,html,math
from pathlib import Path
p=Path(__file__).resolve().parent
rows=list(csv.DictReader((p/'matrix.csv').open()))
order={s:i for i,s in enumerate(['dd','td','qd','mpfr512','mpfr256','mpfr1024'])}
def key(r):return order[r['type']],int(r['n'])
def draw(name,title,labels,values,legends,target):
 width=1080;left=220;right=85;top=80;step=38;height=top+step*len(labels)+55
 maximum=math.ceil(max([target]+[v for a in values for v in a])*1.12)
 scale=(width-left-right)/maximum
 s=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',f'<title>{html.escape(title)}</title>','<rect width="100%" height="100%" fill="white"/>','<g font-family="sans-serif" font-size="12" fill="#172b4d">',f'<text x="20" y="27" font-size="18">{html.escape(title)}</text>','<text x="20" y="48">Explicit GPU entry; one warmup, best of three. Ratios are not hardware occupancy.</text>']
 colors=['#2563eb'] if len(legends)==1 else ['#94a3b8','#2563eb']
 for i,legend in enumerate(legends):
  x=left+200*i;s += [f'<rect x="{x}" y="58" width="12" height="10" fill="{colors[i]}"/>',f'<text x="{x+18}" y="68">{html.escape(legend)}</text>']
 for tick in range(maximum+1):
  if maximum>12 and tick%2:continue
  x=left+tick*scale;s += [f'<line x1="{x}" x2="{x}" y1="{top}" y2="{height-40}" stroke="#e2e8f0"/>',f'<text x="{x}" y="{height-20}" text-anchor="middle">{tick}</text>']
 for i,(label,vals) in enumerate(zip(labels,values)):
  y=top+i*step+8;s.append(f'<text x="{left-12}" y="{y+12}" text-anchor="end">{html.escape(label)}</text>')
  for j,v in enumerate(vals):
   yy=y+j*12;s += [f'<rect x="{left}" y="{yy}" width="{v*scale:.3f}" height="10" fill="{colors[j]}"/>',f'<text x="{left+v*scale+5:.3f}" y="{yy+10}">{v:.3f}×</text>']
 x=left+target*scale;s.append(f'<line x1="{x}" x2="{x}" y1="{top}" y2="{height-40}" stroke="#b45309" stroke-dasharray="4,3"/>')
 s += ['</g></svg>'];(p/name).write_text('\n'.join(s)+'\n')
for kind in ['real','complex']:
 selected=sorted((r for r in rows if r['kind']==kind and r['routine']=='potrf-U'),key=key)
 draw(f'upper-speedup-{kind}.svg',f'{kind.capitalize()} Upper Cholesky: before / after', [f"{r['type']}  n={r['n']}" for r in selected],[[float(r['speedup'])] for r in selected],['Speedup (higher is better)'],1)
selected=sorted((r for r in rows if r['kind']=='complex' and r['routine']=='potrf-U'),key=key)
lower={(r['type'],r['n']):r for r in rows if r['kind']=='complex' and r['routine']=='potrf-L'}
values=[[float(r[s+'_s'])/float(lower[r['type'],r['n']][s+'_s']) for s in ['before','after']] for r in selected]
draw('upper-lower-ratio.svg','Complex Upper / Lower (lower is better; dashed target1.5×)',[f"{r['type']}  n={r['n']}" for r in selected],values,['Before4b','After4b'],1.5)
