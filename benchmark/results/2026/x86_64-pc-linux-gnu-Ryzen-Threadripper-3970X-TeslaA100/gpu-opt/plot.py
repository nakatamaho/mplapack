#!/usr/bin/env python3
"""Regenerate the SVG figures using only the Python standard library."""
import csv,math,html
from pathlib import Path
p=Path(__file__).resolve().parent
colors={'dd':'#0072b2','td':'#e69f00','qd':'#009e73','mpfr512':'#cc79a7'}
def rows(n):return list(csv.DictReader((p/n).open()))
def text(x,y,s,size=13,anchor='start'):
    return f'<text x="{x:.2f}" y="{y:.2f}" font-size="{size}" text-anchor="{anchor}">{html.escape(str(s))}</text>'
def line(x1,y1,x2,y2,color='#bbb',dash=''):
    return f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" y2="{y2:.2f}" stroke="{color}" stroke-dasharray="{dash}"/>'
def save(name,w,h,parts):
    (p/(name+'.svg')).write_text(f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}"><rect width="100%" height="100%" fill="white"/><g font-family="sans-serif" fill="#222">'+''.join(parts)+'</g></svg>\n')
def chart(parts,x,y,w,h,title,series,xs,ymin,ymax,logy=True,model=None):
    parts.append(text(x+w/2,y-16,title,16,'middle'))
    sx=lambda v:x+w*(math.log2(v)-math.log2(min(xs)))/(math.log2(max(xs))-math.log2(min(xs)))
    sy=lambda v:y+h-h*((math.log2(v)-math.log2(ymin))/(math.log2(ymax)-math.log2(ymin)) if logy else (v-ymin)/(ymax-ymin))
    ticks=[2**i for i in range(round(math.log2(ymin)),round(math.log2(ymax))+1)] if logy else [ymax*i/4 for i in range(5)]
    for v in ticks:
        parts.extend([line(x,sy(v),x+w,sy(v),'#444' if v==1 and logy else '#ddd'),text(x-8,sy(v)+4,f'{v:g}',11,'end')])
    for v in xs:parts.extend([line(sx(v),y,sx(v),y+h,'#eee'),text(sx(v),y+h+20,v,11,'middle')])
    parts.extend([line(x,y,x,y+h,'#333'),line(x,y+h,x+w,y+h,'#333'),text(x+w/2,y+h+39,'n',12,'middle')])
    if model:parts.append(line(x,sy(model),x+w,sy(model),'#777','3 4'))
    for label,color,dash,data in series:
        pts=' '.join(f'{sx(a):.2f},{sy(b):.2f}' for a,b in data)
        parts.append(f'<polyline points="{pts}" fill="none" stroke="{color}" stroke-width="2" stroke-dasharray="{dash}"/>')
        for a,b in data:parts.append(f'<circle cx="{sx(a):.2f}" cy="{sy(b):.2f}" r="3" fill="{color}"/>')
factors=rows('factors.csv')
for metric,title,name in [('gpu_speedup','Original GPU / final GPU','factor-speedup'),('gpu_vs_cpu','CPU32 / final GPU','gpu-vs-cpu')]:
    parts=[text(700,26,title+' (ratio > 1 favors final GPU)',20,'middle')]
    for i,(b,c) in enumerate(colors.items()):parts.extend([line(420+i*150,49,450+i*150,49,c),text(457+i*150,54,b)])
    vals=[float(r[metric]) for r in factors];lo=2**math.floor(math.log2(min(vals+[1])*.9));hi=2**math.ceil(math.log2(max(vals+[1])*1.1))
    for row,kind in enumerate(['real','complex']):
        for col,op in enumerate(['getrf','potrf-L','potrf-U']):
            series=[]
            for b,c in colors.items():
                data=sorted((int(r['n']),float(r[metric])) for r in factors if r['type']==b and r['kind']==kind and r['routine']==op)
                series.append((b,c,'',data))
            chart(parts,70+col*460,105+row*335,365,250,kind+' '+op,series,[1024,2048,4096],lo,hi)
    save(name,1400,760,parts)
rs=rows('kernel-comparison.csv');forms=list(dict.fromkeys((r['kernel'],r['form']) for r in rs if r['type']=='dd' and r['n']=='1024'))
parts=[text(750,28,'Kernel benchmark: n=1024, original / final GPU time',20,'middle')];lo=.5;hi=128;y=75;h=320
sy=lambda v:y+h-h*(math.log2(v)-math.log2(lo))/(math.log2(hi)-math.log2(lo))
for v in [.5,1,2,4,8,16,32,64,128]:parts.extend([line(55,sy(v),1470,sy(v),'#555' if v==1 else '#ddd'),text(48,sy(v)+4,v,11,'end')])
for j,b in enumerate(['dd','qd','mpfr512']):
    data={(r['kernel'],r['form']):float(r['speedup']) for r in rs if r['type']==b and r['n']=='1024'}
    parts.extend([line(500+j*170,52,530+j*170,52,colors[b]),text(535+j*170,57,b)])
    for i,f in enumerate(forms):
        x=65+i*94+j*24;v=data[f];parts.append(f'<rect x="{x}" y="{sy(v):.2f}" width="22" height="{sy(lo)-sy(v):.2f}" fill="{colors[b]}"/>')
for i,(k,f) in enumerate(forms):parts.extend([text(100+i*94,417,k,11,'middle'),text(100+i*94,433,f,10,'middle')])
parts.append(text(750,466,'MPFR potrf2 above nb measures blocked Rpotrf; BF potrf2 is recursive/unblocked.',12,'middle'));save('kernel-speedup',1500,485,parts)
rs=rows('mpfr-gemm.csv');peaks={r['type']:2*float(r['ops_per_second'])/1e9 for r in csv.DictReader((p.parent/'roofline/peaks.csv').open()) if r['op']=='madd'}
parts=[text(700,25,'MPFR GEMM: nominal GFLOP/s (kernel timing)',20,'middle')]
for i,(label,c,d) in enumerate([('NN before','#0072b2','6 4'),('NN after','#0072b2',''),('TN before','#d55e00','6 4'),('TN after','#d55e00','')]):parts.extend([line(270+i*210,49,300+i*210,49,c,d),text(307+i*210,54,label)])
for col,prec in enumerate([256,512,1024]):
    series=[];xs=set()
    for form,c in [('NN','#0072b2'),('TN','#d55e00')]:
        for stage,dash in [('before','6 4'),('after','')]:
            data=sorted((int(r['n']),float(r['gflops'])) for r in rs if r['type']==f'mpfr{prec}' and r['form']==form and r['stage']==stage);xs.update(n for n,_ in data);series.append((form+' '+stage,c,dash,data))
    peak=peaks[f'mpfr{prec}'];chart(parts,70+460*col,105,365,250,f'{prec} bits',series,sorted(xs),0,math.ceil(peak*1.1),False,peak)
parts.append(text(700,420,'Dotted line: full-width mixed-path cu_freal model, not measured kernel occupancy/stalls.',12,'middle'));save('mpfr-gemm',1400,440,parts)
