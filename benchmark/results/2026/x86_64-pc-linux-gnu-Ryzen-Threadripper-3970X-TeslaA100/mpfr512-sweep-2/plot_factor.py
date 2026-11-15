# MPFR 512 real Rgetrf / Rpotrf L / Rpotrf U through the GPU entries (whole entry:
# packing + transfers + device factorization). Points: goal 10 (2026-11-07,
# ../gpu-opt/upper/matrix.csv after_s, same entry, warmup 1, best of 3).
import csv, sys, os
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
d=os.path.dirname(os.path.abspath(__file__))
def flops(r,n): return (2/3 if r=='getrf' else 1/3)*n**3
old={}
for x in csv.DictReader(open(f'{d}/../gpu-opt/upper/matrix.csv')):
    if x['type']=='mpfr512' and x['kind']=='real':
        old.setdefault(x['routine'],{})[int(x['n'])]=flops(x['routine'],int(x['n']))/float(x['after_s'])/1e9
col={'getrf':'#2a6fdb','potrf-L':'#1a9a5a','potrf-U':'#c0392b'}
files={'getrf':'getrf.csv','potrf-L':'potrf_L.csv','potrf-U':'potrf_U.csv'}
fig,ax=plt.subplots(figsize=(9,5.5)); rows=[]
for r,f in files.items():
    data=[(int(x['n']),float(x['gflops'])) for x in csv.DictReader(open(f'{d}/{f}'))]
    if data: ax.plot(*zip(*data),'-',color=col[r],lw=1.5,label=f'R{r} now')
    if r in old: ax.plot(list(old[r]),list(old[r].values()),'o',mfc='none',color=col[r],ms=8,label=f'R{r} goal 10 (11/07)')
    nd=dict(data)
    for n in (512,1024,2048,4096): rows.append((r,n,old.get(r,{}).get(n),nd.get(n)))
ax.set_xlabel('n'); ax.set_ylabel('GFLOPS (getrf 2n³/3, potrf n³/3)'); ax.grid(alpha=.3); ax.set_xlim(0,4100); ax.set_ylim(0)
ax.set_title('MPFR 512-bit real factorizations, A100, GPU entry incl. transfers\nnow: '+open(f'{d}/commit.txt').read()[:9])
ax.legend(fontsize=8,loc='upper left'); fig.tight_layout(); fig.savefig(f'{d}/factor_mpfr512.png',dpi=130)
for r,n,o,nw in rows: print(f'{r:8s} {n:5d} before={o and round(o,3)} now={nw and round(nw,3)} '+(f'x{nw/o:.2f}' if o and nw else ''))
