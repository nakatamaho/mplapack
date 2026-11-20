# MPFR 512 real Rgetrf / Rpotrf L / Rpotrf U through the GPU entries (whole entry:
# packing + transfers + device factorization): V100 (this dir) vs A100 (mpfr512-sweep-2).
import csv, sys, os
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
d=os.path.dirname(os.path.abspath(__file__))
a=f'{d}/../../x86_64-pc-linux-gnu-Ryzen-Threadripper-3970X-TeslaA100/mpfr512-sweep-2'
col={'getrf':'#2a6fdb','potrf-L':'#1a9a5a','potrf-U':'#c0392b'}
files={'getrf':'getrf.csv','potrf-L':'potrf_L.csv','potrf-U':'potrf_U.csv'}
def load(p): return dict((int(x['n']),float(x['gflops'])) for x in csv.DictReader(open(p)))
fig,ax=plt.subplots(figsize=(9,5.5))
for r,f in files.items():
    A=load(f'{a}/{f}'); V=load(f'{d}/{f}')
    ax.plot(list(A),list(A.values()),'--',color=col[r],lw=1,alpha=.7,label=f'R{r} A100')
    ax.plot(list(V),list(V.values()),'-',color=col[r],lw=1.5,label=f'R{r} V100')
    for n in (512,1024,2048,4096): print(f'{r:8s} {n:5d} V100={V[n]:.3f} A100={A[n]:.3f} x{V[n]/A[n]:.2f}')
ax.set_xlabel('n'); ax.set_ylabel('GFLOPS (getrf 2n³/3, potrf n³/3)'); ax.grid(alpha=.3); ax.set_xlim(0,4100); ax.set_ylim(0)
ax.set_title('MPFR 512-bit real factorizations, GPU entry incl. transfers\nV100 (sm_70, CUDA 12.6): '+open(f'{d}/commit.txt').read()[:9]+', A100: mpfr512-sweep-2')
ax.legend(fontsize=8,loc='upper left'); fig.tight_layout(); fig.savefig(f'{d}/factor_mpfr512_v100_a100.png',dpi=130)
