# MPFR 512 Rgemm kernel-only GFLOPS: V100 (this dir) vs A100 (mpfr512-sweep-2).
import csv, sys, os
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
d=os.path.dirname(os.path.abspath(__file__))
a=f'{d}/../../x86_64-pc-linux-gnu-Ryzen-Threadripper-3970X-TeslaA100/mpfr512-sweep-2'
def load(p):
    r=list(csv.DictReader(open(p))); return [int(x['n']) for x in r],[float(x['gflops']) for x in r]
col={'NN':'#2a6fdb','NT':'#e07b00','TN':'#1a9a5a','TT':'#c0392b'}
fig,ax=plt.subplots(figsize=(9,5.5))
for f in col:
    n,g=load(f'{a}/mpfr512_{f}.csv'); ax.plot(n,g,'--',color=col[f],lw=1,alpha=.7,label=f'{f} A100')
    n,g=load(f'{d}/mpfr512_{f}.csv'); ax.plot(n,g,'-',color=col[f],lw=1.6,label=f'{f} V100')
ax.set_xlabel('n (m=n=k)'); ax.set_ylabel('GFLOPS'); ax.grid(alpha=.3); ax.set_xlim(0,4100); ax.set_ylim(0)
ax.set_title('Rgemm MPFR 512-bit, kernel only (CUDA events)\nV100 (sm_70, CUDA 12.6): '+open(f'{d}/commit.txt').read()[:9]+', A100: mpfr512-sweep-2')
ax.legend(ncol=2,fontsize=8,loc='lower right'); fig.tight_layout(); fig.savefig(f'{d}/gemm_mpfr512_v100_a100.png',dpi=130)
