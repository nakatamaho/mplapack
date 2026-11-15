# Old (gemm-sweep, 2026-11-04) vs new (this dir) MPFR 512 Rgemm kernel-only GFLOPS.
import csv, sys, os
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
d=os.path.dirname(os.path.abspath(__file__))
def load(p):
    r=list(csv.DictReader(open(p))); return [int(x['n']) for x in r],[float(x['gflops']) for x in r]
col={'NN':'#2a6fdb','NT':'#e07b00','TN':'#1a9a5a','TT':'#c0392b'}
fig,ax=plt.subplots(figsize=(9,5.5))
for f in col:
    n,g=load(f'{d}/../gemm-sweep/mpfr512_{f}.csv'); ax.plot(n,g,'--',color=col[f],lw=1,alpha=.7,label=f'{f} before (11/04)')
    n,g=load(f'{d}/mpfr512_{f}.csv'); ax.plot(n,g,'-',color=col[f],lw=1.6,label=f'{f} now')
ax.set_xlabel('n (m=n=k)'); ax.set_ylabel('GFLOPS'); ax.grid(alpha=.3); ax.set_xlim(0,4100); ax.set_ylim(0)
ax.set_title('Rgemm MPFR 512-bit, A100, kernel only (CUDA events)\nbefore: 11/04 gemm-sweep, now: '+open(f'{d}/commit.txt').read()[:9])
ax.legend(ncol=2,fontsize=8,loc='lower right'); fig.tight_layout(); fig.savefig(sys.argv[1] if len(sys.argv)>1 else f'{d}/gemm_mpfr512_compare.png',dpi=130)
