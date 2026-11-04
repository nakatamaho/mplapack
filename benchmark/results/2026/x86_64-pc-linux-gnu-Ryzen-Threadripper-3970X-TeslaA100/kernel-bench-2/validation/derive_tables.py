#!/usr/bin/env python3
"""Derive tables from actual accepted CSV rows; no benchmark execution."""
import collections,csv,json,pathlib,sys
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[5]))
import kernel_new as kn
out=pathlib.Path(__file__).resolve().parents[1]
rows=kn.read_rows(out/'kernels.csv')
accepted=[r for r in rows if r['phase']=='sweep' and r['seconds'] and r['exit']=='0']
index={}
for r in accepted:
    key=(r['type'],r['kernel'],r['method'],int(r['n']))
    assert key not in index,key
    assert float(r['seconds'])>0,key
    index[key]=r

def save(name,records):
    with (out/name).open('w',newline='') as f:
        w=csv.DictWriter(f,records[0].keys());w.writeheader();w.writerows(records)

coverage=[]
for typ,k,m in kn.series(['dd','td','qd','mpfr512'],['gpu','cpu']):
    rr=[r for r in accepted if (r['type'],r['kernel'],r['method'])==(typ,k,m)]
    nn=sorted(int(r['n']) for r in rr)
    assert nn,(typ,k,m)
    baseline=kn.grid(k)
    dense=kn.grid(k,True)
    coverage.append(dict(type=typ,kernel=k,method=m,points=len(nn),min_n=min(nn),max_n=max(nn),
                         baseline_points=len(baseline),dense_points=len(dense),
                         missing_baseline=';'.join(str(n) for n in baseline if n not in nn),
                         unmeasured_dense=len(set(dense)-set(nn)),
                         reported_over10s=';'.join(r['n'] for r in rr if float(r['seconds'])>10),
                         loaded=sum(r['loaded']=='1' for r in rr),gpu_busy=sum(r['gpu_busy']=='1' for r in rr),
                         recorded_wall_seconds=sum(float(r['elapsed']) for r in rr),
                         measured_n=';'.join(map(str,nn))))
save('coverage.csv',coverage)
ratios=[]
for typ,k in sorted({(r['type'],r['kernel']) for r in accepted}):
    gn={n for t,kk,m,n in index if (t,kk,m)==(typ,k,'gpu')}
    cn={n for t,kk,m,n in index if (t,kk,m)==(typ,k,'cpu')}
    common=gn&cn
    assert common,(typ,k)
    powers={n for n in common if n&(n-1)==0}
    n=next((n for n in (2048,1024) if n in common),max(powers or common))
    choices=[('representative',n)]
    if k in kn.NAMES[:11] and 1<<21 in common:
        choices.append(('vector_2^21',1<<21))
    for purpose,n in choices:
        g=index[typ,k,'gpu',n];c=index[typ,k,'cpu',n];d=index.get((typ,k,'cpu-default',n))
        ge=float(g['seconds']);ce=float(c['seconds'])
        refs=[g,c]+([d] if d else [])
        ratios.append(dict(type=typ,kernel=k,purpose=purpose,n=n,gpu_seconds=g['seconds'],cpu_seconds=c['seconds'],
                           cpu_default_seconds=d['seconds'] if d else '',
                           gpu_speedup=ce/ge,gpu_speedup_vs_default=float(d['seconds'])/ge if d else '',
                           gpu_rate=(n/ge*1e-9 if k in kn.NAMES[1:4] else float(g['gflops'])),
                           rate_unit='G elements/s' if k in kn.NAMES[1:4] else 'GFLOPS',
                           flagged=int(any(r['loaded']=='1' or r['gpu_busy']=='1' for r in refs)),
                           graph=f'{typ}/{k}.png'))
save('ratios.csv',ratios)
extra=json.loads((out/'budget-overhead.json').read_text())['seconds']
checks={p.parent.name:dict(count=len(rr:=kn.read_rows(p)),statuses=dict(collections.Counter(r['status'] for r in rr)),
                          max_absolute_error=max(float(r['max_absolute_error']) for r in rr))
        for p in sorted(out.glob('checks*/checks.csv'))}
summary=dict(records=len(rows),phases=dict(collections.Counter(r['phase'] for r in rows)),
             accepted_sweep_points=len(accepted),series=len(coverage),graphs=len({(r['type'],r['kernel']) for r in accepted}),
             failed_records=sum(not r['seconds'] or r['exit']!='0' for r in rows),
             recorded_wall_seconds=sum(float(r['elapsed']) for r in rows),extra_charged_seconds=extra,
             charged_measurement_seconds=sum(float(r['elapsed']) for r in rows)+extra,
             accepted_load_min=min(float(r['load_before']) for r in accepted),
             accepted_load_max=max(float(r['load_max']) for r in accepted),
             accepted_loaded_points=sum(r['loaded']=='1' for r in accepted),
             accepted_other_gpu_points=sum(r['gpu_busy']=='1' for r in accepted),
             accepted_gpu_points=sum(r['method']=='gpu' for r in accepted),checks=checks)
(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
primary=[r for r in ratios if r['purpose']=='representative']
md=['# Representative measured times and graph index', '',
    'CPU/GPU is CPU seconds divided by GPU seconds: >1 means faster GPU.',
    'CPU is FAST for the seven reductions, default otherwise; OMP_NUM_THREADS=32.',
    'Prefer n=2048, then1024, otherwise the largest common measured power of two.',
    '† marks a sampled load>=16 or other-GPU-app flag. See README for scope and omissions.',
    'The CSV also contains the CPU-default reduction times and vector n=2^21 comparisons.', '',
    '| Type | Routine | n | GPU ms | CPU ms | CPU/GPU | Flags | Graph |',
    '|---|---|---:|---:|---:|---:|---|---|']
for r in primary:
    md.append(f"| {r['type']} | {r['kernel']} | {r['n']} | {float(r['gpu_seconds'])*1000:.4g} | {float(r['cpu_seconds'])*1000:.4g} | {r['gpu_speedup']:.3g}× | {'†' if r['flagged'] else ''} | [PNG]({r['graph']}) |")
(out/'ratios.md').write_text('\n'.join(md)+'\n')
md=['# Additional omissions from the reduced23-point matrix /57-point vector grids', '',
    'All other series completed their baseline grid. Dense matrix grids would have712 points;',
    'the small powers and interior multiples were deliberately reduced before the baseline.',
    'Remaining MPFR omissions below follow reported>10s cutoffs; see coverage.csv and reductions.log.', '',
    '| Type | Routine | Method | Points | Max n | Omitted baseline sizes |', '|---|---|---|---:|---:|---|']
for r in coverage:
    if r['missing_baseline']:
        md.append(f"| {r['type']} | {r['kernel']} | {r['method']} | {r['points']} | {r['max_n']} | {r['missing_baseline']} |")
(out/'grid-omissions.md').write_text('\n'.join(md)+'\n')
print(json.dumps(summary,indent=2))
print('SLOWEST REPRESENTATIVES')
for r in sorted(primary,key=lambda r:r['gpu_speedup'])[:16]:print(r['type'],r['kernel'],r['n'],r['gpu_speedup'],r['flagged'])
print('VECTOR 2^21 RANGES')
for t in ['dd','td','qd','mpfr512']:
    rr=[r for r in ratios if r['type']==t and r['purpose']=='vector_2^21']
    print(t,min(r['gpu_speedup'] for r in rr),max(r['gpu_speedup'] for r in rr))
