#!/usr/bin/env python3
"""Reproduce the supplementary resident/public checks on an otherwise idle GPU."""
import argparse,filecmp,hashlib,json,os,subprocess
from pathlib import Path
ap=argparse.ArgumentParser()
for name in ['source-before','source-after','build','libs-before','libs-after','output']:
 ap.add_argument('--'+name,type=Path,required=True)
ap.add_argument('--dependencies',type=Path)
ap.add_argument('--nvcc',default='/usr/local/cuda/bin/nvcc')
ap.add_argument('--arch',default='sm_80')
ap.add_argument('--timings',action='store_true')
ap.add_argument('--memcheck-updates',action='store_true')
a=ap.parse_args();here=Path(__file__).resolve().parent
for name in ['source_before','source_after','build','libs_before','libs_after','output']:
 setattr(a,name,getattr(a,name).resolve())
deps=(a.dependencies or a.source_after/'external/i').resolve()
a.output.mkdir(parents=True,exist_ok=True)
env=dict(os.environ,OMP_NUM_THREADS='1',CUDA_MODULE_LOADING='EAGER')
def run(name,cmd,extra=None,gpu=False):
 e=dict(env);e.update(extra or {})
 if gpu:
  q=subprocess.check_output(['nvidia-smi','--query-compute-apps=pid,process_name','--format=csv,noheader'],text=True)
  (a.output/(name+'.before')).write_text(q);assert not q.strip(),'GPU is not idle'
 with (a.output/(name+'.log')).open('w') as f:r=subprocess.run([str(x) for x in cmd],env=e,stdout=f,stderr=subprocess.STDOUT)
 (a.output/(name+'.exit')).write_text(str(r.returncode)+'\n')
 if gpu:
  q=subprocess.check_output(['nvidia-smi','--query-compute-apps=pid,process_name','--format=csv,noheader'],text=True)
  (a.output/(name+'.after')).write_text(q);assert not q.strip(),'foreign GPU process appeared'
 assert r.returncode==0,(name,r.returncode)
cu=a.build/'external/i/MPC_CUDA'
for stage,src in [('before',a.source_before),('after',a.source_after)]:
 cmd=[a.nvcc,'-U_GNU_SOURCE','-D_XOPEN_SOURCE=700','-std=c++17','-O3','-arch='+a.arch,'-rdc=true','-fmad=false','-DMPLAPACK_BUILD_WITH_MPFR']
 if stage=='after':cmd+=['-DUPPER_AFTER']
 cmd+=['-I'+str(src/'mpblas/optimized/mpfr/cuda'),'-I'+str(cu/'include'),'-I'+str(a.build/'include'),'-I'+str(a.source_after/'include'),here/'resident-dots.cu',cu/'lib/libmpc_cuda.a','-o',a.output/('resident-'+stage)]
 run('build-resident-'+stage,cmd)
# The reference routine is renamed in the helper, not obtained through the
# accelerator's shadowing dispatcher. Generated reference sources are unchanged.
assert (a.source_before/'mpblas/reference/Cherk.cpp').read_bytes()==(a.source_after/'mpblas/reference/Cherk.cpp').read_bytes()
for name in ['herk-boundaries','trsm-canary','ctrsm-canary']:
 cmd=['c++','-O2','-std=c++17','-DMPLAPACK_BUILD_WITH_MPFR','-I'+str(a.build/'include'),'-I'+str(a.source_after),'-I'+str(a.source_after/'include')]
 cmd+=['-I'+str(deps/d/'include') for d in ['GMPFRXX_MKII','MPFR','GMP','MPC']]
 cmd+=[here/(name+'.cpp'),'-L'+str(a.libs_after),'-l:libmplapack_mpfr_opt_cuda.so.3']
 cmd+=['-L'+str(deps/d/'lib') for d in ['MPC','MPFR','GMP']]
 cmd+=['-lmpc','-lmpfr','-lgmp','-o',a.output/name];run('build-'+name,cmd)
checks=[]
for rank in [0,1]:
 for n in [7,63,64,65]:
  for stage in ['before','after']:
   name=f'boundary-{rank}-{n}-{stage}'
   run(name,[a.output/('resident-'+stage),rank,n],{'UPPER_DUMP':str(a.output/(name+'.dump'))},True)
  b=a.output/f'boundary-{rank}-{n}-before.dump';f=a.output/f'boundary-{rank}-{n}-after.dump'
  assert filecmp.cmp(b,f,shallow=False),(rank,n)
  checks.append({'rank':rank,'n':n,'sha256':hashlib.sha256(b.read_bytes()).hexdigest()})
for stage in ['before','after']:
 lib=getattr(a,'libs_'+stage)
 ld=':'.join([str(lib)]+[str(deps/d/'lib') for d in ['MPC','MPFR','GMP']])
 run('herk-'+stage,[a.output/'herk-boundaries'],{'LD_LIBRARY_PATH':ld},True)
 if stage=='after':
  for name in ['trsm-canary','ctrsm-canary']:run(name,[a.output/name],{'LD_LIBRARY_PATH':ld},True)
if a.timings:
 for rank in [0,1]:
  for n in [512,1024,2048]:
   for stage in ['before','after']:run(f'timing-{rank}-{n}-{stage}',[a.output/('resident-'+stage),rank,n],gpu=True)
if a.memcheck_updates:
 for rank in [0,1]:run(f'memcheck-update-{rank}',['/usr/local/cuda/bin/compute-sanitizer','--tool','memcheck','--error-exitcode','99',a.output/'resident-after',rank,65],gpu=True)
(a.output/'boundary-checks.json').write_text(json.dumps(checks,indent=2)+'\n')
print('PASS: eight exact resident boundary dumps; public HERK and TRSM checks. Full MPFR complex memcheck remains UNVERIFIED.')
