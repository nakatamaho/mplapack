#!/usr/bin/env python3
"""Build/run the supplemental exact GEMM check using existing project builds."""
import argparse, os, shlex, subprocess
from pathlib import Path
p=argparse.ArgumentParser()
p.add_argument('--cmake-dir',type=Path,default=Path('/tmp/factor01-cmake'))
p.add_argument('--output',type=Path,default=Path('/tmp/goal10-gemm-check'))
p.add_argument('--before-libs',type=Path)
a=p.parse_args();root=Path(__file__).resolve().parents[6];out=a.output.resolve();out.mkdir(parents=True,exist_ok=True)
cmake=a.cmake_dir.resolve();obj=out/'boundary.o'
def run(args,label,cwd=root,env=None):
    with (out/(label+'.log')).open('w') as f:
        rc=subprocess.run(args,cwd=cwd,env=env,stdout=f,stderr=subprocess.STDOUT).returncode
    (out/(label+'.exit')).write_text(str(rc)+'\n')
    if rc:raise SystemExit(f'{label}: exit {rc}; see {out}')
flags=['g++','-O2','-std=c++17','-DMPLAPACK_BUILD_WITH_MPFR',f'-I{cmake}/include','-Iinclude']
flags += [f'-Iexternal/i/{b}/include' for b in ['GMPFRXX_MKII','MPC_CUDA','MPFR','GMP','MPC']]
run(flags+['-c',str(Path(__file__).with_name('gemm-boundaries.cpp')),'-o',str(obj)],'compile')
run(['g++',str(obj),'-Lmpblas/optimized/mpfr/cuda/.libs','-Lexternal/i/MPFR/lib','-Lexternal/i/GMP/lib','-Lexternal/i/MPC/lib','-lmplapack_mpfr_opt_cuda','-lmpfr','-lmpc','-lgmp','-fopenmp','-o',str(out/'device')],'link-device')
args=shlex.split((cmake/'CMakeFiles/mpfr_cuda_Rgemm_host_test.dir/link.txt').read_text())
args=[x for x in args if not x.startswith('-Wl,--dependency-file=')]
args=[str(obj) if x.endswith('/test/Rgemm_cuda_test.cpp.o') else x for x in args]
args[args.index('-o')+1]=str(out/'host');run(args,'link-host',cmake)
env=os.environ.copy();env.update(OMP_NUM_THREADS='1',CUDA_MODULE_LOADING='EAGER')
env.pop('MPLAPACK_MPFR_CUDA_FORCE_RUNTIME',None)
ext=':'.join(str(root/f'external/i/{b}/lib') for b in ['MPFR','GMP','MPC'])
env['LD_LIBRARY_PATH']=ext;run([str(out/'host')],'host',env=env)
if a.before_libs:
    env['LD_LIBRARY_PATH']=str(a.before_libs.resolve())+':'+ext
    run([str(out/'device')],'before',env=env)
env['LD_LIBRARY_PATH']=str(root/'mpblas/optimized/mpfr/cuda/.libs')+':'+ext
run([str(out/'device')],'after',env=env)
assert (out/'host.log').read_text()==(out/'after.log').read_text()
if a.before_libs:
    assert (out/'before.log').read_text()==(out/'after.log').read_text()
