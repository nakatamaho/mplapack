#!/usr/bin/env python3
"""Quick dispatch trace checks, not a performance measurement or numerical QA."""
import argparse
import os
from pathlib import Path
import re
import subprocess

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--bin-dir', type=Path, required=True)
p.add_argument('--mpfr-bin-dir', type=Path, required=True)
a = p.parse_args()
checks = 0

def check(backend, routine, n, expected, upper=False, settings=None):
    global checks
    directory = a.mpfr_bin_dir if backend == 'mpfr' else a.bin_dir
    cmd = [str(directory / f'Rfactor.{backend}_cuda_total'), '--n', str(n), '--repeat', '1']
    if routine[0] == 'C':
        cmd.append('--complex')
    if routine[1:] != 'getrf':
        cmd.append('--' + routine[1:])
    if upper:
        cmd.append('--upper')
    env = {k: v for k, v in os.environ.items() if not k.startswith('MPLAPACK_')}
    prefix = 'MPLAPACK_' + backend.upper()
    env.update(OMP_NUM_THREADS='32', OMP_DYNAMIC='FALSE')
    env[prefix + '_CUDA_VERBOSE'] = '1'
    for key, value in (settings or {}).items():
        env[prefix + key] = str(value)
    r = subprocess.run(cmd, env=env, text=True, capture_output=True, timeout=180)
    traces = len(re.findall(r'^' + routine + r' .* CUDA .*upload=', r.stderr, re.M))
    assert r.returncode == 0 and traces == int(expected), (cmd, settings, r.returncode, traces, r.stderr)
    checks += 1
    print(f'PASS {backend} {routine} n={n} upper={upper} settings={settings} entry={traces}', flush=True)

for backend in ('td', 'qd'):
    for upper in (False, True):
        check(backend, 'Rpotrf', 256, False, upper)
        check(backend, 'Rpotrf', 512, True, upper)
for routine in ('Rgetrf', 'Cgetrf'):
    check('qd', routine, 256, False)
    check('qd', routine, 512, True)
for routine in ('Rgetrf', 'Rgesv'):
    check('mpfr', routine, 512, False)
    check('mpfr', routine, 1024, True)
check('mpfr', 'Cgesv', 1024, False)
for backend in ('dd', 'td', 'qd'):
    check(backend, 'Cpotrf', 4096, False, True)
for backend in ('dd', 'td', 'qd', 'mpfr'):
    check(backend, 'Rpotrf', 32, True, settings={'_LAPACK_CUDA_MIN': 0})
    check(backend, 'Rpotrf', 32, False, settings={'_LAPACK_CUDA_MIN': 0, '_LAPACK_CUDA': 0})
    check(backend, 'Cpotrf', 32, True, settings={'_LAPACK_CUDA_MIN': 0})
    check(backend, 'Rgesv', 32, True, settings={'_LAPACK_DRIVER_CUDA_MIN': 0})
print(f'PASS {checks} dispatch checks')
