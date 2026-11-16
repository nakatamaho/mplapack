#!/usr/bin/env python3
"""Fresh-process CPU / per-call GPU BLAS / resident-entry comparison.

Example: python3 benchmark/run_solver_bench.py --bin-dir build/benchmark \
  --backends dd qd --sizes 512 2048 --output solver.jsonl
No correctness tests or other GPU benchmarks should run concurrently.
The driver retains its usual best-of-three default; this runner explicitly
requests one sample per process and takes medians of interleaved rounds.
"""
import argparse
import json
import os
from pathlib import Path
import re
import statistics
import subprocess


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--bin-dir', type=Path, required=True)
    p.add_argument('--mpfr-bin-dir', type=Path)
    p.add_argument('--backends', nargs='+', default=['dd', 'td', 'qd', 'mpfr'])
    p.add_argument('--sizes', nargs='+', type=int, default=[512, 1024, 2048])
    p.add_argument('--routines', nargs='+', choices=['gesv', 'posv', 'getrs', 'potrs'], default=['gesv', 'posv'])
    p.add_argument('--rhs', nargs='+', type=int, default=[1, 0], help='0 means n')
    p.add_argument('--kinds', nargs='+', choices=['real', 'complex'], default=['real', 'complex'])
    p.add_argument('--rounds', type=int, default=3)
    p.add_argument('--threads', type=int, default=32)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    records = []
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('a') as output:
        for backend in args.backends:
            directory = args.mpfr_bin_dir if backend == 'mpfr' and args.mpfr_bin_dir else args.bin_dir
            for routine in args.routines:
                for kind in args.kinds:
                    for n in args.sizes:
                        for rhs in args.rhs:
                            for round_ in range(args.rounds):
                                modes = ['cpu', 'old', 'entry']
                                if round_ % 2:
                                    modes.reverse()
                                for mode in modes:
                                    executable = directory / f'Rfactor.{backend}_{"opt" if mode == "cpu" else "cuda_total"}'
                                    command = [str(executable), '--' + routine, '--n', str(n), '--nrhs', str(rhs), '--repeat', '1']
                                    if kind == 'complex':
                                        command.append('--complex')
                                    if mode == 'entry':
                                        command.append('--entry')
                                    elif mode == 'old':
                                        command.append('--legacy')
                                    env = {k: v for k, v in os.environ.items() if not k.startswith('MPLAPACK_')}
                                    env.update(OMP_NUM_THREADS=str(args.threads), MPFRXX_DEFAULT_PRECISION_BITS='512')
                                    prefix = 'MPLAPACK_' + backend.upper()
                                    env[prefix + '_LAPACK_CUDA'] = '0' if mode == 'old' else '1'
                                    env[prefix + '_LAPACK_CUDA_MIN'] = '0'
                                    env[prefix + '_LAPACK_SOLVE_CUDA_MIN'] = '0'
                                    env[prefix + '_LAPACK_DRIVER_CUDA_MIN'] = '0'
                                    env[prefix + '_LAPACK_SOLVE_CUDA_MIN_RHS'] = '0'
                                    env[prefix + '_CUDA_VERBOSE'] = '1' if mode == 'entry' else '0'
                                    result = subprocess.run(command, env=env, text=True, capture_output=True, check=True)
                                    row = result.stdout.strip().splitlines()[-1].split(',')
                                    record = dict(backend=backend, routine=('C' if kind == 'complex' else 'R') + routine,
                                                  n=n, nrhs=n if rhs == 0 else rhs, round=round_, mode=mode,
                                                  seconds=float(row[4]), stdout=result.stdout, stderr=result.stderr)
                                    for field in ('upload', 'factor_solve', 'download'):
                                        values = re.findall(r'\b' + field + r'=([\deE+.\-]+)', result.stderr)
                                        if values:
                                            record[field] = float(values[-1])
                                    if mode == 'entry' and 'factor_solve' not in record:
                                        raise RuntimeError('missing direct-entry timing trace')
                                    output.write(json.dumps(record) + '\n')
                                    output.flush()
                                    records.append(record)
                                    print(f'{backend} {record["routine"]} n={n} rhs={record["nrhs"]} round={round_} {mode} {record["seconds"]:.6f}', flush=True)
    keys = sorted({(r['backend'], r['routine'], r['n'], r['nrhs']) for r in records})
    print('backend,routine,n,nrhs,cpu,old,entry,upload,factor_solve,download')
    for key in keys:
        selected = [r for r in records if (r['backend'], r['routine'], r['n'], r['nrhs']) == key]
        medians = [statistics.median(r['seconds'] for r in selected if r['mode'] == mode) for mode in ('cpu', 'old', 'entry')]
        split = [statistics.median(r[field] for r in selected if r['mode'] == 'entry') for field in ('upload', 'factor_solve', 'download')]
        print(','.join(map(str, (*key, *medians, *split))))


if __name__ == '__main__':
    main()
