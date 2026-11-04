#!/usr/bin/env python3
"""Interleaved CPU-library / default CUDA dispatcher / forced CUDA comparison.

Use --bin-dir build/benchmark --mpfr-bin-dir other-build/benchmark. Each
process supplies one benchmark result; alternate order each round, then
take medians. The benchmark's normal best-of-three policy for n<=2000 is
retained (one call above 2000); --repeat overrides it for cold-call studies.
Do not run correctness tests, builds, or other benchmarks concurrently.
--routines gesv posv --rhs 1 0 also measures the driver gates (0 means n).
--samples retains each call (first cold, later warm) instead of a best-of-N;
its summary reports cold and warm medians separately. CUDA cold calls include
the initial runtime/device probe. Inverses use synthetic triangular/LU factors,
not a timed preceding factorization. MPFR GPU inverses are not implemented.
The output is append-only and completed samples are skipped on restart.
Use a new output file after rebuilding or changing gates.
"""
import argparse
import csv
import io
import hashlib
import json
import os
import re
from pathlib import Path
import statistics
import subprocess
import time


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--bin-dir', type=Path, required=True)
    p.add_argument('--mpfr-bin-dir', type=Path)
    p.add_argument('--backends', nargs='+', default=['dd', 'td', 'qd', 'mpfr'])
    p.add_argument('--sizes', nargs='+', type=int, default=[512, 1024, 2048, 4096])
    p.add_argument('--routines', nargs='+', choices=['getrf', 'potrf', 'gesv', 'posv', 'trtri', 'lauum', 'potri', 'getri'], default=['getrf', 'potrf'])
    p.add_argument('--rhs', nargs='+', type=int, default=[1])
    p.add_argument('--kinds', nargs='+', choices=['real', 'complex'], default=['real', 'complex'])
    p.add_argument('--rounds', type=int, default=3)
    p.add_argument('--threads', type=int, default=32)
    p.add_argument('--repeat', type=int, default=0, help='0 keeps the benchmark repeat policy')
    p.add_argument('--upper', action='store_true')
    p.add_argument('--samples', action='store_true', help='retain cold and every warm call; default repeat=3, require >=2')
    p.add_argument('--modes', nargs='+', choices=['cpu', 'default', 'forced'], default=['cpu', 'default', 'forced'])
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    if args.rounds < 1 or args.threads < 1 or min(args.sizes) < 1 or min(args.rhs) < 0 or not 0 <= args.repeat <= 100:
        p.error('invalid round/thread/size/RHS count')
    if args.samples:
        args.repeat = args.repeat or 3
        if args.repeat < 2:
            p.error('--samples requires at least two calls')
    inverses = {'trtri', 'lauum', 'potri', 'getri'}
    triangular = {'potrf', 'posv', 'trtri', 'lauum', 'potri'}
    if 'mpfr' in args.backends and inverses.intersection(args.routines) and 'forced' in args.modes:
        p.error('MPFR inverse GPU entries do not exist; measure --modes cpu default separately')
    records = [json.loads(line) for line in args.output.read_text().splitlines()] if args.output.exists() else []
    def group_key(r):
        uplo = r.get('uplo', 'L' if r['routine'][1:] in triangular else '-')
        return tuple(r[k] for k in ('backend', 'routine', 'n', 'nrhs')) + (uplo,)
    def key(r):
        return group_key(r) + (r['round'], r['mode'], r.get('per_call', False))
    done = {key(r) for r in records if r.get('returncode') == 0}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    binary_hashes = {}
    with args.output.open('a') as output:
        for backend in args.backends:
            directory = args.mpfr_bin_dir if backend == 'mpfr' and args.mpfr_bin_dir else args.bin_dir
            for n in args.sizes:
                for routine in args.routines:
                    driver = routine in ('gesv', 'posv')
                    for selected_kind in args.kinds:
                        kind = 'R' if selected_kind == 'real' else 'C'
                        for rhs in args.rhs if driver else [-1]:
                            for round_ in range(args.rounds):
                                modes = list(args.modes)
                                if round_ % 2:
                                    modes.reverse()
                                for mode in modes:
                                    record = dict(backend=backend, routine=kind+routine, n=n,
                                                  nrhs=n if rhs == 0 else rhs, round=round_, mode=mode,
                                                  uplo=('U' if args.upper else 'L') if routine in triangular else '-',
                                                  per_call=args.samples)
                                    if key(record) in done:
                                        continue
                                    exe = directory / f'Rfactor.{backend}_{"opt" if mode == "cpu" else "cuda_total"}'
                                    command = [str(exe), '--n', str(n)]
                                    if args.repeat:
                                        command += ['--repeat', str(args.repeat)]
                                    if args.samples:
                                        command.append('--samples')
                                    if kind == 'C':
                                        command.append('--complex')
                                    if args.upper:
                                        command.append('--upper')
                                    if driver:
                                        command += ['--'+routine, '--nrhs', str(rhs)]
                                    elif routine != 'getrf':
                                        command.append('--'+routine)
                                    env = {k: v for k, v in os.environ.items() if not k.startswith('MPLAPACK_')}
                                    env.update(OMP_NUM_THREADS=str(args.threads), OMP_DYNAMIC='FALSE')
                                    env['MPLAPACK_'+backend.upper()+'_CUDA_VERBOSE'] = '1'
                                    if mode == 'forced':
                                        env['MPLAPACK_'+backend.upper()+'_LAPACK_CUDA_MIN'] = '0'
                                        if driver or routine in inverses:
                                            command.append('--entry')
                                    if exe not in binary_hashes:
                                        binary_hashes[exe] = hashlib.sha256(exe.read_bytes()).hexdigest()
                                    record.update(command=command, binary_sha256=binary_hashes[exe], threads=args.threads,
                                                  started_unix=time.time())
                                    result = subprocess.run(command, env=env, text=True, capture_output=True)
                                    record.update(returncode=result.returncode, stdout=result.stdout, stderr=result.stderr,
                                                  finished_unix=time.time())
                                    record['entry_trace_count'] = len(re.findall(r'^'+kind+routine+r' .* CUDA .*upload=', result.stderr, re.M))
                                    record['entry_trace'] = record['entry_trace_count'] > 0
                                    missing_entry = mode == 'forced' and (not record['entry_trace'] or
                                                    (args.samples and record['entry_trace_count'] != args.repeat))
                                    if missing_entry and not result.returncode:
                                        record['returncode'] = 77
                                    if not result.returncode:
                                        if args.samples:
                                            calls = list(csv.DictReader(io.StringIO(result.stdout)))
                                            if len(calls) != args.repeat or any(c['routine'] != kind+routine or int(c['n']) != n for c in calls):
                                                raise RuntimeError(f'invalid per-call output: {result.stdout}')
                                            record['calls'] = [dict(call=int(c['call']), temperature=c['temperature'], seconds=float(c['seconds'])) for c in calls]
                                            record['cold_seconds'] = record['calls'][0]['seconds']
                                            record['seconds'] = statistics.median(c['seconds'] for c in record['calls'][1:])
                                        else:
                                            record['seconds'] = float(result.stdout.strip().splitlines()[-1].split(',')[-2])
                                    output.write(json.dumps(record)+'\n')
                                    output.flush()
                                    if missing_entry and not result.returncode:
                                        raise RuntimeError(f'forced entry rejected: {command}: {result.stderr}')
                                    if result.returncode:
                                        raise RuntimeError(f'{command}: {result.returncode}: {result.stderr}')
                                    records.append(record)
                                    done.add(key(record))
                                    print(backend, kind+routine, n, rhs, round_, mode, record['seconds'], flush=True)
    records = [r for r in records if r.get('per_call', False) == args.samples]
    groups = sorted({group_key(r) for r in records})
    columns = [m+'_'+t for m in args.modes for t in ['cold', 'warm']] if args.samples else args.modes
    print(','.join(['backend', 'routine', 'n', 'nrhs', 'uplo', *columns]))
    for group in groups:
        selected = [r for r in records if group_key(r) == group and r.get('returncode') == 0]
        if args.samples:
            values = [statistics.median(r[field] for r in selected if r['mode'] == mode)
                      for mode in args.modes for field in ['cold_seconds', 'seconds']]
        else:
            values = [statistics.median(r['seconds'] for r in selected if r['mode'] == mode) for mode in args.modes]
        print(','.join(map(str, (*group, *values))))


if __name__ == '__main__':
    main()
