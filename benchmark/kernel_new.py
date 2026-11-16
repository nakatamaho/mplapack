#!/usr/bin/env python3
"""Goal11 only: called by kernel_bench.py new. Sequential, resumable sweeps.
No legacy series or measurement policy is changed. GPU runs require the
maintainer's marker; this script never creates it. See README.kernel_bench.md.
"""
import argparse
import csv
import json
import math
import os
from pathlib import Path
import subprocess
import time

NAMES = ['Raxpy', 'Rcopy', 'iRamax', 'iCamax', 'Rdot', 'Cdotc', 'Cdotu',
         'Rasum', 'RCasum', 'Rnrm2', 'RCnrm2', 'Rtrsv', 'Rtbsv', 'Rtpsv',
         'Ctrsv', 'Ctbsv', 'Ctpsv', 'Rtrtri', 'Rlauum', 'Rpotri', 'Rgetri',
         'Ctrtri', 'Clauum', 'Cpotri', 'Cgetri']
SUMS = set(NAMES[4:11])
FIELDS = ['type', 'kernel', 'form', 'n', 'method', 'seconds', 'gflops', 'runs',
          'load_before', 'load_after', 'load_max', 'loaded', 'gpu_busy', 'elapsed', 'exit', 'phase']
GATE = Path('/tmp/gemm-sweep-forms.done')
ROOT = Path(__file__).resolve().parent.parent


def grid(kernel, dense=False):
    if kernel in NAMES[:11]:
        return sorted({(1 << e) + q * (1 << e) // 4 for e in range(10, 24) for q in range(4)} | {1 << 24})
    if dense:
        return sorted(set(range(7, 4097, 7)) | set(range(31, 4097, 31)) | {1 << e for e in range(13)})
    # Reduced grid retains all three requested families, including1024/2048.
    return sorted({7, 31} | {1 << e for e in range(6, 13)} |
                  {(1 << e) // k * k for e in range(6, 13) for k in (7, 31)})


def series(types, methods):
    for typ in types:
        for kernel in NAMES if typ == 'mpfr512' else NAMES[:17]:
            for method in methods:
                yield typ, kernel, method
                if method == 'cpu' and kernel in SUMS:
                    yield typ, kernel, 'cpu-default'


def read_rows(path):
    with path.open(newline='') as f:
        return list(csv.DictReader(f))


def append_csv(path, row):
    exists = path.exists()
    with path.open('a') as f:
        w = csv.DictWriter(f, FIELDS)
        if not exists:
            w.writeheader()
        w.writerow(row)


def snapshot(path, tag, gpu=False, own_pid=None):
    load = os.getloadavg()[0]
    commands = [['uptime']]
    if gpu:
        commands.append(['nvidia-smi', '--query-compute-apps=pid,process_name', '--format=csv,noheader'])
    busy = False
    with path.open('a') as f:
        f.write(time.strftime('%Y-%m-%dT%H:%M:%S') + ' ' + tag + '\n')
        for cmd in commands:
            p = subprocess.run(cmd, capture_output=True, text=True)
            f.write('$ ' + ' '.join(cmd) + '\n' + p.stdout + p.stderr + f'exit={p.returncode}\n')
            if cmd[0] == 'nvidia-smi':
                # Compare PIDs, not names: another user's job can be the same
                # benchmark executable. Before launch, every compute app is other;
                # afterwards exclude our just-reaped PID if NVML still lists it.
                busy = p.returncode != 0 or any(line.strip() and line.split(',')[0].strip() != str(own_pid)
                                               for line in p.stdout.splitlines())
    return load, busy


class Runner:
    def __init__(self, a):
        self.a = a
        self.out = Path(a.out)
        (self.out / 'raw').mkdir(parents=True, exist_ok=True)
        self.csv = self.out / 'kernels.csv'
        self.rows = read_rows(self.csv) if self.csv.exists() else []
        # Cumulative point wall time includes setup and condition probes, not
        # merely best kernel times. Optional receipts charge interrupted/older
        # runs' otherwise-unrecorded wall time; restarting cannot reset the cap.
        self.spent = sum(float(r['elapsed']) for r in self.rows)
        overhead = self.out / 'budget-overhead.json'
        if overhead.exists():
            self.spent += float(json.loads(overhead.read_text())['seconds'])
        self.limit = a.hours * 3600
        self.env = {k: v for k, v in os.environ.items() if not k.startswith(('MPLAPACK_', 'MPFRXX_', 'OMP_'))}
        self.env['LD_LIBRARY_PATH'] = ':'.join(str(ROOT / 'external/i' / x / 'lib') for x in ('QD', 'MPFR', 'GMP', 'MPC')) + ':/usr/local/cuda/lib64:' + os.environ.get('LD_LIBRARY_PATH', '')

    def measure(self, typ, kernel, method, n, phase):
        matches = [r for r in self.rows if (r['type'], r['kernel'], r['method'], int(r['n']), r['phase']) == (typ, kernel, method, n, phase)]
        if matches:
            return matches[-1]
        if self.spent >= self.limit:
            raise RuntimeError('cumulative measurement budget exhausted')
        gpu = method == 'gpu'
        if gpu:
            while not GATE.is_file():
                with (self.out / 'gate.log').open('a') as f:
                    f.write(time.strftime('%Y-%m-%dT%H:%M:%S') + ' marker absent; sleeping120; no GPU execution\n')
                time.sleep(120)
        b = 'mpfr' if typ == 'mpfr512' else typ
        bindir = self.a.mpfr_bindir if b == 'mpfr' else self.a.bindir
        suffix = ('cuda_resident' if b == 'mpfr' else 'cuda_kernel') if gpu else 'opt'
        cmd = [str(Path(bindir).resolve() / f'Rkernel_bench.{b}_{suffix}'), kernel,
               'default' if method == 'cpu-default' else 'fast', str(n)]
        if phase == 'pilot':
            cmd += ['--reps', '1', '--warmup', '1']
        env = dict(self.env, OMP_NUM_THREADS='1' if gpu else '32', OMP_DYNAMIC='FALSE')
        if b == 'mpfr':
            env['MPFRXX_DEFAULT_PRECISION_BITS'] = '512'
        tag = f'{typ}_{kernel}_{method}_{n}_{phase}'
        log = self.out / 'raw' / (tag + '.log')
        start = time.monotonic()
        before, busy = snapshot(self.out / 'conditions.log', 'before ' + tag, gpu)
        loads = [before]
        process_start = time.monotonic()
        with log.open('w') as f:
            f.write('$ ' + ' '.join(cmd) + '\nOMP_NUM_THREADS=' + env['OMP_NUM_THREADS'] + '\n')
            f.flush()
            proc = subprocess.Popen(cmd, env=env, stdout=subprocess.PIPE, stderr=f, text=True)
            timed_out = False
            try:
                while True:
                    remaining = self.limit - self.spent - (time.monotonic() - start)
                    if remaining <= 0:
                        proc.kill()
                        text, _ = proc.communicate()
                        timed_out = True
                        break
                    try:
                        text, _ = proc.communicate(timeout=min(30, remaining))
                        break
                    except subprocess.TimeoutExpired:
                        load, active = snapshot(self.out / 'conditions.log', 'during ' + tag, gpu, proc.pid)
                        loads.append(load)
                        busy |= active
            finally:
                if proc.poll() is None:
                    proc.kill()
                    proc.wait()
            process_elapsed = time.monotonic() - process_start
            f.write(text + f'\nexit={proc.returncode}, process_elapsed={process_elapsed:.6f}\n')
        (log.with_suffix('.log.exit')).write_text(str(proc.returncode) + '\n')
        after, active = snapshot(self.out / 'conditions.log', 'after ' + tag, gpu, proc.pid)
        loads.append(after)
        elapsed = time.monotonic() - start
        with log.open('a') as f:
            f.write(f'charged_point_wall={elapsed:.6f}\n')
        row = dict.fromkeys(FIELDS, '')
        row.update(type=typ, kernel=kernel, n=n, method=method, phase=phase,
                   load_before=before, load_after=after, load_max=max(loads), loaded=int(max(loads) >= 16),
                   gpu_busy=int(busy or active), elapsed=elapsed, exit=proc.returncode)
        lines = [line.split(',') for line in text.splitlines() if len(line.split(',')) == 8]
        if proc.returncode == 0 and len(lines) == 1:
            _, _, form, _, _, sec, gf, runs = lines[0]
            row.update(form=form, seconds=sec, gflops=gf, runs=runs)
        append_csv(self.csv, row)
        self.rows.append(row)
        self.spent += elapsed
        print(tag, row['seconds'] or 'FAILED', f'budget={self.spent/3600:.3f}h', flush=True)
        if timed_out:
            raise RuntimeError('measurement deadline; process killed, reaped and charged to CSV')
        return row

    def extend_mpfr(self):
        """After baseline sweeps, replace conservative MPFR caps by actual10s stops.

        Reserve twice the next point's size-scaled forecast (at least60s),
        never reset the shared cap. An inverse-only phase is also available.
        """
        inverse_only = getattr(self.a,'phase','extend-inverses') == 'extend-inverses'
        for typ, k, method in series(self.a.types.split(','), self.a.methods.split(',')):
            if typ != 'mpfr512' or (inverse_only and k not in NAMES[17:]) or (self.a.kernels and k not in self.a.kernels.split(',')):
                continue
            prior = [r for r in self.rows if (r['type'],r['kernel'],r['method'],r['phase']) == (typ,k,method,'sweep')]
            if not prior or any(not r['seconds'] for r in prior):
                continue
            last = max(prior, key=lambda r: int(r['n']))
            vector = k in NAMES[:11]
            power = 1 if vector else 2 if k in NAMES[11:17] else 3
            for n in grid(k, self.a.dense):
                if n <= int(last['n']):
                    continue
                if float(last['seconds']) > 10:
                    break
                runs = 4 if n <= ((1<<21) if vector else 2000) else (2 if method == 'gpu' else 1)
                reserve = max(60., 2*(float(last['seconds'])*(n/int(last['n']))**power*runs + float(last['elapsed'])))
                with (self.out/'reductions.log').open('a') as f:
                    f.write(f'actual MPFR extension {typ} {k} {method}: previous n={last["n"]} seconds={last["seconds"]}; next={n}; reserve={reserve:.3f}s; remaining={self.limit-self.spent:.3f}s\n')
                    if self.limit-self.spent < reserve:
                        f.write('extension stopped: insufficient budget reserve\n')
                        break
                last = self.measure(typ,k,method,n,'sweep')
                if not last['seconds'] or float(last['seconds']) > 10:
                    with (self.out/'reductions.log').open('a') as f:
                        f.write(f'extension actual10s/error stop {typ} {k} {method} n={n}: seconds={last["seconds"]} exit={last["exit"]}\n')
                    break

    def run(self):
        if self.a.phase in ('extend-inverses','extend-mpfr'):
            self.extend_mpfr()
            return
        jobs = list(series(self.a.types.split(','), self.a.methods.split(',')))
        if self.a.kernels:
            jobs = [j for j in jobs if j[1] in self.a.kernels.split(',')]
        for typ, k, m in jobs:
            for n in ([1024, 32768] if k in NAMES[:11] else [64, 128]):
                self.measure(typ, k, m, n, 'pilot')
        if self.a.phase == 'pilot':
            return
        plan = []
        for typ, k, m in jobs:
            pilots = [r for r in self.rows if (r['type'], r['kernel'], r['method'], r['phase']) == (typ, k, m, 'pilot')]
            if len(pilots) != 2 or any(not r['seconds'] for r in pilots):
                plan.append(dict(type=typ, kernel=k, method=m, sizes=[], reason='pilot failed; not validated'))
                continue
            pilot = max(pilots, key=lambda r: int(r['n']))
            sizes = grid(k, self.a.dense)
            power = 1 if k in NAMES[:11] else 2 if k in NAMES[11:17] else 3
            # Conservative forecast, rechecked against every actual observation.
            cap = 10 if typ == 'mpfr512' else 30
            kept = [n for n in sizes if float(pilot['seconds']) * (n/int(pilot['n']))**power <= cap]
            plan.append(dict(type=typ, kernel=k, method=m, sizes=kept, omitted=[n for n in sizes if n not in kept],
                             reason=f'pilot forecast cap {cap}s; matrix grid reduced to powers and nearby multiples7/31'))
        def cost(job, n):
            pp = [r for r in self.rows if (r['type'],r['kernel'],r['method'],r['phase']) ==
                  (job['type'],job['kernel'],job['method'],'pilot')]
            p = max(pp, key=lambda r: int(r['n']))
            vector = job['kernel'] in NAMES[:11]
            power = 1 if vector else 2 if job['kernel'] in NAMES[11:17] else 3
            ratio = n / int(p['n'])
            computation = float(p['seconds']) * ratio**power
            small = min(pp, key=lambda r: int(r['n']))
            setup = max(0, float(p['elapsed']) - 2*float(p['seconds']))
            setup_small = max(0, float(small['elapsed']) - 2*float(small['seconds']))
            memory_power = 1 if vector else 2
            slope = max(0, setup-setup_small) / (int(p['n'])**memory_power-int(small['n'])**memory_power)
            overhead = setup + slope*max(0,n**memory_power-int(p['n'])**memory_power)
            return computation * (4 if n <= ((1<<21) if vector else 2000) else (2 if job['method']=='gpu' else 1)) + overhead
        done = {(r['type'],r['kernel'],r['method'],int(r['n'])) for r in self.rows if r['phase']=='sweep'}
        def pending(j,n):
            return (j['type'],j['kernel'],j['method'],n) not in done
        def total():
            return sum(cost(j,n) for j in plan for n in j['sizes'] if pending(j,n))
        estimate_before = total()
        available = max(0, self.limit-self.spent)*0.8
        if total() > available:
            for j in plan:
                old = j['sizes']
                if not old:
                    continue
                j['sizes'] = [n for n in old if n & (n-1) == 0 or n in (old[0],old[-1],31,63)]
                j['omitted'] += [n for n in old if n not in j['sizes']]
                j['reason'] += '; global budget: remove quarter/interior multiples first'
        while total() > available:
            candidates = [j for j in plan if len(j['sizes']) > 2]
            if not candidates:
                break
            j = max(candidates, key=lambda j: cost(j,j['sizes'][-1]))
            j['omitted'].append(j['sizes'].pop())
            j['reason'] += '; global budget: reduce upper size'
        for j in plan:
            j['estimate_seconds'] = sum(cost(j,n) for n in j['sizes'] if pending(j,n))
        (self.out / ('forecast-' + self.a.methods.replace(',', '-') + '.json')).write_text(json.dumps(
            dict(before_reduction_seconds=estimate_before, after_reduction_seconds=total(),
                 remaining_seconds=self.limit-self.spent, cumulative_spent_seconds=self.spent), indent=2))
        (self.out / ('plan-' + self.a.methods.replace(',', '-') + '.json')).write_text(json.dumps(plan, indent=2))
        for job in plan:
            for n in job['sizes']:
                row = self.measure(job['type'], job['kernel'], job['method'], n, 'sweep')
                if not row['seconds'] or (job['type'] == 'mpfr512' and float(row['seconds']) > 10):
                    with (self.out / 'reductions.log').open('a') as f:
                        f.write(f"{job['type']} {job['kernel']} {job['method']} stop after n={n}: seconds={row['seconds']} exit={row['exit']}\n")
                    break


def plot(a):
    from plot_sweep import draw
    rows = read_rows(Path(a.out) / 'kernels.csv')
    rows = [r for r in rows if r['seconds'] and r['phase'] == 'sweep']
    keys = sorted({(r['type'], r['kernel']) for r in rows if r['type'] in a.types.split(',')
                   and (not a.kernels or r['kernel'] in a.kernels.split(','))})
    for typ, k in keys:
        rr = [r for r in rows if (r['type'], r['kernel']) == (typ, k)]
        # Copies/max index: no invented FLOPs; use element throughput instead.
        def rate(r):
            return int(r['n']) / float(r['seconds']) * 1e-9 if k in NAMES[1:4] else float(r['gflops'])
        gpu = sorted((int(r['n']), rate(r)) for r in rr if r['method'] == 'gpu')
        extra = []
        for method in ('cpu', 'cpu-default'):
            pts = sorted((int(r['n']), rate(r)) for r in rr if r['method'] == method)
            if pts:
                label = 'CPU default' if method == 'cpu-default' or k not in SUMS else 'CPU fast'
                extra.append((label + ' (32 threads)', pts))
        dirty = [(int(r['n']), rate(r)) for r in rr if r['loaded'] == '1' or r['gpu_busy'] == '1']
        d = Path(a.out) / typ
        d.mkdir(exist_ok=True)
        elapsed = sum(float(r['elapsed']) for r in rr)
        draw(gpu, str(d / (k + '.png')), title=f'{typ} {k} — NVIDIA A100 80GB PCIe',
             subtitle=f'Resident operands; CUDA events; no packing/transfers. {len(gpu)} GPU points; all-curve point wall {elapsed/60:.1f} min\n'
                      'GPU warmup1; small best3 / large1 (>10s early stop). CPU OMP_NUM_THREADS=32; small<=2000 / vectors2^21.\n'
                      f'×: load>=16 or other GPU app ({len(dirty)} points); reductions/omissions in README and plan JSON.',
             subtitle_lines=3, xlabel='n (vector length)' if k in NAMES[:11] else
             ('n (square; band k=64)' if k in ('Rtbsv','Ctbsv') else 'n (square matrix)'),
             ylabel='G elements/s (no FLOP count)' if k in NAMES[1:4] else 'GFLOPS (conventional algorithm count)',
             series='GPU fast (resident)' if k in SUMS else 'GPU (resident)',
             extra_cpu=extra, contaminated=dirty, logx=k in NAMES[:11], annotate_min=1, compact_annotations=True)
    print(len(keys), 'graphs')


def bundle(a):
    """Consolidate owned point logs, retaining every command/output/exit receipt.

    Safe after sweeps (not concurrently): CSV is authoritative; unrelated files
    are never touched. Repeating this operation appends only new point files.
    """
    out = Path(a.out)
    rows = read_rows(out / 'kernels.csv')
    for r in rows:
        tag = f"{r['type']}_{r['kernel']}_{r['method']}_{r['n']}_{r['phase']}"
        src = out / 'raw' / (tag + '.log')
        if not src.exists():
            continue
        receipt = src.with_suffix('.log.exit')
        if receipt.read_text().strip() != str(r['exit']):
            raise RuntimeError('exit receipt disagrees with CSV: ' + tag)
        target = out / 'raw' / f"{r['type']}_{r['kernel']}_{r['method']}.log"
        with target.open('a') as f:
            f.write('\n### ' + tag + '\n' + src.read_text() + '\nexit receipt: ' + receipt.read_text())
        src.unlink()
        receipt.unlink()
    for directory in sorted(out.glob('checks*')):
        if not directory.is_dir():
            continue
        for src in sorted(directory.glob('*.log')):
            receipt = src.with_suffix('.log.exit')
            if src.name == 'runs.log' or not receipt.exists():
                continue
            with (directory/'runs.log').open('a') as f:
                f.write('\n### '+src.name+'\n'+src.read_text()+'\nexit receipt: '+receipt.read_text())
            src.unlink()
            receipt.unlink()
    print('point/check logs consolidated; commands/output/exits retained')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', required=True)
    p.add_argument('--bindir', default=str(ROOT / 'benchmark'))
    p.add_argument('--mpfr-bindir', default=str(ROOT / 'benchmark'))
    p.add_argument('--types', default='dd,td,qd,mpfr512')
    p.add_argument('--methods', default='cpu,gpu')
    p.add_argument('--kernels', help='optional comma-separated routine filter')
    p.add_argument('--phase', choices=['pilot', 'sweep', 'extend-inverses', 'extend-mpfr', 'plot', 'bundle'], default='sweep')
    p.add_argument('--hours', type=float, default=5.8)
    p.add_argument('--dense', action='store_true', help='all matrix multiples7/31; usually exceeds budget')
    a = p.parse_args()
    if not 0 < a.hours <= 6:
        p.error('budget must be in (0,6] hours')
    if any(t not in ('dd','td','qd','mpfr512') for t in a.types.split(',')) or any(m not in ('cpu','gpu') for m in a.methods.split(',')):
        p.error('unknown type/method')
    if a.kernels and any(k not in NAMES for k in a.kernels.split(',')):
        p.error('unknown routine in --kernels')
    if a.phase == 'plot':
        plot(a)
    elif a.phase == 'bundle':
        bundle(a)
    else:
        Runner(a).run()

if __name__ == '__main__':
    main()
