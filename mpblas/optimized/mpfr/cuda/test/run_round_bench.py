#!/usr/bin/env python3
"""Interleaved, separate-process Rtrsm old/new timing; stdout is raw CSV."""
import csv
import os
from pathlib import Path
import re
import subprocess
import sys

exe = Path(sys.argv[1]) / "mpfr_cuda_Rtrsm_round_bench"
rounds = int(sys.argv[2]) if len(sys.argv) > 2 else 5
cases = [(p, n, flags) for p in (256, 512, 1024)
         for n in (256, 1024) for flags in ("LUN", "RLT")]
cases += [(512, n, flags) for n in (512, 2048) for flags in ("LUN", "RLT")]
if len(sys.argv) > 3:
    if sys.argv[3] != "acceptance":
        raise SystemExit("optional third argument: acceptance")
    cases = [(512, n, flags) for n in (512, 1024) for flags in ("LUN", "RLT")]
w = csv.writer(sys.stdout)
w.writerow(("PB", "n", "flags", "round", "runtime", "seconds"))
for p, n, flags in cases:
    for repeat in range(rounds):
        for runtime in ((0, 1) if repeat % 2 == 0 else (1, 0)):
            env = dict(os.environ, MPLAPACK_MPFR_CUDA_FORCE_RUNTIME=str(runtime),
                       MPLAPACK_MPFR_CUDA_VERBOSE="1", OMP_NUM_THREADS="1")
            result = subprocess.run([str(exe), str(p), str(n), *flags],
                                    env=env, check=True, text=True,
                                    stdout=subprocess.PIPE)
            match = re.search(r"seconds=([0-9.]+) ok=1", result.stdout)
            if not match:
                raise RuntimeError(result.stdout)
            w.writerow((p, n, flags, repeat, runtime, match[1]))
            sys.stdout.flush()
