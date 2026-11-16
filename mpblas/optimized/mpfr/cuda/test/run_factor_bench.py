#!/usr/bin/env python3
"""Interleaved, fresh-process measurements. Append raw records as JSON lines.

Usage: run_factor_bench.py CPU_BINARY OLD_CUDA_BINARY NEW_CUDA_BINARY OUT [ROUNDS] [complex]
Inputs have full-P significands, deterministic seed, and strict diagonal dominance.
CPU/old/entry order reverses on alternate rounds; report medians, never estimates.
"""
import json
import os
from pathlib import Path
import re
import statistics
import subprocess
import sys

cpu, old, new, output = sys.argv[1:5]
rounds = int(sys.argv[5]) if len(sys.argv) > 5 else 3
complex_mode = len(sys.argv) > 6 and sys.argv[6] == "complex"
cases = [(512, n, routine) for n in ((512, 1024, 2048) if complex_mode else (512, 1024, 2048, 4096))
         for routine in ("getrf", "potrf")]
if not complex_mode:
    cases += [(p, 512, r) for p in (256, 1024, 200) for r in ("getrf", "potrf")]
records = []
# Compare default algorithms, not tuning inherited from the invoking shell.
base_env = {k: v for k, v in os.environ.items()
            if not k.startswith(("MPLAPACK_LAPACK_", "MPLAPACK_MPFR_"))
            and k not in ("MPFRXX_DEFAULT_EMIN", "MPFRXX_DEFAULT_EMAX")}
with Path(output).open("a") as log:
    for turn in range(rounds):
        order = [("cpu", cpu), ("old", old), ("entry", new)]
        if turn % 2:
            order.reverse()
        for prec, n, routine in cases:
            for mode, binary in order:
                result = subprocess.run([binary, str(prec), str(n), routine, mode],
                                        env={**base_env, "OMP_NUM_THREADS": "32",
                                             "MPFRXX_DEFAULT_PRECISION_BITS": str(prec),
                                             "MPFRXX_DEFAULT_ROUNDING_MODE": "RNDN"},
                                        capture_output=True, text=True)
                row = {"round": turn, "prec": prec, "n": n, "routine": routine,
                       "mode": mode, "returncode": result.returncode,
                       "stdout": result.stdout, "stderr": result.stderr}
                fields = result.stdout.strip().splitlines()[-1].split(",") if result.stdout.strip() else []
                if len(fields) == 7:
                    row["seconds"] = float(fields[4])
                for key in ("upload", "factor", "download"):
                    match = re.search(rf"{key}=([0-9.eE+-]+)", result.stderr)
                    if match:
                        row[key] = float(match.group(1))
                log.write(json.dumps(row) + "\n")
                log.flush()
                if result.returncode:
                    raise RuntimeError(row)
                records.append(row)
for prec, n, routine in cases:
    line = [prec, n, routine]
    for mode in ("cpu", "old", "entry"):
        line.append(statistics.median(r["seconds"] for r in records
                    if (r["prec"], r["n"], r["routine"], r["mode"]) == (prec, n, routine, mode)))
    print(*line, flush=True)
