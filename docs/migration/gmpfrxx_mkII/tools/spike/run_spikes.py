#!/usr/bin/env python3
"""Compile and run the P0 wrapper and one-way interoperability spikes."""

from __future__ import annotations

import argparse
import json
import os
import shlex
import subprocess
import sys
from pathlib import Path


CASES = (
    ("default-precision-512", "default_precision.cpp", None, True),
    ("same-family-arithmetic", "same_family.cpp", "P2C", False),
    ("utility-print-hex", "utility_print_hex.cpp", "P2C", False),
    ("mplapack-mpfr-surface", "mplapack_mpfr_surface.cpp", "P2C/P3", False),
    ("mplapack-gmp-surface", "mplapack_gmp_surface.cpp", "P2C/P4", False),
    ("legacy-name-collision", "name_collision.cpp", "P2C", False),
    ("expression-lifetime", "expression_lifetime.cpp", "P2C", False),
    ("REAL-DOUBLE", "real_double.cpp", "P2A", False),
    ("REAL-DD", "real_dd.cpp", "P2A", False),
    ("REAL-QD", "real_qd.cpp", "P2A", False),
    ("REAL-BINARY80", "real_binary80.cpp", "P2A", False),
    ("REAL-BINARY128", "real_binary128.cpp", "P2A", False),
    ("REAL-GMP", "real_gmp.cpp", "P2A", False),
    ("COMPLEX-DOUBLE", "complex_double.cpp", "P2B", False),
    ("COMPLEX-DD", "complex_dd.cpp", "P2B", False),
    ("COMPLEX-QD", "complex_qd.cpp", "P2B", False),
    ("COMPLEX-BINARY80", "complex_binary80.cpp", "P2B", False),
    ("COMPLEX-BINARY128", "complex_binary128.cpp", "P2B", False),
    ("COMPLEX-GMP", "complex_gmp.cpp", "P2B", False),
)


def execute(command: list[str], env: dict[str, str]) -> subprocess.CompletedProcess:
    return subprocess.run(
        command,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mplapack-prefix", required=True, type=Path)
    parser.add_argument("--mplapack-source", required=True, type=Path)
    parser.add_argument("--gmpfrxx-prefix", required=True, type=Path)
    parser.add_argument("--work-dir", required=True, type=Path)
    parser.add_argument("--results", required=True, type=Path)
    args = parser.parse_args()

    source_dir = Path(__file__).resolve().parent
    work_dir = args.work_dir.resolve()
    work_dir.mkdir(parents=True, exist_ok=True)
    log_dir = args.results.resolve().parent / "spike-logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    mplapack_prefix = args.mplapack_prefix.resolve()
    gmpfrxx_prefix = args.gmpfrxx_prefix.resolve()
    mplapack_source = args.mplapack_source.resolve()
    env = os.environ.copy()
    prior = env.get("LD_LIBRARY_PATH", "")
    env["LD_LIBRARY_PATH"] = (
        str(mplapack_prefix / "lib") + (f":{prior}" if prior else "")
    )

    common = [
        "c++",
        "-std=gnu++17",
        "-fopenmp",
        "-ffp-contract=off",
        f"-I{gmpfrxx_prefix / 'include'}",
        "-D___MPLAPACK_INTERNAL___",
        f"-I{mplapack_prefix / 'include' / 'mplapack'}",
        f"-I{mplapack_prefix / 'include'}",
        f"-L{mplapack_prefix / 'lib'}",
        "-Wl,-rpath," + str(mplapack_prefix / "lib"),
        f"-I{mplapack_source / 'include'}",
        "-lmpc",
        "-lmpfr",
        "-lgmp",
        "-lqd",
        "-lquadmath",
    ]
    results: list[dict[str, object]] = []
    hard_failures: list[str] = []
    for case_id, filename, owner, must_pass in CASES:
        source = source_dir / filename
        executable = work_dir / case_id
        command = common[:11] + [str(source), "-o", str(executable)] + common[11:]
        compiled = execute(command, env)
        log = [
            f"$ {shlex.join(command)}",
            compiled.stdout.rstrip(),
            f"compile_exit={compiled.returncode}",
        ]
        outcome = "COMPILE_FAIL"
        run_exit: int | None = None
        if compiled.returncode == 0:
            executed = execute([str(executable)], env)
            run_exit = executed.returncode
            log.extend(
                (
                    f"$ {executable}",
                    executed.stdout.rstrip(),
                    f"run_exit={executed.returncode}",
                )
            )
            outcome = "PASS" if executed.returncode == 0 else "RUN_FAIL"
        log_path = log_dir / f"{case_id}.log"
        log_path.write_text("\n".join(log) + "\n", encoding="utf-8")
        result = {
            "id": case_id,
            "source": filename,
            "outcome": outcome,
            "compile_exit": compiled.returncode,
            "run_exit": run_exit,
            "later_owner": owner,
            "must_pass_p0": must_pass,
            "log": str(log_path.relative_to(args.results.resolve().parent)),
        }
        results.append(result)
        print(f"{case_id}: {outcome}")
        if must_pass and outcome != "PASS":
            hard_failures.append(case_id)
        if outcome != "PASS" and not owner:
            hard_failures.append(f"{case_id} (unowned)")

    payload = {
        "schema": 1,
        "default_precision_bits": 512,
        "semantic_macros_enabled": [],
        "results": results,
    }
    args.results.parent.mkdir(parents=True, exist_ok=True)
    args.results.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    if hard_failures:
        print(f"spike gate failures: {hard_failures}", file=sys.stderr)
        return 1
    print("spike execution: PASS (all failures have later owners)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
