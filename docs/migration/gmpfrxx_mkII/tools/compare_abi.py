#!/usr/bin/env python3
"""Compare ABI/install manifests with phase-scoped backend allowances."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


BACKENDS = {"gmp", "mpfr", "qd", "dd", "double", "binary80", "binary128"}
ALLOWED_BACKEND = {"P0": None, "P3": "mpfr", "P4": "gmp"}


class AbiComparisonError(RuntimeError):
    pass


def load(path: Path) -> dict:
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise AbiComparisonError(f"cannot load {path}: {exc}") from exc
    if not isinstance(data, dict) or data.get("schema") != 1:
        raise AbiComparisonError(f"{path}: expected ABI manifest schema 1")
    if not isinstance(data.get("backends"), dict):
        raise AbiComparisonError(f"{path}: missing backends object")
    if set(data["backends"]) != BACKENDS:
        raise AbiComparisonError(
            f"{path}: backend set is {sorted(data['backends'])}, "
            f"expected {sorted(BACKENDS)}"
        )
    return data


def compare(baseline_path: Path, candidate_path: Path, phase: str) -> None:
    if phase not in ALLOWED_BACKEND:
        raise AbiComparisonError(f"unsupported phase {phase}")
    baseline = load(baseline_path)
    candidate = load(candidate_path)
    baseline_global = {key: value for key, value in baseline.items() if key != "backends"}
    candidate_global = {
        key: value for key, value in candidate.items() if key != "backends"
    }
    failures: list[str] = []
    if baseline_global != candidate_global:
        failures.append("global ABI/install metadata changed")

    allowed = ALLOWED_BACKEND[phase]
    for backend in sorted(BACKENDS):
        if backend == allowed:
            continue
        if baseline["backends"][backend] != candidate["backends"][backend]:
            failures.append(f"non-migrated backend changed: {backend}")
    if failures:
        raise AbiComparisonError("\n".join(failures))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("baseline", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--phase", choices=sorted(ALLOWED_BACKEND), required=True)
    args = parser.parse_args()
    try:
        compare(args.baseline, args.candidate, args.phase)
    except AbiComparisonError as exc:
        print(f"ABI comparison: FAIL\n{exc}", file=sys.stderr)
        return 1
    print(f"ABI comparison: PASS ({args.phase})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
