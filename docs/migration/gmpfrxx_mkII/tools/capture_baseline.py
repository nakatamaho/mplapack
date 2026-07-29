#!/usr/bin/env python3
"""Normalize the locked MPFR/GMP Autotools test results for P0."""

from __future__ import annotations

import argparse
import gzip
import json
import re
from pathlib import Path


BASE_SHA = "b875e74d4b927282c907c3a29e6cadda48a7d57b"
NUMBER_RE = re.compile(
    r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][+-]?\d+)?"
)
MAX_ERROR_RE = re.compile(r"\bmaxerror:\s*(" + NUMBER_RE.pattern + r")")
TESTS_RUN_RE = re.compile(r"\(\s*(\d+)\s+tests?\s+run\s*\)", re.IGNORECASE)
FAILED_COUNT_RE = re.compile(
    r"\b(\d+)\s+out\s+of\s+(\d+).*?\bfailed\b", re.IGNORECASE
)
LARGEST_ERROR_RE = re.compile(
    r"\b(?:value|ratio) of largest test error\s*=\s*(" + NUMBER_RE.pattern + r")",
    re.IGNORECASE,
)
INFO_FAILURES_RE = re.compile(
    r"\bnumber of examples where info is not 0\s*=\s*(\d+)", re.IGNORECASE
)
TOTAL_EXAMPLES_RE = re.compile(
    r"\btotal number of examples tested\s*=\s*(\d+)", re.IGNORECASE
)
DMD_PASS_RE = re.compile(r"\btest\s+PASSED\b", re.IGNORECASE)
DMD_FAIL_RE = re.compile(r"\btest\s+FAILED\b", re.IGNORECASE)


class BaselineError(RuntimeError):
    pass


def source_line(relative: Path, line_no: int, text: str) -> str:
    return f"{relative.as_posix()}:{line_no}:{text.rstrip()}"


def read_lines(path: Path) -> list[str]:
    try:
        return path.read_text(encoding="utf-8", errors="strict").splitlines()
    except (OSError, UnicodeError) as exc:
        raise BaselineError(f"cannot read {path}: {exc}") from exc


class Decimal:
    """Comparable decimal with an arbitrary-size base-10 exponent."""

    def __init__(self, token: str):
        match = NUMBER_RE.fullmatch(token)
        if match is None:
            raise ValueError(token)
        normalized = token.replace("D", "E").replace("d", "e")
        mantissa, marker, exponent_text = normalized.lower().partition("e")
        exponent = int(exponent_text) if marker else 0
        sign = -1 if mantissa.startswith("-") else 1
        mantissa = mantissa.lstrip("+-")
        integer, dot, fraction = mantissa.partition(".")
        if not dot:
            fraction = ""
        combined = integer + fraction
        leading = len(combined) - len(combined.lstrip("0"))
        significant = combined.lstrip("0")
        if not significant:
            self.sign = 0
            self.adjusted_exponent = 0
            self.digits = "0"
            self.text = "0"
            return
        self.sign = sign
        self.adjusted_exponent = exponent + len(integer) - leading - 1
        self.digits = significant.rstrip("0")
        self.text = normalized.lstrip("+")

    def _absolute_compare(self, other: "Decimal") -> int:
        if self.adjusted_exponent != other.adjusted_exponent:
            return -1 if self.adjusted_exponent < other.adjusted_exponent else 1
        width = max(len(self.digits), len(other.digits))
        left = self.digits.ljust(width, "0")
        right = other.digits.ljust(width, "0")
        return (left > right) - (left < right)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Decimal):
            return NotImplemented
        if self.sign != other.sign:
            return False
        return self.sign == 0 or self._absolute_compare(other) == 0

    def __lt__(self, other: "Decimal") -> bool:
        if self.sign != other.sign:
            return self.sign < other.sign
        if self.sign == 0:
            return False
        comparison = self._absolute_compare(other)
        return comparison < 0 if self.sign > 0 else comparison > 0

    def __abs__(self) -> "Decimal":
        return Decimal(self.text.lstrip("-"))

    def __str__(self) -> str:
        return self.text


def decimal_value(token: str) -> Decimal:
    try:
        return Decimal(token)
    except ValueError as exc:
        raise BaselineError(f"malformed numeric value: {token}") from exc


def result(
    test_id: str,
    metric: str,
    unit: str,
    direction: str,
    value: str,
    status: str,
    line: str,
) -> dict[str, str]:
    return {
        "test_id": test_id,
        "metric": metric,
        "unit": unit,
        "direction": direction,
        "value": value,
        "status": status,
        "source_line": line,
    }


def parse_trs_status(path: Path) -> tuple[str, int, str]:
    lines = read_lines(path)
    for line_no, line in enumerate(lines, 1):
        if line.startswith(":test-result:"):
            status = line.split(":", 2)[2].strip()
            if status not in {"PASS", "FAIL", "SKIP", "XFAIL", "XPASS", "ERROR"}:
                raise BaselineError(f"unknown status in {path}:{line_no}: {status}")
            return status, line_no, line
    raise BaselineError(f"missing :test-result: in {path}")


def parse_mpblas(
    build_root: Path, backend: str
) -> tuple[list[dict[str, str]], list[Path]]:
    directory = build_root / "mpblas" / "test" / backend
    logs = sorted(directory.glob("run_*.log"))
    if not logs:
        raise BaselineError(f"no MPBLAS {backend} wrapper logs in {directory}")

    normalized: list[dict[str, str]] = []
    raw_paths: list[Path] = []
    for log_path in logs:
        trs_path = log_path.with_suffix(".trs")
        if not trs_path.is_file():
            raise BaselineError(f"missing Automake result file: {trs_path}")
        status, trs_line_no, trs_line = parse_trs_status(trs_path)
        log_lines = read_lines(log_path)
        relative_log = log_path.relative_to(build_root)
        relative_trs = trs_path.relative_to(build_root)
        test_name = log_path.name.removeprefix("run_").removesuffix(".log")
        test_id = f"mpblas/{backend}/{test_name}"

        maxima: list[tuple[Decimal, int, str]] = []
        for line_no, line in enumerate(log_lines, 1):
            match = MAX_ERROR_RE.search(line)
            if match:
                maxima.append((abs(decimal_value(match.group(1))), line_no, line))
        if maxima:
            maximum, line_no, line = max(maxima, key=lambda item: item[0])
            normalized.append(
                result(
                    test_id,
                    "max-error",
                    "absolute-error",
                    "lower-is-better",
                    str(maximum),
                    status,
                    source_line(relative_log, line_no, line),
                )
            )

        normalized.append(
            result(
                test_id,
                "exit-status",
                "status",
                "equal",
                status,
                status,
                source_line(relative_trs, trs_line_no, trs_line),
            )
        )
        raw_paths.extend((log_path, trs_path))
    return normalized, raw_paths


def output_status(lines: list[str]) -> tuple[str, int, int, int, str]:
    failures = 0
    executed = 0
    status_line_no = 0
    status_line = ""
    saw_pass = False
    info_failures: int | None = None
    total_examples: int | None = None
    dmd_passes = 0
    dmd_failures = 0
    dmd_completed = False

    for line_no, line in enumerate(lines, 1):
        for match in TESTS_RUN_RE.finditer(line):
            executed += int(match.group(1))
        for match in FAILED_COUNT_RE.finditer(line):
            failures += int(match.group(1))
            executed += int(match.group(2))
            status_line_no, status_line = line_no, line
        info_match = INFO_FAILURES_RE.search(line)
        if info_match:
            info_failures = int(info_match.group(1))
            status_line_no, status_line = line_no, line
        total_match = TOTAL_EXAMPLES_RE.search(line)
        if total_match:
            total_examples = int(total_match.group(1))
        dmd_passes += len(DMD_PASS_RE.findall(line))
        dmd_failures += len(DMD_FAIL_RE.findall(line))
        if "test completed" in line.lower():
            dmd_completed = True
            status_line_no, status_line = line_no, line
        lowered = line.lower()
        if "passed the threshold" in lowered or "passed the tests" in lowered:
            saw_pass = True
            status_line_no, status_line = line_no, line

    if info_failures is not None or total_examples is not None:
        if info_failures is None or total_examples is None:
            raise BaselineError("incomplete largest-error example summary")
        status = "FAIL" if info_failures else "PASS"
        return status, info_failures, total_examples, status_line_no, status_line
    if dmd_passes or dmd_failures or dmd_completed:
        if not dmd_completed or not (dmd_passes or dmd_failures):
            raise BaselineError("incomplete DMD summary")
        status = "FAIL" if dmd_failures else "PASS"
        return status, dmd_failures, dmd_passes + dmd_failures, status_line_no, status_line
    if failures:
        return "FAIL", failures, executed, status_line_no, status_line
    if not saw_pass:
        raise BaselineError("result output has neither pass summary nor failure count")
    return "PASS", 0, executed, status_line_no, status_line


def parse_driver_outputs(
    build_root: Path, family: str, backend: str
) -> tuple[list[dict[str, str]], list[Path]]:
    result_root = build_root / "mplapack" / "test" / family / "results"
    pattern = f"*/{backend}/*.out"
    paths = sorted(result_root.glob(pattern))
    if backend == "mpfr":
        paths = [path for path in paths if path.name.endswith(".default.out")]
    if not paths:
        raise BaselineError(
            f"no {family.upper()} {backend} result outputs matching {result_root / pattern}"
        )

    normalized: list[dict[str, str]] = []
    for path in paths:
        lines = read_lines(path)
        try:
            status, failures, executed, line_no, line = output_status(lines)
        except BaselineError as exc:
            raise BaselineError(f"{path}: {exc}") from exc
        relative = path.relative_to(build_root)
        profile_name = path.name.removesuffix(".out")
        test_id = f"{family}/{backend}/{profile_name}"
        summary = source_line(relative, line_no, line)
        largest_errors = []
        for metric_line_no, metric_line in enumerate(lines, 1):
            match = LARGEST_ERROR_RE.search(metric_line)
            if match:
                largest_errors.append(
                    (abs(decimal_value(match.group(1))), metric_line_no, metric_line)
                )
        if largest_errors:
            maximum, metric_line_no, metric_line = max(
                largest_errors, key=lambda item: item[0]
            )
            normalized.append(
                result(
                    test_id,
                    "test-error-ratio",
                    "ratio",
                    "lower-is-better",
                    str(maximum),
                    status,
                    source_line(relative, metric_line_no, metric_line),
                )
            )
        normalized.extend(
            (
                result(
                    test_id,
                    "failure-count",
                    "tests",
                    "lower-is-better",
                    str(failures),
                    status,
                    summary,
                ),
                result(
                    test_id,
                    "executed-tests",
                    "tests",
                    "equal",
                    str(executed),
                    status,
                    summary,
                ),
                result(
                    test_id,
                    "exit-status",
                    "status",
                    "equal",
                    status,
                    status,
                    summary,
                ),
            )
        )
    return normalized, paths


def write_raw_archive(build_root: Path, paths: list[Path], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("wb") as raw_file:
        with gzip.GzipFile(
            filename="", mode="wb", fileobj=raw_file, compresslevel=9, mtime=0
        ) as compressed:
            for path in sorted(set(paths)):
                relative = path.relative_to(build_root).as_posix()
                compressed.write(f"===== BEGIN {relative} =====\n".encode())
                data = path.read_bytes()
                compressed.write(data)
                if data and not data.endswith(b"\n"):
                    compressed.write(b"\n")
                compressed.write(f"===== END {relative} =====\n".encode())


def capture(build_root: Path) -> tuple[list[dict[str, str]], list[Path]]:
    results: list[dict[str, str]] = []
    raw_paths: list[Path] = []
    for backend in ("gmp", "mpfr"):
        parsed, raw = parse_mpblas(build_root, backend)
        results.extend(parsed)
        raw_paths.extend(raw)
        for family in ("lin", "eig"):
            parsed, raw = parse_driver_outputs(build_root, family, backend)
            results.extend(parsed)
            raw_paths.extend(raw)

    results.sort(key=lambda item: (item["test_id"], item["metric"]))
    keys = [(item["test_id"], item["metric"]) for item in results]
    duplicates = sorted({key for key in keys if keys.count(key) > 1})
    if duplicates:
        raise BaselineError(f"duplicate normalized result keys: {duplicates}")
    return results, raw_paths


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--build-root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--raw-output", required=True, type=Path)
    args = parser.parse_args()

    build_root = args.build_root.resolve()
    results, raw_paths = capture(build_root)
    payload = {
        "schema": 1,
        "mplapack_base_sha": BASE_SHA,
        "default_precision_bits": 512,
        "test_scope": {
            "backends": ["gmp", "mpfr"],
            "mpfr_profiles": ["default"],
            "autotools": "complete discovered MPBLAS, LIN, and EIG sets",
        },
        "raw_source_count": len(set(raw_paths)),
        "results": results,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_raw_archive(build_root, raw_paths, args.raw_output)
    print(
        f"captured {len(results)} normalized results from "
        f"{len(set(raw_paths))} raw files"
    )
    print(f"baseline: {args.output}")
    print(f"raw archive: {args.raw_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
