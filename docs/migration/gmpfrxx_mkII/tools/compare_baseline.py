#!/usr/bin/env python3
"""Compare a candidate numerical result set with the immutable P0 baseline."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from re import fullmatch
from pathlib import Path


MODES = {
    "exact",
    "status-equal",
    "upper-bound",
    "lower-bound",
    "nonincreasing-error",
}
REQUIRED_FIELDS = {
    "test_id",
    "metric",
    "unit",
    "direction",
    "value",
    "status",
    "source_line",
}


class ComparisonError(RuntimeError):
    pass


def load_results(path: Path) -> dict[tuple[str, str], dict[str, str]]:
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ComparisonError(f"cannot load {path}: {exc}") from exc
    if not isinstance(data, dict) or not isinstance(data.get("results"), list):
        raise ComparisonError(f"{path}: expected object with a results array")

    indexed: dict[tuple[str, str], dict[str, str]] = {}
    for position, item in enumerate(data["results"]):
        if not isinstance(item, dict) or set(item) != REQUIRED_FIELDS:
            raise ComparisonError(
                f"{path}: result {position} must have exactly "
                f"{sorted(REQUIRED_FIELDS)}"
            )
        if not all(isinstance(value, str) and value for value in item.values()):
            raise ComparisonError(f"{path}: result {position} has an empty field")
        key = (item["test_id"], item["metric"])
        if key in indexed:
            raise ComparisonError(f"{path}: duplicate result {key}")
        indexed[key] = item
    if not indexed:
        raise ComparisonError(f"{path}: empty results array")
    return indexed


def load_rules(path: Path) -> dict[str, dict[str, str]]:
    try:
        with path.open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
    except (OSError, UnicodeError, csv.Error) as exc:
        raise ComparisonError(f"cannot load rules {path}: {exc}") from exc
    expected = {"metric", "unit", "direction", "mode"}
    if not rows or set(rows[0]) != expected:
        raise ComparisonError(f"{path}: expected TSV columns {sorted(expected)}")
    rules: dict[str, dict[str, str]] = {}
    for row in rows:
        metric = row["metric"]
        if not metric or metric in rules:
            raise ComparisonError(f"{path}: empty or duplicate metric rule {metric!r}")
        if row["mode"] not in MODES:
            raise ComparisonError(f"{path}: unsupported mode {row['mode']!r}")
        rules[metric] = row
    return rules


def numeric_parts(value: str, key: tuple[str, str]) -> tuple[int, int, str]:
    match = fullmatch(
        r"([+-]?)(?:(\d+)(?:\.(\d*))?|\.(\d+))(?:[eE]([+-]?\d+))?",
        value,
    )
    if match is None:
        raise ComparisonError(f"{key}: expected numeric value, got {value!r}")
    sign_text, integer, fraction, fraction_only, exponent_text = match.groups()
    if fraction_only is not None:
        integer, fraction = "", fraction_only
    integer = integer or ""
    fraction = fraction or ""
    exponent = int(exponent_text or "0")
    combined = integer + fraction
    leading = len(combined) - len(combined.lstrip("0"))
    significant = combined.lstrip("0")
    if not significant:
        return 0, 0, "0"
    sign = -1 if sign_text == "-" else 1
    adjusted_exponent = exponent + len(integer) - leading - 1
    return sign, adjusted_exponent, significant.rstrip("0")


def numeric_compare(left: str, right: str, key: tuple[str, str]) -> int:
    left_sign, left_exponent, left_digits = numeric_parts(left, key)
    right_sign, right_exponent, right_digits = numeric_parts(right, key)
    if left_sign != right_sign:
        return (left_sign > right_sign) - (left_sign < right_sign)
    if left_sign == 0:
        return 0
    if left_exponent != right_exponent:
        absolute = (left_exponent > right_exponent) - (left_exponent < right_exponent)
        return absolute if left_sign > 0 else -absolute
    width = max(len(left_digits), len(right_digits))
    left_padded = left_digits.ljust(width, "0")
    right_padded = right_digits.ljust(width, "0")
    absolute = (left_padded > right_padded) - (left_padded < right_padded)
    return absolute if left_sign > 0 else -absolute


def compare_value(
    key: tuple[str, str], baseline: str, candidate: str, mode: str
) -> bool:
    if mode == "status-equal":
        return candidate == baseline
    comparison = numeric_compare(candidate, baseline, key)
    if mode == "exact":
        return comparison == 0
    if mode in {"upper-bound", "nonincreasing-error"}:
        return comparison <= 0
    if mode == "lower-bound":
        return comparison >= 0
    raise AssertionError(mode)


def compare(baseline_path: Path, candidate_path: Path, rules_path: Path) -> None:
    baseline = load_results(baseline_path)
    candidate = load_results(candidate_path)
    rules = load_rules(rules_path)

    missing = sorted(set(baseline) - set(candidate))
    unknown = sorted(set(candidate) - set(baseline))
    if missing:
        raise ComparisonError(f"missing candidate results: {missing}")
    if unknown:
        raise ComparisonError(f"unknown candidate results: {unknown}")

    failures: list[str] = []
    for key in sorted(baseline):
        expected = baseline[key]
        actual = candidate[key]
        metric = key[1]
        rule = rules.get(metric)
        if rule is None:
            failures.append(f"{key}: no metric-specific rule")
            continue
        for field in ("unit", "direction"):
            if expected[field] != rule[field]:
                failures.append(
                    f"{key}: baseline {field}={expected[field]!r}, "
                    f"rule requires {rule[field]!r}"
                )
            if actual[field] != expected[field]:
                failures.append(
                    f"{key}: candidate changed {field} from "
                    f"{expected[field]!r} to {actual[field]!r}"
                )
        if actual["status"] != expected["status"]:
            failures.append(
                f"{key}: status changed from {expected['status']} "
                f"to {actual['status']}"
            )
        if not compare_value(key, expected["value"], actual["value"], rule["mode"]):
            failures.append(
                f"{key}: {actual['value']} violates {rule['mode']} "
                f"against {expected['value']}"
            )
    if failures:
        raise ComparisonError("\n".join(failures))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("baseline", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--rules", required=True, type=Path)
    args = parser.parse_args()
    try:
        compare(args.baseline, args.candidate, args.rules)
    except ComparisonError as exc:
        print(f"baseline comparison: FAIL\n{exc}", file=sys.stderr)
        return 1
    print("baseline comparison: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
