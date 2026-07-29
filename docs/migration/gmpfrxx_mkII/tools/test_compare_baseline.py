#!/usr/bin/env python3

from __future__ import annotations

import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location(
    "compare_baseline", HERE / "compare_baseline.py"
)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def entry(test_id: str = "mpblas/mpfr/Rgemm", value: str = "1e-100") -> dict:
    return {
        "test_id": test_id,
        "metric": "max-error",
        "unit": "absolute-error",
        "direction": "lower-is-better",
        "value": value,
        "status": "PASS",
        "source_line": "raw.log:1:maxerror",
    }


class CompareBaselineTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.rules = self.root / "rules.tsv"
        self.rules.write_text(
            "metric\tunit\tdirection\tmode\n"
            "max-error\tabsolute-error\tlower-is-better\tnonincreasing-error\n",
            encoding="utf-8",
        )

    def tearDown(self) -> None:
        self.temp.cleanup()

    def write(self, name: str, results: list[dict]) -> Path:
        path = self.root / name
        path.write_text(json.dumps({"results": results}), encoding="utf-8")
        return path

    def assert_fails(self, baseline: Path, candidate: Path) -> None:
        with self.assertRaises(MODULE.ComparisonError):
            MODULE.compare(baseline, candidate, self.rules)

    def test_identical(self) -> None:
        baseline = self.write("baseline.json", [entry()])
        candidate = self.write("candidate.json", [entry()])
        MODULE.compare(baseline, candidate, self.rules)

    def test_allowed_directional_improvement(self) -> None:
        baseline = self.write("baseline.json", [entry(value="1e-100")])
        candidate = self.write("candidate.json", [entry(value="1e-101")])
        MODULE.compare(baseline, candidate, self.rules)

    def test_worsening(self) -> None:
        baseline = self.write("baseline.json", [entry(value="1e-100")])
        candidate = self.write("candidate.json", [entry(value="1e-99")])
        self.assert_fails(baseline, candidate)

    def test_missing(self) -> None:
        baseline = self.write("baseline.json", [entry()])
        candidate = self.write("candidate.json", [])
        self.assert_fails(baseline, candidate)

    def test_duplicate(self) -> None:
        baseline = self.write("baseline.json", [entry()])
        candidate = self.write("candidate.json", [entry(), entry()])
        self.assert_fails(baseline, candidate)

    def test_unknown(self) -> None:
        baseline = self.write("baseline.json", [entry()])
        candidate = self.write(
            "candidate.json", [entry(), entry(test_id="mpblas/mpfr/unknown")]
        )
        self.assert_fails(baseline, candidate)

    def test_arbitrary_large_exponent(self) -> None:
        huge = "1.425e+1388255822130839110"
        baseline = self.write("baseline.json", [entry(value=huge)])
        equal = self.write("equal.json", [entry(value=huge)])
        improved = self.write(
            "improved.json", [entry(value="9.99e+1388255822130839109")]
        )
        worse = self.write(
            "worse.json", [entry(value="1.426e+1388255822130839110")]
        )
        MODULE.compare(baseline, equal, self.rules)
        MODULE.compare(baseline, improved, self.rules)
        self.assert_fails(baseline, worse)

    def test_numeric_significand_ordering(self) -> None:
        key = ("test", "max-error")
        self.assertGreater(MODULE.numeric_compare("1.2", "1.19", key), 0)
        self.assertEqual(MODULE.numeric_compare("1.20", "1.2", key), 0)
        self.assertLess(MODULE.numeric_compare("-1.2", "-1.19", key), 0)

    def test_malformed_input(self) -> None:
        baseline = self.write("baseline.json", [entry()])
        candidate = self.root / "candidate.json"
        candidate.write_text("{not-json", encoding="utf-8")
        self.assert_fails(baseline, candidate)


if __name__ == "__main__":
    unittest.main()
