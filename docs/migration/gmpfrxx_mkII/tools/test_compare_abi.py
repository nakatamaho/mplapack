#!/usr/bin/env python3

from __future__ import annotations

import copy
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("compare_abi", HERE / "compare_abi.py")
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def manifest() -> dict:
    return {
        "schema": 1,
        "mplapack_base_sha": "a" * 40,
        "platform": {"system": "Linux"},
        "cmake_package_dependencies": [],
        "backends": {
            backend: {"abi": f"locked-{backend}"} for backend in MODULE.BACKENDS
        },
    }


class CompareAbiTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)

    def tearDown(self) -> None:
        self.temp.cleanup()

    def write(self, name: str, value: dict) -> Path:
        path = self.root / name
        path.write_text(json.dumps(value), encoding="utf-8")
        return path

    def test_self_comparison(self) -> None:
        data = manifest()
        path = self.write("same.json", data)
        MODULE.compare(path, path, "P0")

    def test_p3_allows_only_mpfr(self) -> None:
        baseline_data = manifest()
        candidate_data = copy.deepcopy(baseline_data)
        candidate_data["backends"]["mpfr"]["abi"] = "changed"
        baseline = self.write("baseline.json", baseline_data)
        candidate = self.write("candidate.json", candidate_data)
        MODULE.compare(baseline, candidate, "P3")
        with self.assertRaises(MODULE.AbiComparisonError):
            MODULE.compare(baseline, candidate, "P0")

    def test_p4_allows_only_gmp(self) -> None:
        baseline_data = manifest()
        candidate_data = copy.deepcopy(baseline_data)
        candidate_data["backends"]["gmp"]["abi"] = "changed"
        baseline = self.write("baseline.json", baseline_data)
        candidate = self.write("candidate.json", candidate_data)
        MODULE.compare(baseline, candidate, "P4")

    def test_non_migrated_backend_is_immutable(self) -> None:
        baseline_data = manifest()
        candidate_data = copy.deepcopy(baseline_data)
        candidate_data["backends"]["dd"]["abi"] = "changed"
        baseline = self.write("baseline.json", baseline_data)
        candidate = self.write("candidate.json", candidate_data)
        with self.assertRaises(MODULE.AbiComparisonError):
            MODULE.compare(baseline, candidate, "P3")


if __name__ == "__main__":
    unittest.main()
