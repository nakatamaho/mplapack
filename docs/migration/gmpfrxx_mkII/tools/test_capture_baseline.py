#!/usr/bin/env python3

import importlib.util
import sys
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).with_name("capture_baseline.py")
SPEC = importlib.util.spec_from_file_location("p0_capture_baseline", MODULE_PATH)
assert SPEC and SPEC.loader
capture = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = capture
SPEC.loader.exec_module(capture)


class CaptureBaselineTest(unittest.TestCase):
    def test_standard_threshold_summary(self):
        status = capture.output_status([
            "All tests for RGE routines passed the threshold ( 42 tests run)",
        ])
        self.assertEqual(status[:3], ("PASS", 0, 42))

    def test_largest_error_example_summary(self):
        lines = [
            "value of largest test error = +6.25e-03",
            "number of examples where info is not 0 = 0",
            "total number of examples tested = 7",
            "End of tests",
        ]
        self.assertEqual(capture.output_status(lines)[:3], ("PASS", 0, 7))
        match = capture.LARGEST_ERROR_RE.search(lines[0])
        self.assertIsNotNone(match)
        self.assertEqual(capture.decimal_value(match.group(1)), capture.Decimal("0.00625"))

    def test_largest_error_info_failure(self):
        status = capture.output_status([
            "ratio of largest test error = +1.0e+02",
            "number of examples where info is not 0 = 2",
            "total number of examples tested = 10",
            "End of tests",
        ])
        self.assertEqual(status[:3], ("FAIL", 2, 10))

    def test_dmd_summary(self):
        status = capture.output_status([
            "Z - U*V test PASSED.",
            "A*U test PASSED.",
            "Rgedmd :: ALL TESTS PASSED.",
            "Test completed.",
        ])
        self.assertEqual(status[:3], ("PASS", 0, 2))

    def test_arbitrary_large_exponent(self):
        smaller = capture.Decimal("9.99e+1388255822130839109")
        larger = capture.Decimal("1.425e+1388255822130839110")
        self.assertLess(smaller, larger)
        self.assertEqual(str(larger), "1.425e+1388255822130839110")

    def test_incomplete_summary_is_rejected(self):
        with self.assertRaises(capture.BaselineError):
            capture.output_status(["partial output"])
        with self.assertRaises(capture.BaselineError):
            capture.output_status([
                "number of examples where info is not 0 = 0",
                "End of tests",
            ])


if __name__ == "__main__":
    unittest.main()
