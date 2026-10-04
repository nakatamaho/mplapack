#!/usr/bin/env python3
"""Test Fable comparison policy and DMD-safe source-list regeneration."""
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parent.parent
COMPARATOR = ROOT / "release/docker/common/compare-fable-outputs.py"


class FableComparisonTest(unittest.TestCase):
    def compare(self, *, added=False, removed=False, changed=False):
        with tempfile.TemporaryDirectory(prefix="mplapack-fable-policy-") as tmp:
            baseline, generated = (Path(tmp) / name for name in ("baseline", "generated"))
            for root in (baseline, generated):
                directory = root / "mpblas/reference"
                directory.mkdir(parents=True)
                (directory / "Keep.cpp").write_text("unchanged\n")
                (directory / "Check.cpp").write_text("baseline\n")
            if added:
                (generated / "mpblas/reference/New.cpp").write_text("new\n")
            if removed:
                (generated / "mpblas/reference/Check.cpp").unlink()
            if changed:
                (generated / "mpblas/reference/Check.cpp").write_text("changed\n")
            return subprocess.run(
                [sys.executable, str(COMPARATOR), str(baseline), str(generated)],
                capture_output=True, text=True,
            )

    def test_identical_passes(self):
        result = self.compare()
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_new_file_is_allowed_and_reported(self):
        result = self.compare(added=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("INFO: added mpblas/reference/New.cpp (allowed)", result.stdout)

    def test_removed_file_fails(self):
        result = self.compare(removed=True, added=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("FAIL: removed", result.stdout)

    def test_changed_file_fails(self):
        result = self.compare(changed=True, added=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("FAIL: changed", result.stdout)

    def test_regeneration_excludes_eig_test_imlaenv(self):
        # Exercise the generator's actual embedded source-list program.
        script = (ROOT / "fable/go_testing.sh").read_text()
        function = script.split("regen_testing_makefile_sources() {", 1)[1]
        program = function.split("<<'PY'\n", 1)[1].split("\nPY", 1)[0]
        with tempfile.TemporaryDirectory(prefix="mplapack-fable-sources-") as tmp:
            directory = Path(tmp)
            for name in ("Alahdg.cpp", "iMlaenv.cpp"):
                (directory / name).write_text("// fixture\n")
            for variable in ("EIG_SOURCES", "LIN_SOURCES"):
                makefile = directory / "Makefile.am"
                makefile.write_text(f"{variable} = old.cpp\n")
                subprocess.run(
                    [sys.executable, "-", str(makefile), str(directory), variable, "common/"],
                    input=program, text=True, check=True,
                )
                text = makefile.read_text()
                self.assertIn("common/Alahdg.cpp", text)
                self.assertEqual("common/iMlaenv.cpp" in text, variable == "LIN_SOURCES")


if __name__ == "__main__":
    unittest.main()
