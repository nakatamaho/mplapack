#!/usr/bin/env python3

import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).with_name("inventory.py")
SPEC = importlib.util.spec_from_file_location("p0_inventory", MODULE_PATH)
assert SPEC and SPEC.loader
inventory = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = inventory
SPEC.loader.exec_module(inventory)


class InventoryFixtureTest(unittest.TestCase):
    def test_generated_build_and_history_categories(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fixtures = {
                Path("mplapack/reference/Rfoo.cpp"): "mpreal x;\n",
                Path("cmake/MplapackPackage.cmake"): "find_package(GMP)\n",
                Path("docs/migration/gmpfrxx_mkII/INTEROP_PROTOTYPE.md"):
                    "cast2mpf_class(x)\n",
            }
            for relative, content in fixtures.items():
                path = root / relative
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(content)
            hits = inventory.scan(root, fixtures)
            categories = {(hit.path, hit.category) for hit in hits}
            self.assertIn(
                ("mplapack/reference/Rfoo.cpp", "generated"),
                categories,
            )
            self.assertIn(
                ("cmake/MplapackPackage.cmake", "build-metadata"),
                categories,
            )
            self.assertIn(
                (
                    "docs/migration/gmpfrxx_mkII/INTEROP_PROTOTYPE.md",
                    "history-whitelist",
                ),
                categories,
            )


if __name__ == "__main__":
    unittest.main()
