#!/usr/bin/env python3
"""Create the immutable P0 source inventory for the gmpfrxx_mkII migration."""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable


MIGRATION_DIR = Path("docs/migration/gmpfrxx_mkII")
OUTPUT_PATHS = {
    MIGRATION_DIR / "INVENTORY.md",
    MIGRATION_DIR / "inventory.json",
    MIGRATION_DIR / "REPORT-P0.md",
}
PROTOTYPE_PATH = MIGRATION_DIR / "INTEROP_PROTOTYPE.md"

PATTERNS: tuple[tuple[str, re.Pattern[str]], ...] = (
    (
        "legacy-mpfr-type-or-header",
        re.compile(r"\bmpreal\b|\bmpcomplex\b|mpreal\.h|mpcomplex\.h|using\s+namespace\s+mpfr"),
    ),
    (
        "legacy-gmp-type-or-header",
        re.compile(r"gmpxx\.h|(?<![\w:])mpf_class\b|(?<![\w:])mpc_class\b|include/mpc_class\.h"),
    ),
    (
        "gmpxx-link-or-build",
        re.compile(r"-lgmpxx\b|libgmpxx\b|--enable-cxx\b|gmpxx\.pc"),
    ),
    (
        "precision-default",
        re.compile(
            r"set_default_prec|set_default_precision|set_default_mpf_precision|"
            r"mpf_set_default_prec|mpfr_set_default_prec|default_prec|"
            r"default_precision_bits|default_mpf_precision_bits|"
            r"\bprecision\(\)|\bget_prec\("
        ),
    ),
    (
        "comparison-interop",
        re.compile(
            r"cast2(?:mpf_class|mpc_class|double|qd_real|qd_complex|dd_real|"
            r"dd_complex|binary80|binary128|complex_binary)|"
            r"\binfnorm(?:_mat)?\b|\breldiff\b|REAL_REF|COMPLEX_REF"
        ),
    ),
    ("reverse-cast2", re.compile(r"\bcast2[A-Za-z0-9_]*\s*\(")),
    (
        "utility-math-print",
        re.compile(
            r"\b(?:Mlsame|Mxerbla|Msign|Mmax|Mmin|Mutils|Rlamch|"
            r"printnum(?:_short)?|printvec|printmat|printnum_hex|"
            r"urandomb|urandom_c|abs|sqrt|exp|log|log2|log10|sin|cos|"
            r"atan2|pow|ceil|floor|nextafter|ldexp)\s*\("
        ),
    ),
    (
        "fable-codegen",
        re.compile(r"\bfable\b|gen_include_mplapack|generate\.sh|FABLE|patch-[A-ZRC]"),
    ),
    (
        "test-executable-or-metric",
        re.compile(
            r"(?:check_PROGRAMS|TESTS|_PROGRAMS)\s*=|"
            r"\b(?:maxdiff|EPSILON(?:[0-9]+)?|infnorm|reldiff|PASS|FAIL|metric)\b"
        ),
    ),
    (
        "abi-install-package-consumer",
        re.compile(
            r"version-info|SOVERSION|SONAME|install_name|pkgconfig|pkg-config|"
            r"add_library|target_link_libraries|INTERFACE_LINK_LIBRARIES|"
            r"install\(TARGETS|configure_package_config_file|find_package|"
            r"mplapackConfig|mplapack\.pc|Makefile\.(?:linux|macos|mingw|freebsd)"
        ),
    ),
)

PRECISION_SETTER = re.compile(
    r"(?:set_default_prec|set_default_precision_bits|"
    r"set_default_mpf_precision_bits|mpf_set_default_prec|mpfr_set_default_prec)"
    r"\s*\(([^)]*)\)"
)
PRECISION_ASSIGNMENT = re.compile(
    r"(?:default_prec|default_real_prec|default_imag_prec)\s*=\s*([^;]+)"
)
PRECISION_DECLARATION = re.compile(
    r"^(?:inline\s+|static\s+)*(?:void|bool|int|mpfr_prec_t|mp_bitcnt_t)\s+"
    r"(?:[A-Za-z_][A-Za-z0-9_]*::)?(?:set_default_prec|"
    r"set_default_precision_bits|set_default_mpf_precision_bits)\s*\("
)
PRECISION_512 = re.compile(
    r"(?:^|\D)512(?:\D|$)|___MPREAL_DEFAULT_PRECISION___|"
    r"MPLAPACK_(?:GMP|MPFR)_DEFAULT_PRECISION"
)
PRECISION_QUERY = re.compile(
    r"(?:get_default_prec|default_precision_bits|default_mpf_precision_bits|"
    r"get_prec|precision)\s*\("
)
PER_OBJECT_PRECISION = re.compile(r"(?:with_precision|mp[fr]*_init2|set_prec)\s*\(")


@dataclass(frozen=True)
class Hit:
    path: str
    line: int
    category: str
    generated: bool
    owning_build_targets: list[str]
    precision_classification: str | None
    text: str


def is_generated(path: Path) -> bool:
    value = path.as_posix()
    if value.startswith(("mpblas/reference/", "mplapack/reference/")):
        return True
    if re.match(r"include/mplapack_(?:eig_|lin_|matgen_)?(?:gmp|mpfr|qd|dd|double|binary80|binary128)\.h$", value):
        return True
    if "/test/" in value and "/common/" not in value and path.suffix in {".cpp", ".h"}:
        return True
    return False


def path_category(path: Path) -> str:
    value = path.as_posix()
    if value == "docs/migration/gmpfrxx_mkII/INTEROP_PROTOTYPE.md":
        return "history-whitelist"
    if path.name in {
        "Makefile.am",
        "Makefile.in",
        "CMakeLists.txt",
        "configure.ac",
    } or value.startswith(("cmake/", "packaging/", "release/", ".github/")):
        return "build-metadata"
    if is_generated(path):
        return "generated"
    return "source"


def owning_targets(path: Path) -> list[str]:
    value = path.as_posix()
    backends = [
        backend
        for backend in ("gmp", "mpfr", "qd", "dd", "double", "binary80", "binary128")
        if re.search(rf"(?:^|[_.\/]){backend}(?:[_.\/]|$)", value)
    ]
    if value.startswith("mpblas/reference/"):
        return [f"libmplapack_{backend}" for backend in backends] or ["all-backend-libraries"]
    if value.startswith("mplapack/reference/"):
        return [f"libmplapack_{backend}" for backend in backends] or ["all-backend-libraries"]
    if "/test/" in value:
        return [f"{backend}-test-suite" for backend in backends] or ["shared-test-support"]
    if value.startswith("include/"):
        return [f"installed-{backend}-headers" for backend in backends] or ["installed-public-headers"]
    if value.startswith(("fable/", "misc/")):
        return ["maintainer-codegen-tools"]
    if value.startswith(("examples/", "benchmark/")):
        root = value.split("/", 1)[0]
        return [f"{root}-{backend}" for backend in backends] or [root]
    if path_category(path) == "build-metadata":
        return ["autotools-and-cmake-build-graph"]
    if path_category(path) == "history-whitelist":
        return ["p0-provenance-only"]
    return ["source-tree"]


def precision_classification(line: str) -> str | None:
    stripped = line.strip()
    if PRECISION_DECLARATION.search(stripped):
        return None
    setter = PRECISION_SETTER.search(line)
    if setter:
        return "sets-512" if PRECISION_512.search(setter.group(1)) else "sets-non-512"
    assignment = PRECISION_ASSIGNMENT.search(line)
    if assignment:
        return "sets-512" if PRECISION_512.search(assignment.group(1)) else "sets-non-512"
    if PER_OBJECT_PRECISION.search(line):
        return "per-object"
    if PRECISION_QUERY.search(line):
        return "query-only"
    return None


def list_files(root: Path) -> list[Path]:
    proc = subprocess.run(
        ["git", "ls-files", "-z", "--cached", "--others", "--exclude-standard"],
        cwd=root,
        check=True,
        stdout=subprocess.PIPE,
    )
    paths = []
    for raw in proc.stdout.split(b"\0"):
        if not raw:
            continue
        path = Path(os.fsdecode(raw))
        if path in OUTPUT_PATHS:
            continue
        if path.is_relative_to(MIGRATION_DIR) and path != PROTOTYPE_PATH:
            continue
        if path.parts[:2] == (".git", "objects"):
            continue
        paths.append(path)
    return sorted(set(paths), key=lambda item: item.as_posix())


def scan_file(root: Path, path: Path) -> list[Hit]:
    absolute = root / path
    data = absolute.read_bytes()
    if b"\0" in data[:8192]:
        return []
    text = data.decode("utf-8", errors="surrogateescape")
    hits: list[Hit] = []
    base_category = path_category(path)
    for line_number, line in enumerate(text.splitlines(), start=1):
        categories = [category for category, pattern in PATTERNS if pattern.search(line)]
        if not categories:
            continue
        if base_category in {"generated", "build-metadata", "history-whitelist"}:
            categories.append(base_category)
        for category in sorted(set(categories)):
            hits.append(
                Hit(
                    path=path.as_posix(),
                    line=line_number,
                    category=category,
                    generated=is_generated(path),
                    owning_build_targets=owning_targets(path),
                    precision_classification=(
                        precision_classification(line)
                        if category == "precision-default"
                        else None
                    ),
                    text=line.strip()[:500],
                )
            )
    return hits


def scan(root: Path, files: Iterable[Path] | None = None) -> list[Hit]:
    selected = list(files) if files is not None else list_files(root)
    hits: list[Hit] = []
    errors: list[str] = []
    for path in selected:
        try:
            hits.extend(scan_file(root, path))
        except (OSError, UnicodeError) as exc:
            errors.append(f"{path.as_posix()}: {exc}")
    if errors:
        raise RuntimeError("unreadable inventory files:\n" + "\n".join(errors))
    return sorted(
        hits,
        key=lambda hit: (hit.path, hit.line, hit.category, hit.text),
    )


def render_markdown(hits: list[Hit], root: Path) -> str:
    counts: dict[str, int] = {}
    for hit in hits:
        counts[hit.category] = counts.get(hit.category, 0) + 1
    lines = [
        "# P0 Source Inventory",
        "",
        f"Repository root: `{root}`",
        "",
        "Generated by `tools/inventory.py`; do not edit by hand.",
        "",
        "## Category counts",
        "",
        "| Category | Hits |",
        "| --- | ---: |",
    ]
    for category, count in sorted(counts.items()):
        lines.append(f"| `{category}` | {count} |")
    lines.extend(
        [
            "",
            "## Hits",
            "",
            "| Path | Line | Category | Generated | Owning build targets | Precision | Text |",
            "| --- | ---: | --- | --- | --- | --- | --- |",
        ]
    )
    for hit in hits:
        escaped = hit.text.replace("|", "\\|").replace("`", "\\`")
        targets = ", ".join(f"`{target}`" for target in hit.owning_build_targets)
        precision = hit.precision_classification or ""
        lines.append(
            f"| `{hit.path}` | {hit.line} | `{hit.category}` | "
            f"{'yes' if hit.generated else 'no'} | {targets} | `{precision}` | "
            f"`{escaped}` |"
        )
    lines.append("")
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--json", type=Path)
    parser.add_argument("--markdown", type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    json_path = args.json or root / MIGRATION_DIR / "inventory.json"
    markdown_path = args.markdown or root / MIGRATION_DIR / "INVENTORY.md"
    try:
        hits = scan(root)
    except (RuntimeError, subprocess.CalledProcessError) as exc:
        print(f"inventory: error: {exc}", file=sys.stderr)
        return 1
    payload = {
        "schema": 1,
        "repository_root": str(root),
        "hit_count": len(hits),
        "hits": [asdict(hit) for hit in hits],
    }
    json_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    markdown_path.write_text(render_markdown(hits, root))
    print(f"inventory: {len(hits)} hits")
    print(f"inventory: wrote {json_path}")
    print(f"inventory: wrote {markdown_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
