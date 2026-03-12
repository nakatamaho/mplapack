#!/usr/bin/env python3
"""
summarize_mplapack_tests.py
============================
Parse MPLAPACK *.out test result files and produce a structured summary.

Precision naming rules
-----------------------
  - Non-mpfr backends       -> precision = backend name  (e.g. "double", "qd")
  - mpfr with variant       -> precision = "mpfr (binary128)" / "mpfr (binary64)"
  - mpfr with "default"     -> precision = "mpfr (default)"

Supported output formats: text (default), csv, json
"""

import argparse
import csv
import io
import json
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import NamedTuple, List, Dict, Tuple

# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------

class TestRecord(NamedTuple):
    category:  str   # "eig" or "lin"
    precision: str   # e.g. "double", "mpfr (binary64)", "qd", ...
    file:      str   # path relative to <root>
    suite:     str   # e.g. "DSB", "DGS drivers", ...
    failed:    int
    total:     int
    status:    str   # "PASS", "FAIL", or "UNKNOWN"

    @property
    def fail_rate(self) -> float:
        return self.failed / self.total if self.total > 0 else 0.0


# ---------------------------------------------------------------------------
# Regex patterns
# ---------------------------------------------------------------------------

# Pattern A - passed:
#   "All tests for DSB passed the threshold (   810 tests run)"
#   "All tests for DST drivers passed the threshold ( 13464 tests run)"
RE_PASS = re.compile(
    r"All\s+tests\s+for\s+(.+?)\s+passed\s+the\s+threshold\s*\(\s*(\d+)\s+tests\s+run\)",
    re.IGNORECASE,
)

# Pattern B - failed:
#   "DES:    2 out of  3822 tests failed to pass the threshold"
#   "DGS drivers:      1 out of   1555 tests failed to pass the threshold"
#
# NOTE: suite name must start with a letter ([A-Za-z]).
# When a test suite fails, LAPACK also emits per-matrix-type lines such as:
#   "  2:   1 out of    1 tests failed to pass the threshold"
# where "2" is a matrix-type index, not a suite name.  Anchoring the capture
# group to start with [A-Za-z] rejects those bare-number lines while still
# matching all real suite names (DES, ZST drivers, ZBD routines, etc.).
RE_FAIL = re.compile(
    r"^[ \t]*([A-Za-z][\w\s]*?):\s+(\d+)\s+out\s+of\s+(\d+)\s+tests\s+failed\s+to\s+pass\s+the\s+threshold",
    re.IGNORECASE | re.MULTILINE,
)

# mpfr variant suffix: <stem>.<variant>.out
RE_MPFR_VARIANT = re.compile(r"^.+\.(default|binary64|binary128)\.out$", re.IGNORECASE)


# ---------------------------------------------------------------------------
# Path helpers
# ---------------------------------------------------------------------------

def detect_precision(out_path: Path, search_root: Path) -> str:
    """
    Derive a single precision label from an .out file path.

    New layout:  <search_root>/results/<outdir>/<backend>/<file>.out
    Old layout:  <search_root>/<backend>/<outdir>/<file>.out
    For mpfr the variant is embedded in the filename suffix.
    """
    try:
        rel = out_path.relative_to(search_root)
    except ValueError:
        return "unknown"

    # New layout: results / <outdir> / <backend> / <file>.out  -> parts[2]
    # Old layout: <backend> / <outdir> / <file>.out            -> parts[0]
    if len(rel.parts) >= 3 and rel.parts[0] == "results":
        backend = rel.parts[2]
    else:
        backend = rel.parts[0] if rel.parts else "unknown"

    if backend != "mpfr":
        return backend

    # mpfr: extract variant from filename
    m = RE_MPFR_VARIANT.match(out_path.name)
    variant = m.group(1) if m else "default"
    return f"mpfr ({variant})"


def _infer_category(search_root: Path) -> str:
    """
    Infer eig/lin category from the resolved path of the search root.
    Falls back to the directory's own name, then 'unknown'.
    """
    for part in reversed(search_root.resolve().parts):
        if part in ("eig", "lin"):
            return part
    return search_root.resolve().name or "unknown"


# ---------------------------------------------------------------------------
# Parser
# ---------------------------------------------------------------------------

def parse_out_file(out_path: Path, category: str, search_root: Path) -> List[TestRecord]:
    """
    Parse a single .out file and return zero or more TestRecord entries.
    Returns one UNKNOWN record if no recognizable result lines are found.
    """
    precision = detect_precision(out_path, search_root)
    rel_file  = str(out_path.relative_to(search_root))

    try:
        text = out_path.read_text(errors="replace")
    except OSError as exc:
        print(f"[WARNING] Cannot read {out_path}: {exc}", file=sys.stderr)
        return [TestRecord(category, precision, rel_file, "READ_ERROR", 0, 0, "UNKNOWN")]

    records: List[TestRecord] = []

    # Collect PASS matches
    for m in RE_PASS.finditer(text):
        suite = m.group(1).strip()
        total = int(m.group(2))
        records.append(TestRecord(category, precision, rel_file, suite, 0, total, "PASS"))

    # Collect FAIL matches.
    # NOTE: Do NOT skip suites already seen in pass_suites.
    # A single .out file runs the same suite multiple times with different NB settings.
    # Some NB runs may PASS while others FAIL; every result line must be recorded
    # independently.  Filtering by pass_suites would silently drop FAIL records and
    # under-count both total_tests and total_failed.
    for m in RE_FAIL.finditer(text):
        suite  = m.group(1).strip()
        failed = int(m.group(2))
        total  = int(m.group(3))
        records.append(TestRecord(category, precision, rel_file, suite, failed, total, "FAIL"))

    if not records:
        records.append(TestRecord(category, precision, rel_file, "UNKNOWN", 0, 0, "UNKNOWN"))

    return records


def collect_all_records(outdir: str, search_root: Path) -> List[TestRecord]:
    """
    Walk all backend subdirectories and parse every *.out file.

    New layout (preferred):
        <search_root>/results/<outdir>/<backend>/*.out

    Old layout (fallback):
        <search_root>/<backend>/<outdir>/*.out
    """
    all_records: List[TestRecord] = []
    category = _infer_category(search_root)

    found_any = False

    # --- New layout: results/<outdir>/<backend>/*.out ---
    new_root = search_root / "results" / outdir
    if new_root.is_dir():
        for backend_dir in sorted(new_root.iterdir()):
            if not backend_dir.is_dir():
                continue
            for out_path in sorted(backend_dir.glob("*.out")):
                found_any = True
                all_records.extend(parse_out_file(out_path, category, search_root))

    # --- Old layout fallback: <backend>/<outdir>/*.out ---
    if not found_any:
        for backend_dir in sorted(search_root.iterdir()):
            if not backend_dir.is_dir() or backend_dir.name == "results":
                continue
            target_dir = backend_dir / outdir
            if not target_dir.is_dir():
                continue
            for out_path in sorted(target_dir.glob("*.out")):
                found_any = True
                all_records.extend(parse_out_file(out_path, category, search_root))

    if not found_any:
        print(
            f"[INFO] No *.out files found for OUTDIR='{outdir}' under: {search_root}\n"
            f"       Tried: {search_root}/results/{outdir}/<backend>/*.out\n"
            f"       and:  {search_root}/<backend>/{outdir}/*.out",
            file=sys.stderr,
        )

    return all_records


# ---------------------------------------------------------------------------
# Sorting
# ---------------------------------------------------------------------------

def sort_records(records: List[TestRecord]) -> List[TestRecord]:
    """Default sort: fail_rate descending, then category / precision / file / suite."""
    return sorted(
        records,
        key=lambda r: (-r.fail_rate, r.category, r.precision, r.file, r.suite),
    )


# ---------------------------------------------------------------------------
# Summary statistics
# ---------------------------------------------------------------------------

def build_summary(records: List[TestRecord]) -> Dict:
    """
    Aggregate statistics grouped by (category, precision).
    Returns an ordered dict keyed by (category, precision).
    """
    groups: Dict[Tuple[str, str], Dict] = defaultdict(
        lambda: {"total_tests": 0, "total_failed": 0, "files_with_fail": set(), "file_count": set()}
    )

    for r in records:
        key = (r.category, r.precision)
        g   = groups[key]
        g["total_tests"]  += r.total
        g["total_failed"] += r.failed
        # Count only files that produced at least one real test result.
        # UNKNOWN records (total=0) come from files with no recognizable output
        # (e.g. Cbak.out, Cbal.out) and should not inflate file_count.
        if r.status != "UNKNOWN":
            g["file_count"].add(r.file)
        if r.failed > 0:
            g["files_with_fail"].add(r.file)

    result = {}
    for key, g in sorted(groups.items()):
        result[key] = {
            "total_tests":     g["total_tests"],
            "total_failed":    g["total_failed"],
            "fail_rate":       g["total_failed"] / g["total_tests"] if g["total_tests"] > 0 else 0.0,
            "files_with_fail": len(g["files_with_fail"]),
            "file_count":      len(g["file_count"]),
        }
    return result


# ---------------------------------------------------------------------------
# Output formatters
# ---------------------------------------------------------------------------

DETAIL_HEADERS = ["category", "precision", "file", "suite",
                  "failed", "total", "fail_rate", "status"]


def format_text(records: List[TestRecord], summary: Dict) -> str:
    """Plain-text fixed-width table with dynamic column widths."""
    lines: List[str] = []

    # Dynamic widths
    prec_w  = max((len(r.precision) for r in records), default=10)
    prec_w  = max(prec_w, len("PRECISION"))
    file_w  = max((len(r.file)      for r in records), default=50)
    file_w  = max(file_w, len("FILE"))
    suite_w = max((len(r.suite)     for r in records), default=22)
    suite_w = max(suite_w, len("SUITE"))

    def row(cat, prec, file_, suite, failed, total, fail_pct, status):
        return (
            f"{cat:<8}  {prec:<{prec_w}}  {file_:<{file_w}}  "
            f"{suite:<{suite_w}}  {str(failed):>7}  {str(total):>9}  "
            f"{fail_pct:>8}  {status:<7}"
        )

    header = row("CATEGORY", "PRECISION", "FILE", "SUITE",
                 "FAILED", "TOTAL", "FAIL%", "STATUS")
    sep = "-" * len(header)

    lines += [sep, "MPLAPACK Test Summary - Detail", sep, header, sep]

    for r in records:
        fail_pct = f"{r.fail_rate * 100:.2f}%" if r.total > 0 else "N/A"
        lines.append(row(r.category, r.precision, r.file, r.suite,
                         r.failed, r.total, fail_pct, r.status))

    lines += [sep, f"Total rows: {len(records)}", ""]

    # --- Summary table ---
    prec_ws = max((len(k[1]) for k in summary), default=10)
    prec_ws = max(prec_ws, len("PRECISION"))

    def srow(cat, prec, tt, tf, fp, wf, fc):
        return (
            f"{cat:<8}  {prec:<{prec_ws}}  "
            f"{str(tt):>12}  {str(tf):>10}  {fp:>8}  {str(wf):>12}  {str(fc):>10}"
        )

    sheader = srow("CATEGORY", "PRECISION",
                   "TOTAL_TESTS", "TOTAL_FAIL", "FAIL%", "FILES_W_FAIL", "FILE_COUNT")
    ssep = "-" * len(sheader)

    lines += [ssep,
              "MPLAPACK Test Summary - Aggregated by (category, precision)",
              ssep, sheader, ssep]

    for (cat, prec), stats in summary.items():
        fp = f"{stats['fail_rate'] * 100:.2f}%"
        lines.append(srow(cat, prec,
                          stats["total_tests"], stats["total_failed"], fp,
                          stats["files_with_fail"], stats["file_count"]))

    lines.append(ssep)
    return "\n".join(lines)


def format_csv(records: List[TestRecord], summary: Dict) -> str:
    buf    = io.StringIO()
    writer = csv.writer(buf)

    writer.writerow(DETAIL_HEADERS)
    for r in records:
        writer.writerow([
            r.category, r.precision, r.file, r.suite,
            r.failed, r.total, f"{r.fail_rate:.6f}", r.status,
        ])

    writer.writerow([])
    writer.writerow(["# Summary by (category, precision)"])
    writer.writerow(["category", "precision", "total_tests", "total_failed",
                     "fail_rate", "files_with_fail", "file_count"])
    for (cat, prec), stats in summary.items():
        writer.writerow([
            cat, prec,
            stats["total_tests"], stats["total_failed"],
            f"{stats['fail_rate']:.6f}",
            stats["files_with_fail"], stats["file_count"],
        ])
    return buf.getvalue()


def format_json(records: List[TestRecord], summary: Dict) -> str:
    detail = [
        {
            "category":  r.category,
            "precision": r.precision,
            "file":      r.file,
            "suite":     r.suite,
            "failed":    r.failed,
            "total":     r.total,
            "fail_rate": round(r.fail_rate, 6),
            "status":    r.status,
        }
        for r in records
    ]
    agg = [
        {
            "category":  cat,
            "precision": prec,
            **{k: (round(v, 6) if isinstance(v, float) else v) for k, v in stats.items()},
        }
        for (cat, prec), stats in summary.items()
    ]
    return json.dumps({"detail": detail, "summary": agg}, indent=2)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Summarize MPLAPACK test *.out files into a structured report.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument(
        "outdir",
        metavar="OUTDIR",
        help=(
            "Output subdirectory name (e.g. Ryzen_Threadripper_3970X_gcc_14_2_0_ubuntu22_04). "
            "The script first tries SEARCH_ROOT/results/OUTDIR/<backend>/*.out (new layout), "
            "then falls back to SEARCH_ROOT/<backend>/OUTDIR/*.out (old layout)."
        ),
    )
    p.add_argument(
        "search_root",
        metavar="SEARCH_ROOT",
        nargs="?",
        default=".",
        help=(
            "The eig/ or lin/ directory. "
            "For the new layout this must contain a results/ subdirectory. "
            "Defaults to the current directory."
        ),
    )
    p.add_argument(
        "--only-fail",
        action="store_true",
        help="Show only records where status == FAIL.",
    )
    p.add_argument(
        "--format",
        choices=["text", "csv", "json"],
        default="text",
        help="Output format (default: text).",
    )
    p.add_argument(
        "--no-sort",
        action="store_true",
        help="Disable default fail_rate descending sort; keep discovery order.",
    )
    return p


def main() -> int:
    parser = build_parser()
    args   = parser.parse_args()

    search_root = Path(args.search_root).expanduser().resolve()
    if not search_root.is_dir():
        print(f"[ERROR] SEARCH_ROOT does not exist or is not a directory: {search_root}",
              file=sys.stderr)
        return 1

    records = collect_all_records(args.outdir, search_root)

    if not records:
        print(
            f"[WARNING] No records extracted. "
            f"Check that <backend>/{args.outdir}/*.out files exist under {search_root}.",
            file=sys.stderr,
        )
        return 0

    if not args.no_sort:
        records = sort_records(records)

    # Build summary from ALL records before any display filtering,
    # so aggregate totals are always accurate regardless of --only-fail.
    summary = build_summary(records)

    if args.only_fail:
        records = [r for r in records if r.status == "FAIL"]

    if args.format == "text":
        output = format_text(records, summary)
    elif args.format == "csv":
        output = format_csv(records, summary)
    else:
        output = format_json(records, summary)

    print(output)
    return 0


if __name__ == "__main__":
    sys.exit(main())


# =============================================================================
# README / Usage Examples
# =============================================================================
#
# BASIC USAGE  (run from within an eig/ or lin/ results-containing directory)
#   python summarize_mplapack_tests.py <OUTDIR>
#   python summarize_mplapack_tests.py Ryzen_Threadripper_3970X_gcc_14_2_0_ubuntu22_04
#
# WITH EXPLICIT SEARCH ROOT
#   python summarize_mplapack_tests.py <OUTDIR> /path/to/eig
#
# The script first looks for the new layout:
#   SEARCH_ROOT/results/<OUTDIR>/<backend>/*.out
# and falls back to the old layout if results/ is absent:
#   SEARCH_ROOT/<backend>/<OUTDIR>/*.out
#
# SHOW ONLY FAILURES
#   python summarize_mplapack_tests.py <OUTDIR> . --only-fail
#
# CSV / JSON
#   python summarize_mplapack_tests.py <OUTDIR> . --format csv  > results.csv
#   python summarize_mplapack_tests.py <OUTDIR> . --format json > results.json
#
# DISABLE SORT (keep filesystem discovery order)
#   python summarize_mplapack_tests.py <OUTDIR> . --no-sort
#
# =============================================================================
# SAMPLE AGGREGATED OUTPUT
# =============================================================================
#
# CATEGORY  PRECISION          TOTAL_TESTS  TOTAL_FAIL     FAIL%  FILES_W_FAIL  FILE_COUNT
# ------------------------------------------------------------------------------------
# eig       binary128              1232900           0     0.00%             0          40
# eig       binary80               1232900           0     0.00%             0          40
# eig       dd                     1206500           3     0.00%             1          40
# eig       double                 1228460           0     0.00%             0          40
# eig       gmp                    1205048           2     0.00%             1          38
# eig       mpfr (binary128)       1219436           0     0.00%             0          40
# eig       mpfr (binary64)        1231340           2     0.00%             1          40
# eig       mpfr (default)         1219436           0     0.00%             0          40
# eig       qd                     1206486          37     0.00%             3          40
# lin       binary128               865775          12     0.00%             2           4
# lin       binary80                865775           1     0.00%             1           4
# lin       dd                      865775          12     0.00%             2           4
# lin       double                  865775           1     0.00%             1           4
# lin       gmp                     865775           1     0.00%             1           4
# lin       mpfr (binary128)        865775          12     0.00%             2           4
# lin       mpfr (binary64)         865775           4     0.00%             1           4
# lin       mpfr (default)          865775          12     0.00%             2           4
# lin       qd                      865775          12     0.00%             2           4
# =============================================================================
