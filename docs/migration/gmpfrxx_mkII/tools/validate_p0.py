#!/usr/bin/env python3
"""Validate committed-data portions of the MPLAPACK P0 migration gate."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
MIGRATION = ROOT / "docs/migration/gmpfrxx_mkII"
BASE = "b875e74d4b927282c907c3a29e6cadda48a7d57b"
PROMPT_HASH = "c4175afa3ec59c385f00a54c7b762b09c0f05b7c8418b33e940ce340d46d2c0e"
PROTOTYPE_HASH = "d14b03964672da96bb73c86960913d9a566a209bcba86379083bd71c56f262a4"
ARCHIVE_HASH = "c0816b3538b6b77009f714bb391cebe11abb2fdb69e07aa3bb305ff822764afb"
BACKENDS = {"gmp", "mpfr", "qd", "dd", "double", "binary80", "binary128"}
REQUIRED = {
    "REAL-DOUBLE", "REAL-DD", "REAL-QD", "REAL-BINARY80",
    "REAL-BINARY128", "REAL-GMP", "COMPLEX-DOUBLE", "COMPLEX-DD",
    "COMPLEX-QD", "COMPLEX-BINARY80", "COMPLEX-BINARY128", "COMPLEX-GMP",
}


class ValidationError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValidationError(message)


def load(path: Path) -> dict:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ValidationError(f"cannot load {path}: {exc}") from exc


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    try:
        with path.open("rb") as stream:
            for chunk in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(chunk)
    except OSError as exc:
        raise ValidationError(f"cannot hash {path}: {exc}") from exc
    return digest.hexdigest()


def validate_lock() -> None:
    lock = load(MIGRATION / "LOCK.json")
    mplapack = lock.get("mplapack", {})
    upstream = lock.get("gmpfrxx_mkII", {})
    prompt = lock.get("prompt_pack", {})
    prototype = lock.get("prototype", {})
    precision = lock.get("precision_policy", {})
    require(lock.get("schema") == 1, "LOCK schema is not 1")
    require(mplapack.get("base_sha") == BASE, "locked MPLAPACK base mismatch")
    require(re.fullmatch(r"[0-9a-f]{40}", upstream.get("adapter_base_sha", "")) is not None,
            "gmpfrxx_mkII adapter base is not a full SHA")
    require(upstream.get("bootstrap_sha256") == ARCHIVE_HASH,
            "bootstrap archive hash mismatch")
    require(prompt.get("version") == "5.4" and prompt.get("sha256") == PROMPT_HASH,
            "canonical prompt lock mismatch")
    require(sha256(ROOT / prompt["path"]) == PROMPT_HASH,
            "canonical prompt content changed")

    source = Path(prototype.get("source_path", ""))
    imported = ROOT / prototype.get("imported_path", "")
    require(prototype.get("source_filename") == "FABLE_PROTOTYPE_PROMPT.v4.md",
            "prototype source filename mismatch")
    require(prototype.get("source_sha256") == PROTOTYPE_HASH,
            "prototype source hash mismatch")
    require(prototype.get("imported_sha256") == PROTOTYPE_HASH,
            "prototype import hash mismatch")
    require(prototype.get("byte_for_byte_identical") is True,
            "prototype identity result is not true")
    require(source.is_file() and sha256(source) == PROTOTYPE_HASH,
            "historical prototype source unavailable or changed")
    require(sha256(imported) == PROTOTYPE_HASH and source.read_bytes() == imported.read_bytes(),
            "historical prototype import is not byte-for-byte identical")

    require(precision.get("default_bits") == 512,
            "default precision is not 512")
    require(precision.get("assignment_precision_policy") == "unchecked",
            "assignment precision policy changed")
    require(precision.get("interop_precision_parameters") == "none",
            "interop precision parameters changed")
    require(precision.get("runtime_compatibility_status")
            == "accepted-retire-variable-precision-behavior",
            "maintainer precision decision is not locked")


def validate_matrix() -> set[str]:
    with (MIGRATION / "interop_requirements.tsv").open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        require(reader.fieldnames == [
            "id", "source_type", "target_type", "form", "requirement", "evidence", "test_id"
        ], "interop matrix columns changed")
        rows = list(reader)
    require(len(rows) == len({row["id"] for row in rows}),
            "interop matrix contains duplicate IDs")
    require(all(row["requirement"] in {"REQUIRED", "FORBIDDEN", "OUT_OF_SCOPE", "NOT_USED"}
                for row in rows), "interop matrix has an unknown requirement")
    required = {row["id"] for row in rows if row["requirement"] == "REQUIRED"}
    require(required == REQUIRED, f"REQUIRED interop set mismatch: {sorted(required ^ REQUIRED)}")
    for row in rows:
        if row["requirement"] != "REQUIRED":
            continue
        require(bool(row["evidence"].strip()) and bool(row["test_id"].strip()),
                f"{row['id']} lacks evidence or test ID")
        for citation in row["evidence"].split(";"):
            evidence_path = citation.rsplit(":", 1)[0]
            require((ROOT / evidence_path).is_file(),
                    f"{row['id']} evidence path missing: {evidence_path}")
    require(not any(row["requirement"] == "REQUIRED" and
                    (row["id"].startswith("REVERSE-") or row["id"] in {
                        "MIXED-ARITHMETIC", "BOTH-OPERAND-ORDERS",
                        "COMPOUND-ASSIGN", "EXPR-COMPOSITION"}) for row in rows),
            "forbidden reverse or mixed operation is REQUIRED")
    require({row["id"] for row in rows if row["requirement"] == "OUT_OF_SCOPE"}
            == {"REAL-EDD", "COMPLEX-EDD", "REAL-TD", "COMPLEX-TD"},
            "edd/td scope rows changed")
    return required


def validate_spikes(required: set[str]) -> None:
    payload = load(MIGRATION / "spike/results.json")
    results = payload.get("results", [])
    by_id = {result.get("id"): result for result in results}
    require(payload.get("schema") == 1 and payload.get("default_precision_bits") == 512,
            "spike schema or precision changed")
    require(payload.get("semantic_macros_enabled") == [],
            "correctness-affecting semantic macro was enabled")
    require(len(by_id) == len(results) and required <= set(by_id),
            "spike result IDs are duplicate or incomplete")
    smoke = by_id.get("default-precision-512", {})
    require(smoke.get("must_pass_p0") is True and smoke.get("outcome") == "PASS"
            and smoke.get("compile_exit") == 0 and smoke.get("run_exit") == 0,
            "fresh-process 512-bit smoke did not pass")
    report = (MIGRATION / "SPIKE.md").read_text(encoding="utf-8")
    for result in results:
        require((MIGRATION / "tools/spike" / result.get("source", "")).is_file(),
                f"spike source missing: {result.get('source')}")
        require((MIGRATION / "spike" / result.get("log", "")).is_file(),
                f"spike log missing: {result.get('log')}")
        if result.get("outcome") != "PASS":
            require(bool(result.get("later_owner")),
                    f"{result.get('id')} failure has no owner")
            require(result.get("id", "") in report and result.get("later_owner", "") in report,
                    f"{result.get('id')} failure/owner absent from SPIKE.md")
    queried = []
    for source in sorted((MIGRATION / "tools/spike").glob("*.cpp")):
        text = source.read_text(encoding="utf-8")
        require("set_default" not in text and "mpf_set_default_prec" not in text,
                f"precision setter found in {source.name}")
        if re.search(r"\bget_prec\s*\(|\bprecision\s*\(", text):
            queried.append(source.name)
    require(queried == ["default_precision.cpp"],
            "precision metadata inspected outside the one default smoke")
    default = (MIGRATION / "tools/spike/default_precision.cpp").read_text(encoding="utf-8")
    require(default.count("get_prec() == 512") == 2,
            "default smoke does not check both wrappers at exactly 512 bits")


def validate_baseline() -> None:
    payload = load(MIGRATION / "baseline.json")
    scope = payload.get("test_scope", {})
    results = payload.get("results", [])
    require(payload.get("schema") == 1 and payload.get("mplapack_base_sha") == BASE,
            "baseline schema or base mismatch")
    require(payload.get("default_precision_bits") == 512,
            "baseline precision is not 512")
    require(set(scope.get("backends", [])) == {"gmp", "mpfr"}
            and scope.get("mpfr_profiles") == ["default"],
            "baseline backend/profile scope changed")
    keys = [(result.get("test_id"), result.get("metric")) for result in results]
    require(bool(results) and len(keys) == len(set(keys)),
            "baseline is empty or has duplicate keys")
    require(all(result.get("status") == "PASS" for result in results),
            "baseline contains a non-PASS result")
    for result in results:
        require(all(result.get(field) not in (None, "") for field in (
            "test_id", "metric", "unit", "direction", "value", "status", "source_line")),
            "baseline result has an empty field")
    try:
        with gzip.open(MIGRATION / "logs/baseline-raw.log.gz", "rt", encoding="utf-8") as stream:
            raw = stream.read()
    except (OSError, UnicodeError) as exc:
        raise ValidationError(f"cannot read baseline raw archive: {exc}") from exc
    require(raw.count("===== BEGIN ") == payload.get("raw_source_count")
            and raw.count("===== END ") == payload.get("raw_source_count"),
            "baseline raw archive is incomplete")


def validate_abi() -> None:
    payload = load(MIGRATION / "abi/manifest.json")
    backends = payload.get("backends", {})
    require(payload.get("schema") == 1 and payload.get("mplapack_base_sha") == BASE,
            "ABI manifest schema or base mismatch")
    require(set(backends) == BACKENDS, "ABI backend set is incomplete")
    require(payload.get("platform", {}).get("shared_library_format") == "ELF",
            "ABI platform format missing")
    require(bool(payload.get("cmake_package_dependencies")),
            "CMake package dependencies missing")
    for backend, record in backends.items():
        require(all(bool(record.get(field)) for field in
                    ("libraries", "public_headers", "pkg_config", "cmake_targets")),
                f"{backend}: ABI/install data incomplete")
        consumers = record.get("consumers", {})
        require(set(consumers) == {"pkg_config", "cmake"}
                and all(item.get("status") == "PASS" for item in consumers.values()),
                f"{backend}: downstream consumers incomplete or failing")
        for library in record["libraries"].values():
            require(bool(library.get("soname")) and re.fullmatch(
                r"[0-9a-f]{64}", library.get("sha256", "")) is not None,
                f"{backend}: library hash or SONAME missing")
            for kind in ("mangled_symbols", "demangled_symbols"):
                symbols = library.get(kind, {})
                path = MIGRATION / "abi" / symbols.get("file", "")
                require(path.is_file() and symbols.get("count", 0) > 0
                        and sha256(path) == symbols.get("sha256"),
                        f"{backend}: {kind} incomplete or changed")


def main() -> int:
    try:
        validate_lock()
        required = validate_matrix()
        validate_spikes(required)
        validate_baseline()
        validate_abi()
    except (OSError, UnicodeError, ValidationError) as exc:
        print(f"validate-p0: FAIL: {exc}", file=sys.stderr)
        return 1
    print("validate-p0: PASS")
    print("  lock/provenance: complete")
    print("  interop matrix: 12 REQUIRED; reverse/mixed forbidden; edd/td out")
    print("  precision smoke: one fresh process; MPFR=512; MPF=512")
    print("  baseline: GMP/MPFR Autotools sets; MPFR default profile only")
    print("  ABI/install: seven backends; pkg-config and CMake consumers pass")
    print("  compile spike: every observed failure has a later-phase owner")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
