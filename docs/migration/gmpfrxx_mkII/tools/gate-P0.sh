#!/usr/bin/env bash
set -euo pipefail

readonly base_sha=b875e74d4b927282c907c3a29e6cadda48a7d57b
readonly branch=topic/gmpfrxx_mkII_migration
readonly migration_rel=docs/migration/gmpfrxx_mkII

root=$(git rev-parse --show-toplevel)
cd "$root"
migration="$root/$migration_rel"

fail() {
    printf 'gate-P0: FAIL: %s\n' "$*" >&2
    exit 1
}

pass() {
    printf 'gate-P0: PASS: %s\n' "$*"
}

test "$(git branch --show-current)" = "$branch" ||
    fail "current branch is not $branch"
git merge-base --is-ancestor "$base_sha" HEAD ||
    fail "locked P0 base is not an ancestor of HEAD"
test "$(git merge-base "$base_sha" HEAD)" = "$base_sha" ||
    fail "branch ancestry does not preserve the locked P0 base"
pass "branch and immutable base"

mapfile -t changed_paths < <(
    {
        git diff --name-only "$base_sha"...HEAD
        git diff --name-only
        git diff --cached --name-only
        git ls-files --others --exclude-standard
    } | sort -u
)
test "${#changed_paths[@]}" -gt 0 || fail "no P0 deliverables found"
for path in "${changed_paths[@]}"; do
    case "$path" in
        "$migration_rel"/*) ;;
        *) fail "path outside migration documentation subtree changed: $path" ;;
    esac
done
pass "changed-path restriction (${#changed_paths[@]} paths)"

python3 "$migration/tools/validate_p0.py"
pass "locked data, interop, precision, baseline, ABI, and spike manifests"

temporary=$(mktemp -d "${TMPDIR:-/tmp}/mplapack-gate-p0.XXXXXX")
trap 'rm -rf "$temporary"' EXIT
python3 "$migration/tools/inventory.py" \
    --root "$root" \
    --json "$temporary/inventory.json" \
    --markdown "$temporary/INVENTORY.md"
cmp "$temporary/inventory.json" "$migration/inventory.json" ||
    fail "inventory.json is not reproducible"
cmp "$temporary/INVENTORY.md" "$migration/INVENTORY.md" ||
    fail "INVENTORY.md is not reproducible"
python3 -m unittest -v "$migration/tools/test_inventory.py"
pass "inventory reproduction and fixture test"

python3 -m unittest -v \
    "$migration/tools/test_capture_baseline.py" \
    "$migration/tools/test_compare_baseline.py"
python3 "$migration/tools/compare_baseline.py" \
    "$migration/baseline.json" \
    "$migration/baseline.json" \
    --rules "$migration/baseline_rules.tsv"
pass "baseline comparator tests and self-comparison"

python3 -m unittest -v "$migration/tools/test_compare_abi.py"
python3 "$migration/tools/compare_abi.py" \
    "$migration/abi/manifest.json" \
    "$migration/abi/manifest.json" \
    --phase P0
pass "ABI comparator tests and self-comparison"

: "${MPLAPACK_P0_PREFIX:=/tmp/mplapack-p0-final}"
: "${MPLAPACK_P0_SOURCE:=/tmp/mplapack-p0-base}"
: "${GMPFRXX_MKII_P0_PREFIX:=/tmp/gmpfrxx-p0-install}"
test -d "$MPLAPACK_P0_PREFIX" || fail "missing Autotools prefix: $MPLAPACK_P0_PREFIX"
test -d "$MPLAPACK_P0_SOURCE" || fail "missing locked source: $MPLAPACK_P0_SOURCE"
test -d "$GMPFRXX_MKII_P0_PREFIX" ||
    fail "missing gmpfrxx_mkII prefix: $GMPFRXX_MKII_P0_PREFIX"
python3 "$migration/tools/spike/run_spikes.py" \
    --mplapack-prefix "$MPLAPACK_P0_PREFIX" \
    --mplapack-source "$MPLAPACK_P0_SOURCE" \
    --gmpfrxx-prefix "$GMPFRXX_MKII_P0_PREFIX" \
    --work-dir "$temporary/spikes" \
    --results "$temporary/spike-results.json"
cmp "$temporary/spike-results.json" "$migration/spike/results.json" ||
    fail "compile-spike results are not reproducible"
pass "compile spikes and fresh-process 512-bit smoke"

test -s "$migration/REPORT-P0.md" || fail "REPORT-P0.md is missing"
test -s "$migration/CMAKE-TEST-LEVEL.md" ||
    fail "CMake test-level record is missing"
pass "P0 report and separate CMake test-level record"

printf 'gate-P0: PASS\n'
