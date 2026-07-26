#!/usr/bin/env bash
set -euo pipefail

readonly p1_sha=2ebf798fd0081ecdc5c1b53fc117431c406bf884
readonly locked_base=b875e74d4b927282c907c3a29e6cadda48a7d57b
readonly branch=topic/gmpfrxx_mkII_migration
readonly upstream_repo=nakatamaho/gmpfrxx_mkII
readonly upstream_tag=v1.1.0
readonly upstream_tag_target=429fd1b35e1927ebaccc9fda5aa2801300b45bf5
readonly upstream_evidence_sha=24849d87a222c096cbbf3cdaa9a295c828be0aa3
readonly archive_name=gmpfrxx_mkII.1.1.0.tar.xz
readonly archive_size=15169540
readonly archive_sha=e0f3b813463b7a45dd493a818c60a17530075e0e647ea02227b75501c1984c73
readonly migration_rel=docs/migration/gmpfrxx_mkII
readonly vendored_rel=external/gmpfrxx_mkII/download/gmpfrxx_mkII.1.1.0.tar.xz
readonly bootstrap_version=1.0.1
readonly bootstrap_rel="external/gmpfrxx_mkII/download/gmpfrxx_mkII.${bootstrap_version}.tar.xz"

root=$(git rev-parse --show-toplevel)
cd "$root"

fail() {
    printf 'gate-P1R: FAIL: %s\n' "$*" >&2
    exit 1
}

pass() {
    printf 'gate-P1R: PASS: %s\n' "$*"
}

for tool in awk autoreconf cmake cmp find gh git grep make patch python3 rg sha256sum stat tar; do
    command -v "$tool" >/dev/null 2>&1 || fail "missing required tool: $tool"
done

jobs=${P1R_JOBS:-32}
case "$jobs" in
    ''|*[!0-9]*) fail "P1R_JOBS must be a positive integer" ;;
esac
test "$jobs" -gt 0 || fail "P1R_JOBS must be greater than zero"

START_HEAD=$(git rev-parse HEAD)
START_STATUS=$(git status --porcelain=v1 --untracked-files=all)
START_TREE=$(git rev-parse HEAD^{tree})
START_INDEX_TREE=$(git write-tree)
test -z "$START_STATUS" || fail "primary worktree must be clean"

RUNTIME_DIR=$(mktemp -d "${TMPDIR:-/tmp}/mplapack-gmpfrxx-p1r-gate.XXXXXX")
cleanup() {
    rm -rf "$RUNTIME_DIR"
}
trap cleanup EXIT INT TERM
exec > >(tee "$RUNTIME_DIR/gate-P1R.log") 2>&1

printf 'P1R runtime: %s\n' "$RUNTIME_DIR"
printf 'P1R jobs: %s\n' "$jobs"
printf 'P1R compiler: %s\n' "$(g++ --version | head -n 1)"
printf 'P1R CMake: %s\n' "$(cmake --version | head -n 1)"

test "$(git branch --show-current)" = "$branch" || fail "wrong branch"
git merge-base --is-ancestor "$locked_base" HEAD || fail "locked base is not an ancestor"
git merge-base --is-ancestor "$p1_sha" HEAD || fail "approved P1 is not an ancestor"
test "$(git merge-base "$locked_base" HEAD)" = "$locked_base" || fail "locked-base merge base changed"

mapfile -t changed_paths < <(git diff --name-only "$p1_sha"..HEAD | LC_ALL=C sort -u)
test "${#changed_paths[@]}" -gt 0 || fail "no P1R changes found"
for path in "${changed_paths[@]}"; do
    case "$path" in
        Makefile.am | \
        external/distfiles.sha256 | \
        external/gmpfrxx_mkII/Makefile.am | \
        external/gmpfrxx_mkII/download/gmpfrxx_mkII.1.0.1.tar.xz | \
        external/gmpfrxx_mkII/download/gmpfrxx_mkII.1.1.0.tar.xz | \
        "$migration_rel/P1R_RELEASE.json" | \
        "$migration_rel/P1R_HISTORICAL_REFERENCES.tsv" | \
        "$migration_rel/REPORT-P1R.md" | \
        "$migration_rel/tools/gate-P1R.sh") ;;
        *) fail "path outside P1R scope changed: $path" ;;
    esac
done
pass "branch, ancestry, and ${#changed_paths[@]} classified P1R paths"

immutable_paths=(
    "$migration_rel/CODEX_PROMPTS.md"
    "$migration_rel/LOCK.json"
    "$migration_rel/PRECONDITIONS.md"
    "$migration_rel/REPORT-P0.md"
    "$migration_rel/REPORT-P1.md"
    "$migration_rel/baseline.json"
    "$migration_rel/baseline_rules.tsv"
    "$migration_rel/abi/manifest.json"
    "$migration_rel/interop_requirements.tsv"
    "$migration_rel/SPIKE.md"
    "$migration_rel/spike/results.json"
    "$migration_rel/tools/gate-P0.sh"
    "$migration_rel/tools/gate-P1.sh"
)
for path in "${immutable_paths[@]}"; do
    test "$(git rev-parse "$p1_sha:$path")" = "$(git rev-parse "HEAD:$path")" ||
        fail "immutable P0/P1 artifact changed: $path"
done
git diff --quiet "$p1_sha"..HEAD -- include mpblas mplapack fable ||
    fail "backend, numerical, generated, or Fable path changed"
pass "immutable P0/P1 provenance, baseline, ABI, interop, spike, and source paths"

test -f "$vendored_rel" || fail "missing vendored 1.1.0 archive"
test ! -e "$bootstrap_rel" || fail "bootstrap archive remains actively vendored"
test "$(stat -c '%s' "$vendored_rel")" = "$archive_size" || fail "vendored archive size mismatch"
test "$(sha256sum "$vendored_rel" | awk '{print $1}')" = "$archive_sha" ||
    fail "vendored archive digest mismatch"
grep -Fx "$archive_sha  gmpfrxx_mkII/download/$archive_name" external/distfiles.sha256 >/dev/null ||
    fail "external/distfiles.sha256 lacks authoritative archive entry"
test "$(grep -Fc 'gmpfrxx_mkII/download/' external/distfiles.sha256)" = 1 ||
    fail "external/distfiles.sha256 contains multiple gmpfrxx_mkII archives"
grep -Fx 'GMPFRXX_MKII_VERSION = 1.1.0' external/gmpfrxx_mkII/Makefile.am >/dev/null ||
    fail "active external package version is not 1.1.0"
grep -F "external/gmpfrxx_mkII/download/$archive_name" Makefile.am >/dev/null ||
    fail "top-level distribution list lacks 1.1.0 archive"
test -z "$(find external/gmpfrxx_mkII/patches -type f -print 2>/dev/null || true)" ||
    fail "gmpfrxx_mkII patches directory is not empty"
pass "active vendored archive identity and empty patch set"

python3 - "$migration_rel/P1R_RELEASE.json" <<'PY'
import json
import sys

with open(sys.argv[1], encoding="utf-8") as handle:
    data = json.load(handle)
expected = {
    "version": "1.1.0",
    "tag": "v1.1.0",
    "tag_target_sha": "429fd1b35e1927ebaccc9fda5aa2801300b45bf5",
    "upstream_evidence_sha": "24849d87a222c096cbbf3cdaa9a295c828be0aa3",
    "archive_filename": "gmpfrxx_mkII.1.1.0.tar.xz",
    "archive_size_bytes": 15169540,
    "archive_sha256": "e0f3b813463b7a45dd493a818c60a17530075e0e647ea02227b75501c1984c73",
    "release_draft": False,
    "release_prerelease": False,
    "published_asset_verified": True,
}
if data != expected:
    raise SystemExit(f"P1R_RELEASE.json mismatch: {data!r}")
print("P1R release JSON: PASS")
PY

python3 - "$root" "$migration_rel/P1R_HISTORICAL_REFERENCES.tsv" <<'PY'
import pathlib
import re
import subprocess
import sys

root = pathlib.Path(sys.argv[1])
whitelist_path = root / sys.argv[2]
lines = whitelist_path.read_text(encoding="utf-8").splitlines()
if not lines or lines[0] != "path\treason":
    raise SystemExit("invalid historical-reference whitelist header")
whitelist = {}
for line in lines[1:]:
    path, reason = line.split("\t", 1)
    if not reason.strip() or path in whitelist:
        raise SystemExit(f"invalid whitelist row: {line!r}")
    if not (root / path).is_file():
        raise SystemExit(f"whitelisted path is not a file: {path}")
    whitelist[path] = reason

tracked = subprocess.check_output(
    ["git", "-C", str(root), "ls-files", "-z"]
).split(b"\0")
patterns = (
    re.compile(rb"gmpfrxx_mkII\.1\.0\.1\.tar\.xz"),
    re.compile(rb"c0816b3538b6b77009f714bb391cebe11abb2fdb69e07aa3bb305ff822764afb"),
    re.compile(rb"GMPFRXX_MKII_VERSION[ \t]*=[ \t]*1\.0\.1"),
    re.compile(rb"gmpfrxx_mkII.{0,96}\b1\.0\.1\b", re.IGNORECASE),
)
historical_version = re.compile(rb"\bv?1\.0\.1\b")
matches = {}
for raw_path in tracked:
    if not raw_path:
        continue
    path = raw_path.decode("utf-8", "surrogateescape")
    data = (root / path).read_bytes()
    if b"\0" in data:
        continue
    hit_count = sum(len(pattern.findall(data)) for pattern in patterns)
    if path.startswith("docs/migration/gmpfrxx_mkII/"):
        hit_count += len(historical_version.findall(data))
    if hit_count:
        matches[path] = hit_count
unexpected = sorted(set(matches) - set(whitelist))
if unexpected:
    raise SystemExit("active stale bootstrap references: " + ", ".join(unexpected))
unused = sorted(set(whitelist) - set(matches))
if unused:
    raise SystemExit("whitelist entries without a matching historical reference: " + ", ".join(unused))
print("historical bootstrap references:")
for path in sorted(matches):
    print(f"  {path}: {matches[path]}")
print("active stale-reference scan: PASS")
PY
pass "P1R provenance and explicit active-reference gate"

gh auth status
gh api "repos/$upstream_repo/releases/tags/$upstream_tag" > "$RUNTIME_DIR/release.json"
python3 - "$RUNTIME_DIR/release.json" <<'PY'
import json
import sys

with open(sys.argv[1], encoding="utf-8") as handle:
    release = json.load(handle)
assert release["tag_name"] == "v1.1.0"
assert release["name"] == "gmpfrxx_mkII 1.1.0"
assert release["draft"] is False
assert release["prerelease"] is False
assets = [asset for asset in release.get("assets", [])
          if asset.get("name") == "gmpfrxx_mkII.1.1.0.tar.xz"]
assert len(release.get("assets", [])) == 1
assert len(assets) == 1
assert assets[0]["size"] == 15169540
assert "e0f3b813463b7a45dd493a818c60a17530075e0e647ea02227b75501c1984c73" in release.get("body", "")
print(f"published release id={release['id']} asset id={assets[0]['id']}")
PY

tag_lines=$(git ls-remote --tags "https://github.com/$upstream_repo.git" 'refs/tags/v1.1.0*')
tag_object=$(printf '%s\n' "$tag_lines" | awk '$2 == "refs/tags/v1.1.0" {print $1}')
tag_target=$(printf '%s\n' "$tag_lines" | awk '$2 == "refs/tags/v1.1.0^{}" {print $1}')
test -n "$tag_object" || fail "upstream annotated tag object is missing"
test "$tag_target" = "$upstream_tag_target" || fail "upstream tag target mismatch"
evidence_sha=$(git ls-remote --heads "https://github.com/$upstream_repo.git" \
    topic/mplapack-release-hardening | awk '{print $1}')
test "$evidence_sha" = "$upstream_evidence_sha" || fail "upstream evidence branch mismatch"

download_dir="$RUNTIME_DIR/published-download"
mkdir -p "$download_dir"
gh release download "$upstream_tag" --repo "$upstream_repo" \
    --pattern "$archive_name" --dir "$download_dir"
published_archive="$download_dir/$archive_name"
test -f "$published_archive" || fail "published archive download is missing"
test "$(stat -c '%s' "$published_archive")" = "$archive_size" || fail "published archive size mismatch"
test "$(sha256sum "$published_archive" | awk '{print $1}')" = "$archive_sha" ||
    fail "published archive digest mismatch"
cmp --silent "$published_archive" "$vendored_rel" || fail "published and vendored archives differ"
tar -tf "$published_archive" > "$RUNTIME_DIR/published-archive-files.txt"
if grep -E '(^|/)\.git(/|$)|/home/docker|prototype-snapshot|forensic' \
    "$RUNTIME_DIR/published-archive-files.txt"; then
    fail "published archive contains a forbidden path"
fi
pass "published release, annotated tag $tag_object, evidence branch, and byte-identical asset"

verify_archive() {
    local archive=$1
    test "$(sha256sum "$archive" | awk '{print $1}')" = "$archive_sha"
}

corruption_dir="$RUNTIME_DIR/corruption"
mkdir -p "$corruption_dir"
corrupt_archive="$corruption_dir/$archive_name"
cp "$vendored_rel" "$corrupt_archive"
printf X >> "$corrupt_archive"
if verify_archive "$corrupt_archive"; then
    fail "corrupted archive passed digest verification"
fi
cp "$published_archive" "$corrupt_archive"
verify_archive "$corrupt_archive" || fail "restored archive failed verification"
cmp --silent "$published_archive" "$corrupt_archive" || fail "archive restoration was not exact"
pass "corruption rejected and exact published archive restored in disposable storage"

source_tree="$RUNTIME_DIR/source"
mkdir -p "$source_tree"
git archive --format=tar HEAD | tar -xf - -C "$source_tree"
cmp --silent "$vendored_rel" "$source_tree/$vendored_rel" ||
    fail "committed source snapshot contains a different archive"

(
    cd "$source_tree"
    autoreconf --force --install
)
pass "disposable Autotools source regenerated with autoreconf --force --install"

export http_proxy=http://127.0.0.1:9
export https_proxy=http://127.0.0.1:9
export HTTP_PROXY=http://127.0.0.1:9
export HTTPS_PROXY=http://127.0.0.1:9
export ALL_PROXY=http://127.0.0.1:9
export NO_PROXY=

autotools_build="$RUNTIME_DIR/autotools-build"
deps_src="$RUNTIME_DIR/dependency-source"
deps_build="$RUNTIME_DIR/dependency-build"
deps_prefix="$autotools_build/external/i"
mkdir -p "$autotools_build" "$deps_src" "$deps_build" "$deps_prefix"

build_dependency() {
    local name=$1
    local archive=$2
    local source_subdir=$3
    local prefix=$4
    shift 4
    local source="$deps_src/$source_subdir"
    local build="$deps_build/$source_subdir"
    tar -xf "$archive" -C "$deps_src"
    if test "$name" = GMP; then
        patch -d "$source" -p0 < "$source_tree/external/gmp/patches/gcc-15.patch"
    fi
    mkdir -p "$build"
    (
        cd "$build"
        "$source/configure" --prefix="$prefix" "$@"
        make -j"$jobs"
        make install
    )
    test -d "$prefix/include" || fail "$name staged include directory missing"
    test -d "$prefix/lib" || fail "$name staged library directory missing"
    pass "$name staged from committed archive"
}

build_dependency GMP \
    "$source_tree/external/gmp/download/gmp-6.3.0.tar.xz" \
    gmp-6.3.0 "$deps_prefix/GMP" --enable-cxx
export LD_LIBRARY_PATH="$deps_prefix/GMP/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
build_dependency MPFR \
    "$source_tree/external/mpfr/download/mpfr-4.2.2.tar.bz2" \
    mpfr-4.2.2 "$deps_prefix/MPFR" \
    --with-gmp-include="$deps_prefix/GMP/include" \
    --with-gmp-lib="$deps_prefix/GMP/lib"
export LD_LIBRARY_PATH="$deps_prefix/MPFR/lib:$LD_LIBRARY_PATH"
build_dependency MPC \
    "$source_tree/external/mpc/download/mpc-1.4.1.tar.xz" \
    mpc-1.4.1 "$deps_prefix/MPC" \
    --with-gmp-include="$deps_prefix/GMP/include" \
    --with-gmp-lib="$deps_prefix/GMP/lib" \
    --with-mpfr-include="$deps_prefix/MPFR/include" \
    --with-mpfr-lib="$deps_prefix/MPFR/lib"
export LD_LIBRARY_PATH="$deps_prefix/MPC/lib:$LD_LIBRARY_PATH"

final_prefix="$RUNTIME_DIR/final-prefix"
(
    cd "$autotools_build"
    "$source_tree/configure" \
        --prefix="$final_prefix" \
        --disable-dependency-tracking \
        --enable-gmpfrxx-mkII=yes \
        --enable-gmp=no \
        --enable-mpfr=no \
        --enable-qd=no \
        --enable-dd=no \
        --enable-double=no \
        --enable-binary80=no \
        --enable-binary128=no \
        --enable-test=no \
        --enable-benchmark=no \
        --disable-examples
)

make -C "$autotools_build/external/gmpfrxx_mkII" \
    GMPFRXX_MKII_BUILD_JOBS="$jobs" all
internal_prefix="$deps_prefix/GMPFRXX_MKII"
test -f "$internal_prefix/include/gmpfrxx_mkII.h" || fail "internal header missing"
test -f "$internal_prefix/lib/cmake/gmpfrxx_mkII/gmpfrxx_mkIIConfig.cmake" ||
    fail "internal package config missing"
pass "Autotools internal installation"

make -C "$autotools_build/external/gmpfrxx_mkII" \
    GMPFRXX_MKII_BUILD_JOBS="$jobs" install
test -f "$final_prefix/include/gmpfrxx_mkII.h" || fail "final header missing"
test -f "$final_prefix/lib/cmake/gmpfrxx_mkII/gmpfrxx_mkIIConfig.cmake" ||
    fail "final package config missing"
pass "Autotools final-prefix installation"

for cache in \
    "$autotools_build/external/gmpfrxx_mkII/work/internal/CMakeCache.txt" \
    "$autotools_build/external/gmpfrxx_mkII/work/install/CMakeCache.txt"; do
    grep -Fx 'GMPFRXX_MKII_DEPS_AUTO_FETCH:BOOL=OFF' "$cache" >/dev/null ||
        fail "dependency auto-fetch is not OFF in $cache"
done

for prefix in "$internal_prefix" "$final_prefix"; do
    cmake_dir="$prefix/lib/cmake/gmpfrxx_mkII"
    test -f "$cmake_dir/gmpfrxx_mkIIConfigVersion.cmake" || fail "version metadata missing"
    grep -F 'set(PACKAGE_VERSION "1.1.0")' "$cmake_dir/gmpfrxx_mkIIConfigVersion.cmake" >/dev/null ||
        fail "installed package version is not 1.1.0"
    grep -F 'add_library(gmpfrxx_mkII::gmpfrxx_mkII INTERFACE IMPORTED)' \
        "$cmake_dir/gmpfrxx_mkIITargets.cmake" >/dev/null || fail "aggregate exported target missing"
    if rg -n -F "$RUNTIME_DIR" "$cmake_dir"; then
        fail "installed metadata contains a runtime build path"
    fi
    if rg -n -F "$source_tree" "$cmake_dir"; then
        fail "installed metadata contains a source path"
    fi
    if rg -n '(^|[^A-Za-z0-9_])(-lgmpxx|libgmpxx\.(so|a|dylib))([^A-Za-z0-9_]|$)' "$prefix"; then
        fail "unexpected libgmpxx dependency in installed prefix"
    fi
    find "$prefix" -type f -printf '%P\n' | LC_ALL=C sort > "$RUNTIME_DIR/$(basename "$prefix").manifest"
done
pass "installed manifests, version, targets, path hygiene, and no libgmpxx dependency"

independent_source="$RUNTIME_DIR/independent-source"
independent_build="$RUNTIME_DIR/independent-build"
independent_prefix="$RUNTIME_DIR/independent-prefix"
mkdir -p "$independent_source"
tar -xf "$source_tree/$vendored_rel" -C "$independent_source"
cmake -S "$independent_source/gmpfrxx_mkII.1.1.0" -B "$independent_build" \
    -DCMAKE_INSTALL_PREFIX="$independent_prefix" \
    -DBUILD_TESTING=OFF \
    -DGMPFRXX_MKII_DEPS_AUTO_FETCH=OFF \
    -DGMPFRXX_MKII_BUILD_EXAMPLES=OFF \
    -DGMPFRXX_MKII_BUILD_BENCHMARKS=OFF \
    -DGMPFRXX_MKII_COMPONENTS=GMP,MPFR,MPC \
    -DCMAKE_PREFIX_PATH="$deps_prefix/GMP;$deps_prefix/MPFR;$deps_prefix/MPC" \
    -DGMP_INCLUDE_DIR="$deps_prefix/GMP/include" \
    -DGMP_LIBRARY="$deps_prefix/GMP/lib/libgmp.so" \
    -DMPFR_INCLUDE_DIR="$deps_prefix/MPFR/include" \
    -DMPFR_LIBRARY="$deps_prefix/MPFR/lib/libmpfr.so" \
    -DMPC_INCLUDE_DIR="$deps_prefix/MPC/include" \
    -DMPC_LIBRARY="$deps_prefix/MPC/lib/libmpc.so"
cmake --build "$independent_build" -j"$jobs"
cmake --install "$independent_build"
grep -Fx 'GMPFRXX_MKII_DEPS_AUTO_FETCH:BOOL=OFF' "$independent_build/CMakeCache.txt" >/dev/null ||
    fail "independent install enabled dependency auto-fetch"
pass "independent published-archive installation"

cmake_discovery() {
    local name=$1
    local package_prefix=$2
    local build="$RUNTIME_DIR/mplapack-cmake-$name"
    cmake -S "$source_tree" -B "$build" \
        -DMPLAPACK_ENABLE_GMPFRXX_MKII=ON \
        -DMPLAPACK_ENABLE_GMP=OFF \
        -DMPLAPACK_ENABLE_MPFR=OFF \
        -DMPLAPACK_ENABLE_QD=OFF \
        -DMPLAPACK_ENABLE_DD=OFF \
        -DMPLAPACK_ENABLE_DOUBLE=ON \
        -DMPLAPACK_ENABLE_BINARY80=OFF \
        -DMPLAPACK_ENABLE_BINARY128=OFF \
        -DMPLAPACK_BUILD_TESTS=OFF \
        -DMPLAPACK_BUILD_EXAMPLES=OFF \
        -DMPLAPACK_BUILD_BENCHMARKS=OFF \
        -DCMAKE_PREFIX_PATH="$package_prefix;$deps_prefix/GMP;$deps_prefix/MPFR;$deps_prefix/MPC"
    grep -F 'gmpfrxx_mkII_DIR:PATH=' "$build/CMakeCache.txt" >/dev/null ||
        fail "$name MPLAPACK CMake discovery did not record package directory"
    pass "$name MPLAPACK CMake discovery"
}

cmake_discovery staged "$internal_prefix"
cmake_discovery independent "$independent_prefix"

consumer_source="$RUNTIME_DIR/consumer-source"
mkdir -p "$consumer_source"
cat > "$consumer_source/CMakeLists.txt" <<'EOF'
cmake_minimum_required(VERSION 3.16)
project(p1r_consumer LANGUAGES CXX)
find_package(gmpfrxx_mkII 1.1.0 EXACT CONFIG REQUIRED)
foreach(required_target
    gmpxx_mkII::gmpxx_mkII
    gmpxx_mkII::default_context_provider
    mpfrxx_mkII::mpfrxx_mkII
    mpcxx_mkII::mpcxx_mkII
    gmpfrxx_mkII::gmpfrxx_mkII)
  if(NOT TARGET ${required_target})
    message(FATAL_ERROR "missing exported target: ${required_target}")
  endif()
endforeach()
if(NOT gmpfrxx_mkII_VERSION VERSION_EQUAL "1.1.0")
  message(FATAL_ERROR "unexpected package version: ${gmpfrxx_mkII_VERSION}")
endif()
add_executable(p1r_consumer main.cpp)
target_link_libraries(p1r_consumer PRIVATE gmpfrxx_mkII::gmpfrxx_mkII)
EOF
cat > "$consumer_source/main.cpp" <<'EOF'
#include <gmpfrxx_mkII.h>

int main() {
    mpfrxx::mpfr_class mpfr_value(1);
    gmpxx::mpf_class mpf_value(2);
    return mpfr_value == 1 && mpf_value == 2 ? 0 : 1;
}
EOF

run_consumer() {
    local name=$1
    local package_prefix=$2
    local build="$RUNTIME_DIR/consumer-$name"
    cmake -S "$consumer_source" -B "$build" \
        -DCMAKE_PREFIX_PATH="$package_prefix;$deps_prefix/GMP;$deps_prefix/MPFR;$deps_prefix/MPC" \
        -DCMAKE_VERBOSE_MAKEFILE=ON
    cmake --build "$build" -j"$jobs" --verbose
    LD_LIBRARY_PATH="$package_prefix/lib:$LD_LIBRARY_PATH" "$build/p1r_consumer"
    if rg -n '(^|[ /])-lgmpxx([ ;]|$)|/libgmpxx\.(so|a|dylib)' "$build"; then
        fail "$name consumer unexpectedly links libgmpxx"
    fi
    pass "$name installed consumer compile, link, and run"
}

run_consumer installed "$final_prefix"
relocated_prefix="$RUNTIME_DIR/relocated-prefix"
mv "$final_prefix" "$relocated_prefix"
run_consumer relocated "$relocated_prefix"

make -C "$autotools_build" dist
mapfile -t dist_archives < <(find "$autotools_build" -maxdepth 1 -type f \
    -name 'mplapack-*.tar.*' -print | LC_ALL=C sort)
test "${#dist_archives[@]}" -gt 0 || fail "make dist produced no archive"
for dist_archive in "${dist_archives[@]}"; do
    tar -tf "$dist_archive" > "$RUNTIME_DIR/dist-files.txt"
    new_entry=$(grep -E "/external/gmpfrxx_mkII/download/${archive_name}$" \
        "$RUNTIME_DIR/dist-files.txt" || true)
    test -n "$new_entry" || fail "distribution archive lacks $archive_name"
    if grep -F "/external/gmpfrxx_mkII/download/gmpfrxx_mkII.${bootstrap_version}.tar.xz" \
        "$RUNTIME_DIR/dist-files.txt"; then
        fail "distribution archive contains bootstrap archive"
    fi
    tar -xOf "$dist_archive" "$new_entry" > "$RUNTIME_DIR/dist-vendored.tar.xz"
    test "$(sha256sum "$RUNTIME_DIR/dist-vendored.tar.xz" | awk '{print $1}')" = "$archive_sha" ||
        fail "distributed vendored archive digest mismatch"
    cmp --silent "$published_archive" "$RUNTIME_DIR/dist-vendored.tar.xz" ||
        fail "distributed vendored archive differs from published asset"
done
pass "make dist includes exact 1.1.0 asset and excludes 1.0.1"

python3 -m unittest -v "$source_tree/$migration_rel/tools/test_compare_baseline.py"
python3 "$source_tree/$migration_rel/tools/compare_baseline.py" \
    "$source_tree/$migration_rel/baseline.json" \
    "$source_tree/$migration_rel/baseline.json" \
    --rules "$source_tree/$migration_rel/baseline_rules.tsv"
pass "accepted P0 baseline comparator and immutable 874-result self-comparison"

if rg -n 'gmpfrxx_mkII' \
    "$source_tree/include" \
    "$source_tree/mpblas/reference" \
    "$source_tree/mplapack/reference" \
    "$source_tree/mpblas/optimized" \
    "$source_tree/mplapack/optimized" \
    "$source_tree/mpblas/test" \
    "$source_tree/mplapack/test" \
    "$source_tree/examples" \
    "$source_tree/benchmark"; then
    fail "MPLAPACK production, backend, test, example, or benchmark path includes gmpfrxx_mkII"
fi
test ! -d "$autotools_build/external/gmpfrxx_mkII/work/internal/_deps" ||
    fail "internal package created an auto-fetch dependency tree"
test ! -d "$autotools_build/external/gmpfrxx_mkII/work/install/_deps" ||
    fail "final package created an auto-fetch dependency tree"
test ! -d "$independent_build/_deps" || fail "independent package created an auto-fetch dependency tree"
pass "no backend binding change, production include, dependency auto-fetch, or network-required build"

test "$(git rev-parse HEAD)" = "$START_HEAD" || fail "HEAD changed during gate"
test "$(git rev-parse HEAD^{tree})" = "$START_TREE" || fail "HEAD tree changed during gate"
test "$(git write-tree)" = "$START_INDEX_TREE" || fail "index changed during gate"
test -z "$(git status --porcelain=v1 --untracked-files=all)" ||
    fail "primary worktree changed during gate"
git diff --check
pass "primary HEAD, tree, index, worktree, and whitespace are unchanged"

printf 'gate-P1R: PASS\n'
