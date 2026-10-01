#!/usr/bin/env bash
# Regenerate the release sources with the Fable pipeline from the same Git ref.
set -euo pipefail

baseline="${1:?source directory is required}"
jobs="${2:?job count is required}"
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
report_dir="${MPLAPACK_TEST_RESULTS_BASE:-/results}/fable-reproduction"
mkdir -p "$report_dir"
qa_root="$(mktemp -d /work/fable-reproduction.XXXXXX)"
source_dir="$qa_root/source"
source_ref="${MPLAPACK_REF:?}"
if [ -e "$baseline/.git" ]; then
    source_ref="$(git -C "$baseline" rev-parse HEAD)"
elif [[ ! "$source_ref" =~ ^[0-9a-fA-F]{40}$ ]]; then
    echo "ERROR: tarball reproduction requires a full Git SHA in MPLAPACK_REF" >&2
    exit 1
fi

echo "=== Fable source reproduction regression ==="
/usr/local/bin/checkout-source.sh "${MPLAPACK_REPO:?}" "$source_ref" "$source_dir"
git -C "$source_dir" rev-parse HEAD | tee "$report_dir/git-sha.txt"
for tool in python python3 parallel ctags clang-format; do
    command -v "$tool" >/dev/null || {
        echo "ERROR: Fable regression requires $tool" >&2
        exit 1
    }
done
ctags --version | head -n 1
clang-format --version
export JOBS="$jobs" LC_ALL=C

if ! (
    cd "$source_dir"
    autoreconf --force --install || exit
    ./configure --disable-dependency-tracking --enable-gmp=no --enable-mpfr=no --enable-qd=no \
        --enable-dd=no --enable-double=yes --enable-binary80=no \
        --enable-binary128=no --disable-test --disable-benchmark --with-openblas=no || exit
    export LAPACK_VERSION
    LAPACK_VERSION="$(sed -n 's/^LAPACK_VER_MAJOR="\([^"]*\)"/\1/p' configure.ac).$(sed -n 's/^LAPACK_VER_MINOR="\([^"]*\)"/\1/p' configure.ac).$(sed -n 's/^LAPACK_VER_PATCH="\([^"]*\)"/\1/p' configure.ac)"
    echo "LAPACK_VERSION=$LAPACK_VERSION"
    bash fable/go.sh || exit
    ROOT="$source_dir" bash fable/go_testing.sh || exit
) > "$report_dir/generation.log" 2>&1; then
    tail -n 60 "$report_dir/generation.log"
    echo "FAIL: Fable regeneration; see $report_dir/generation.log" >&2
    exit 1
fi

python3 "$script_dir/compare-fable-outputs.py" "$baseline" "$source_dir" \
    > "$report_dir/comparison.log" 2>&1 || {
    tail -n 60 "$report_dir/comparison.log"
    echo "FAIL: Fable source reproduction; see $report_dir" >&2
    exit 1
}
cat "$report_dir/comparison.log"
echo "PASS: fable/go.sh and fable/go_testing.sh reproduce the release sources"
