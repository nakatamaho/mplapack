#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/../../../.." && pwd)
cd "$root"
archive=external/gmpfrxx_mkII/download/gmpfrxx_mkII.1.0.1.tar.xz
expected=c0816b3538b6b77009f714bb391cebe11abb2fdb69e07aa3bb305ff822764afb
test "$(sha256sum "$archive" | awk '{print $1}')" = "$expected"
grep -F "$expected  gmpfrxx_mkII/download/gmpfrxx_mkII.1.0.1.tar.xz" external/distfiles.sha256 >/dev/null
cp "$archive" /tmp/gmpfrxx-p1-corrupt.tar.xz
printf X >> /tmp/gmpfrxx-p1-corrupt.tar.xz
test "$(sha256sum /tmp/gmpfrxx-p1-corrupt.tar.xz | awk '{print $1}')" != "$expected"
rm -f /tmp/gmpfrxx-p1-corrupt.tar.xz
test -f external/i/GMPFRXX_MKII/include/gmpfrxx_mkII.h
test -f external/i/GMPFRXX_MKII/lib/cmake/gmpfrxx_mkII/gmpfrxx_mkIIConfig.cmake
test -f /tmp/mplapack-p1-final/include/gmpfrxx_mkII.h
test -f /tmp/mplapack-p1-final/lib/cmake/gmpfrxx_mkII/gmpfrxx_mkIIConfig.cmake
test -z "$(find external/gmpfrxx_mkII/patches -type f -print)"
! rg -n 'gmpfrxx_mkII' include mpblas/reference mplapack/reference mpblas/optimized mplapack/optimized mpblas/test mplapack/test examples benchmark
rm -rf /tmp/mplapack-p1-stage-cmake /tmp/mplapack-p1-final-cmake
cmake -S . -B /tmp/mplapack-p1-stage-cmake -DMPLAPACK_ENABLE_GMPFRXX_MKII=ON -DMPLAPACK_ENABLE_GMP=OFF -DMPLAPACK_ENABLE_MPFR=OFF -DMPLAPACK_ENABLE_QD=OFF -DMPLAPACK_ENABLE_DD=OFF -DMPLAPACK_ENABLE_DOUBLE=ON -DMPLAPACK_ENABLE_BINARY80=OFF -DMPLAPACK_ENABLE_BINARY128=OFF -DCMAKE_PREFIX_PATH="$PWD/external/i/GMPFRXX_MKII;$PWD/external/i/GMP;$PWD/external/i/MPFR;$PWD/external/i/MPC" >/tmp/p1-cmake-stage.log
cmake -S . -B /tmp/mplapack-p1-final-cmake -DMPLAPACK_ENABLE_GMPFRXX_MKII=ON -DMPLAPACK_ENABLE_GMP=OFF -DMPLAPACK_ENABLE_MPFR=OFF -DMPLAPACK_ENABLE_QD=OFF -DMPLAPACK_ENABLE_DD=OFF -DMPLAPACK_ENABLE_DOUBLE=ON -DMPLAPACK_ENABLE_BINARY80=OFF -DMPLAPACK_ENABLE_BINARY128=OFF -DCMAKE_PREFIX_PATH="/tmp/mplapack-p1-final;$PWD/external/i/GMP;$PWD/external/i/MPFR;$PWD/external/i/MPC" >/tmp/p1-cmake-final.log
test -s docs/migration/gmpfrxx_mkII/baseline.json
echo 'P1 gate: PASS'
