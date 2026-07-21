# AGENTS.md — MPLAPACK repository guide for coding agents

## What this repository is

MPLAPACK: multiple-precision BLAS + LAPACK in C++. Backends: gmp, mpfr, qd,
dd, double, binary80, binary128. Two build systems in parallel: autotools
(configure.ac, Makefile.am tree, regenerated via ./gen_configure.sh) and
modern CMake (CMakeLists.txt + cmake/). Both must stay green.

## Hard rules

- Sources under `mpblas/reference/` and `mplapack/reference/` are
  FABLE-generated translations of Fortran LAPACK. NEVER edit routine bodies
  there. Fixes to numerics go through the fable pipeline, not hand edits.
- `include/mplapack_config.h` in the source tree may be a stale autotools
  artifact; the CMake build generates its own under the build dir and must
  keep winning via include order.
- dd backend requires -ffp-contract=off on every TU; removing it silently
  corrupts double-double arithmetic.
- Version numbers live in BOTH configure.ac and CMakeLists.txt; change both.
- libtool -version-info is ABI versioning, independent of package version.

## Ongoing work: unified library layout (branch topic/unified-libs)

The library taxonomy is being consolidated: one self-contained library per
flavor (libmplapack_<b>, libmplapack_<b>_opt, libmplapack_dd_opt_cuda,
libmplapack_binary128_opt_opencl), duplicate public symbols eliminated by a
layered basename-shadowing rule (accelerator > optimized > reference).
Normative spec: .agents/skills/unified-libs-goal/references/00-overview.md.
Execution units: references/01..05 in the same directory, driven by the
repo-local skill `unified-libs-goal` — one goal per session, in order
(invoke: "$unified-libs-goal 01"). Do not act on the goal files outside a
skill-driven session.

## Build quickstarts

autotools:
    ./gen_configure.sh
    ./configure --enable-dd --enable-double   # minimal fast config
    make -j$(nproc) && make check

CMake:
    cmake -B build -DMPLAPACK_ENABLE_GMP=OFF -DMPLAPACK_ENABLE_MPFR=OFF \
          -DMPLAPACK_ENABLE_QD=OFF -DMPLAPACK_ENABLE_BINARY128=OFF \
          -DMPLAPACK_BUILD_TESTS=ON
    cmake --build build -j && ctest --test-dir build

## Verification assets

- misc/check_unique_symbols.sh <lib...> — fails on duplicate global defined
  symbols inside one library. Run it on every archive you produce.
- Test suites live in mpblas/test and mplapack/test; they compare against
  reference results and are the ground truth for numerical correctness.
