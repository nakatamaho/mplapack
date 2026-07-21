# GOAL 04: CMake parity audit for the unified layout

Prerequisite: goals 01-03 merged on topic/unified-libs. Read 00-overview.md.

## Scope

Goals 01-03 changed CMake incrementally. This goal is a parity audit: every
library, flag, test, and install artifact the autotools build produces under
the new layout must have a CMake equivalent, and the exported package must be
consumable.

## Audit checklist (produce a table in the commit message body)

For each of: mplapack_<b> (7 backends), mplapack_<b>_opt (optimized-capable
backends), mplapack_dd_opt_cuda, mplapack_binary128_opt_opencl —

1. Target exists, alias exists, installed, exported in mplapackTargets.cmake.
2. PUBLIC compile definition (backend macro) and include dirs match.
3. Link deps match autotools (GMP/MPFR/MPC/QD/quadmath/OpenMP/CUDA/OpenCL).
4. dd: -ffp-contract=off present on all C++ TUs including _opt and _opt_cuda.
5. VERSION/SOVERSION consistent across all targets.
6. pkg-config file(s): decide and implement ONE convention — either one
   mplapack.pc with Libs listing nothing backend-specific plus per-flavor
   mplapack_<b>[_opt[_cuda|_opencl]].pc files, or per-flavor only. Autotools
   and CMake must emit identical .pc contents (diff them in a test).

## Consumption tests (new)

Add `cmake/tests/consume/` mini-project (a CMakeLists + one TU calling Rgemm
and one LAPACK routine, e.g. Rpotrf) and CTest fixtures that:
1. Build against the build tree via `find_package(mplapack CONFIG)` pointing
   at the export dir, linking `mplapack::mplapack_dd`.
2. Repeat for `mplapack::mplapack_dd_opt` when MPLAPACK_ENABLE_OPT=ON.
3. FetchContent smoke test is optional; skip if it doubles configure time.

## Test-suite wiring

- Ensure the ctest suite runs the same test binaries `make check` runs for at
  least dd and double, for both reference and _opt flavors, and that
  `RunStdinTest.cmake` covers the _opt variants.
- CI: if .github/workflows exists, extend the matrix with a
  `-DMPLAPACK_ENABLE_OPT=ON` job; CUDA/OpenCL jobs configure-only.

## Acceptance criteria

- `cmake --build build --target install` into a scratch prefix, then the
  consume project builds against that prefix for dd (reference and _opt).
- `ctest` green with dd+double+opt enabled.
- No CMake code path still references removed targets (grep `mpblas_`).
- One commit: "CMake parity for unified library layout (goal 04)" with the
  audit table in the body.
