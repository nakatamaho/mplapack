# GOAL 05: GPU Cholesky entry points (Rpotrf / Cpotrf, dd td qd)

Prerequisites: goals 02, 04. Read 00-overview.md (including "Known defects
found by goal 01"; the bug-fix track must be merged first) first.

## Scope

New public entry points Rpotrf_<b>_cuda, Cpotrf_<b>_cuda (b = dd, td, qd),
their dispatchers (Rpotrf.cpp / Cpotrf.cpp under the libraries' cuda/
directories, shadowing mplapack/optimized/common), headers through
misc/gen_accelerator_headers.py, both build systems.

## Tasks

1. Entry point: argument checks (return false on errors and n = 0), memory
   check, upload, blocked Cholesky on the device in reference order (nb from
   iMlaenv, Rpotrf2 on the diagonal block on the device, Rtrsm / Rsyrk /
   Rherk / Rgemm on the device), panel and update streams with events,
   download, info as the reference.
2. Fast mode `<P>_LAPACK_FAST=1`: allowed to change the block size and the
   recursion (and to use `<P>_CUDA_SLOPPY_*` arithmetic only if that
   variable is also set).
3. `<P>_LAPACK_CUDA=0`, `<P>_LAPACK_CUDA_MIN` (default from goal 01).
4. README of the cuda directory (environment, modes, limits).

## Acceptance criteria

- Rfactor_test CUDA variant: Rpotrf / Cpotrf through the entry points are
  bit-identical to the reference (both UPLO, info = j cases), and the entry
  point returns false (outputs untouched) for invalid arguments and n = 0.
- Host emulation variant passes in autotools make check and CMake.
- Timings (commit body), dd and qd, n = 1024, 2048, 4096 (dd also 8192):
  CPU default (goal 02), current _opt_cuda (GPU BLAS per call, goal 01),
  new entry point, with the time split upload / factor / download. The
  entry point must be faster than the per-call path at every n ≥
  `<P>_LAPACK_CUDA_MIN`.
- nm -D shows the new symbols; check_unique_symbols.sh passes.
- One commit: "Rpotrf/Cpotrf_<b>_cuda: device-resident Cholesky (factorization goal 05)".
