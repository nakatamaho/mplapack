# GOAL 01: Tests and benchmarks for getrf / potrf (no new algorithm)

Prerequisites: none. Read 00-overview.md first.

## Scope

Create the verification and measurement assets every later goal uses. No
change to any library routine.

## Tasks

1. Test program mplapack/optimized/common/test/Rfactor_test.cpp (one source,
   built per backend like mpblas/optimized/common/test/*): for Rgetrf,
   Cgetrf, Rpotrf, Cpotrf (both UPLO) compare the routine of the library it
   is linked with against the reference (`<R>_reference_cpu`, or the
   reference library) with memcmp of A, ipiv, info on the matrix set of
   00-overview.md "Verification gates", and compute the backward errors.
   Environment switch for the mode under test: MPLAPACK_LAPACK_FAST=1 turns
   the memcmp checks into backward-error checks. Exit 0 / 1.
2. Register it: the optimized libraries and their tests are built from
   mpblas/optimized/<b>/Makefile.am (mplapack/optimized/<b> only holds a
   placeholder), which misc/gen_backend_makefiles.py generates from
   misc/backend_makefiles/mpblas_optimized.am.in — add the test there, never
   by hand-editing the generated Makefile.am; CMake add_test per backend.
   A CUDA variant linked against libmplapack_<b>_opt_cuda for dd, td, qd,
   mpfr (CMake add_test with SKIP_RETURN_CODE 77; it exercises the GPU BLAS
   the factorizations call today).
3. Benchmarks benchmark/Rgetrf.cpp and benchmark/Rpotrf.cpp exist already
   for the CPU flavours: check what they measure; extend or add programs so
   that the measurement of 00-overview.md "Benchmarks" is available for
   _opt and _opt_cuda builds (size list, best of 3 up to 2000, LAPACK flop
   counts, CSV output), without changing existing benchmark outputs that
   benchmark/results relies on.
4. Measure and record in the commit body (no results files committed):
   dd, qd, mpfr 512 bits, n = 512, 1024, 2048 (and 4096 for dd): reference,
   _opt (look-ahead), _opt_cuda (GPU BLAS per call), MPLAPACK_LAPACK_FAST=1.
   From these choose the default of `<P>_LAPACK_CUDA_MIN` for later goals
   and write it in this goal's commit body.

## Acceptance criteria

- Rfactor_test passes for every enabled backend in both build systems, and
  its CUDA variant passes on a GPU (or exits 77 without one).
- Deliberately breaking one element of a factor in a scratch build makes the
  test fail (state that you checked this; do not commit the breakage).
- Benchmarks build and print the table; numbers in the commit body.
- misc/check_unique_symbols.sh passes on the changed archives.
- One commit: "LAPACK factorizations: tests and benchmarks (factorization goal 01)".
