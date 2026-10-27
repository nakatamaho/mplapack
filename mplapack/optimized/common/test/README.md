# Factorization verification

`Rfactor_test.cpp` is linked once per optimized backend, and with the CUDA
library where available (CMake). Run with `MPLAPACK_LAPACK_FAST=0` for exact
comparison with the existing `*_reference_cpu` blocked routines, or `=1` for
backward-error checks. Both modes check info and non-tied pivots. CUDA builds
return 77 without a device. Tests never enable fast mode themselves. CPU tests
set MPLAPACK_LAPACK_LOOKAHEAD=2 unless explicitly overridden, exercising tasks
at small sizes and for Cholesky (production defaults would bypass them).
CUDA tests leave that switch at its production default: LU tasks require
min(m,n) >= 12*nb and Cholesky tasks are disabled. They also check the final
CUDA runtime status so that a poisoned context followed by silent CPU fallback
cannot make a device test appear to pass.

Known pre-existing failure: explicitly setting MPLAPACK_LAPACK_LOOKAHEAD=2
with a CUDA-linked test produces identity/backward-error failures (and an
incorrect Cholesky info) on dd/td/qd. Goal 01 does not change library routines;
this is retained as an explicit reproducer, not an expected-pass mode:

    OMP_NUM_THREADS=4 MPLAPACK_LAPACK_LOOKAHEAD=2 ./dd_opt_cuda_Rfactor_test

Exact diagnostics vary with scheduling. Check this before enabling GPU task
look-ahead in later goals.

The fixed seed covers real and complex LU, tall/wide matrices, padded leading
dimensions, block boundaries, odd sizes through 599, zero/dependent columns,
forced swaps, tied maxima, and Cholesky with both triangles and a known failing
pivot. Residual products are evaluated at backend precision, not in double.
The bound is 30 in units of n*norm(A)*Rlamch("E"). Complex column sums use
|Re|+|Im| (the BLAS abs1 convention). CMake gives each test a 30-minute budget.

Exact comparisons use memcmp for fixed-size scalar representations (and for
pivots). MPFR/GMP object bytes contain heap pointers, so exact values are
compared instead (MPFR includes precision and signed zero). The six padding
bytes of x87 binary80 are likewise not compared. Complex components are
checked separately, excluding any object padding. The existing reference CPU
symbols use the BLAS of the linked library, as in Rlapack_omp_test; these tests
are not an independent test of BLAS kernels.

For fault injection, copy the source outside the tree, add `g[0] += T(1);`
after the tested factorization, build the copy with the same flags/library,
and verify exit 1 in both modes. Do not modify the library or commit the fault.
