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

## Device-resident Cholesky (goal 05)

`Rpotrf_gpu_test.cpp` exercises the six dd/td/qd entry points directly;
`Rpotrf_gpu_device_test.cu` is its nvcc wrapper. CMake runs
`{dd,td,qd}_cuda_Rpotrf_{host,device}` (device absence exits 77).
Autotools runs `Rpotrf_gpu_host_test_<b>` in each optimized backend's
`make check`. Each variant performs 398 checks: memcmp of the entire padded
matrix and info against the reference, both triangles, block boundaries,
odd n through 599, non-positive pivots, imaginary complex diagonals,
invalid/disabled/quick/overflow rejection without publication, and normalized
backward error <=30. Backend LAPACK_FAST=1 currently uses the same arithmetic
and is tested too. Device oracles disable GPU BLAS, not LAPACK entry points.
The test fixes OpenMP threads to one to avoid tiny-oracle parallel overhead.

## Device-resident LU (goal 06)

`Rgetrf_gpu_test.cpp` and its nvcc wrapper exercise all six LU entries
without a CPU fallback. CMake registers `{dd,td,qd}_cuda_Rgetrf_{host,device}`;
autotools `make check` runs `Rgetrf_gpu_host_test_<b>`. Each variant performs
890 host / 892 device checks: whole padded A, ipiv (including sentinels) and
info memcmp against
the CPU-only reference oracle, square/tall/wide matrices, block boundaries,
odd sizes through 599, zero/duplicate columns, forced swaps and tied maxima.
Tiny/subnormal pivots cover reciprocal scaling and direct division; singular
factorizations complete. A forced nb=1 recursive real/complex check also
verifies bounded scratch for wide matrices. Invalid/quick/disabled/threshold/overflow calls
preserve all outputs. Both modes check normalized backward error <=30 on
small matrices. Device tests exit 77 without a device. Set OMP_NUM_THREADS=1
for the host tests (as CMake does) to avoid tiny-oracle thread overhead.

For fault injection, copy the source outside the tree, add `g[0] += T(1);`
after the tested factorization, build the copy with the same flags/library,
and verify exit 1 in both modes. Do not modify the library or commit the fault.
