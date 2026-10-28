# Device-resident Cholesky (dd)

`Rpotrf_dd_cuda` / `Cpotrf_dd_cuda` (in `mpblas_dd_cuda.h`) return true
when the GPU computed A and info. Invalid arguments, n=0, ineligible calls,
allocation/transfer/kernel failures return false without publishing any output.
The public Rpotrf/Cpotrf dispatchers then run the existing CPU path.

Environment (read on every call):
- `MPLAPACK_DD_LAPACK_CUDA=0`: disable these entry points.
- `MPLAPACK_DD_LAPACK_CUDA_MIN=N`: minimum n (default 4096).
- `MPLAPACK_DD_LAPACK_FAST=1`: currently uses the same reference-order path;
  no alternate rounding or sloppy arithmetic is selected implicitly.
- `MPLAPACK_DD_CUDA_VERBOSE=1`: upload/factor/download seconds on stderr;
  packing/unpacking is included, allocation/free is excluded.

One packed device matrix, O(n*nb) device scratch, and a transactional host
staging matrix are used. Opposite triangle and lda padding are not written.
The default uses iMlaenv's block size and the reference potrf2 recursion.
Only scalar status words return to the host while factoring. A panel stream
and an update stream are connected by ready/done events: each block row's
next-panel segment is solved first, then the far segment overlaps the next
diagonal factorization. Both streams are joined before publishing results.

# Device-resident LU (dd)

`Rgetrf_dd_cuda` / `Cgetrf_dd_cuda`, declared in `mpblas_dd_cuda.h`, use
one upload and download of the rectangular matrix and download the pivots.
They return false without changing A, ipiv or info on invalid arguments,
quick returns, disabled/below-threshold calls, or runtime errors. The CUDA
Rgetrf/Cgetrf dispatchers then use the goal-03 CPU LU path.

The environment variables above also control LU; the threshold applies to
min(m,n), with a LU-specific default of 2048 (Cholesky remains 4096).
Interleaved three-round real/complex measurements at 2048 and 4096 found
resident LU faster than per-call transfers at both sizes. These are
conservative measured eligibility defaults, not an extrapolated crossover.
FAST currently retains reference arithmetic. The reference block
size, recursive panels, first-index pivot ties, and singular-info semantics
are preserved. Singular pivots do not stop factorization. Panels, swaps,
scaling, TRSM and GEMM run on the device. Next-panel updates precede far
updates on the second stream; left swaps wait for the preceding far update.
The matrix plus O(max(m,n)*nb) workspace is checked against available device
memory. Padding is preserved; all host outputs are committed only after both
streams and all downloads succeed. Verbose times include packing/pivot
publication, but not allocation/initialization/free.
