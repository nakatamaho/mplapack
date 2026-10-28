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
