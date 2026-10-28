# Device-resident Cholesky (qd)

`Rpotrf_qd_cuda` / `Cpotrf_qd_cuda` are declared in `mpblas_qd_cuda.h`.
See [the dd implementation notes](../../dd/cuda/README.md) for the contract,
reference ordering, memory limits and timing definitions.

The variables are `MPLAPACK_QD_LAPACK_CUDA` (0 disables),
`MPLAPACK_QD_LAPACK_CUDA_MIN` (default 2048),
`MPLAPACK_QD_LAPACK_FAST` (currently the same reference-order path), and
`MPLAPACK_QD_CUDA_VERBOSE` (1 prints upload/factor/download seconds).
They are read on every call. No implicit sloppy arithmetic is used.
The provisional per-BLAS-call threshold 512 was retuned to 2048: interleaved
three-round A100 measurements found the resident path slower at 512/1024
and faster at 2048/4096 (full-call times, including allocation and packing).
