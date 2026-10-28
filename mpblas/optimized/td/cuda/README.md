# Device-resident Cholesky (td)

`Rpotrf_td_cuda` / `Cpotrf_td_cuda` are declared in `mpblas_td_cuda.h`.
See [the dd implementation notes](../../dd/cuda/README.md) for the contract,
reference ordering, memory limits and timing definitions.

The variables are `MPLAPACK_TD_LAPACK_CUDA` (0 disables),
`MPLAPACK_TD_LAPACK_CUDA_MIN` (default 2048, QD proxy; TD not measured),
`MPLAPACK_TD_LAPACK_FAST` (currently the same reference-order path), and
`MPLAPACK_TD_CUDA_VERBOSE` (1 prints upload/factor/download seconds).
They are read on every call. No implicit sloppy arithmetic is used.
