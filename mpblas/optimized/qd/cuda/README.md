# Device-resident Cholesky (qd)

`Rpotrf_qd_cuda` / `Cpotrf_qd_cuda` are declared in `mpblas_qd_cuda.h`.
See [the dd implementation notes](../../dd/cuda/README.md) for the contract,
reference ordering, memory limits and timing definitions.

The variables are `MPLAPACK_QD_LAPACK_CUDA` (0 disables),
`MPLAPACK_QD_LAPACK_CUDA_MIN` (Rpotrf default 2048, Cpotrf default 4096),
`MPLAPACK_QD_LAPACK_FAST` (currently the same reference-order path), and
`MPLAPACK_QD_CUDA_VERBOSE` (1 prints upload/factor/download seconds).
They are read on every call. No implicit sloppy arithmetic is used.
The provisional per-BLAS-call threshold 512 was retuned to 2048: interleaved
three-round A100 measurements found the resident path slower at 512/1024
and faster at 2048/4096 (full-call times, including allocation and packing).
Goal 06 measured complex Cholesky with both UPLOs: Upper lost by about 2%
at 2048, but both triangles won at 4096 (three alternating rounds). Cpotrf
therefore uses 4096 for both triangles; Rpotrf retains 2048.

`Rgetrf_qd_cuda` / `Cgetrf_qd_cuda` in the same header implement rectangular,
device-resident LU under the same variables (minimum applies to min(m,n),
LU default 1024). Three alternating rounds at 1024/2048 measured faster
resident LU for real and complex. The Cholesky defaults above are separate.
See the dd LU notes for reference ordering, pivot semantics, two-stream
scheduling and transactional A/ipiv/info publication. The dispatchers use
the optimized CPU LU when an entry returns false.
