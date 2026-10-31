# Device-resident Cholesky (td)

Resident getrs/potrs/gesv/posv (real and complex):
[API, controls, verification and timings](../../../../mplapack/optimized/common/Rsolver_cuda.md).

`Rpotrf_td_cuda` / `Cpotrf_td_cuda` are declared in `mpblas_td_cuda.h`.
See [the dd implementation notes](../../dd/cuda/README.md) for the contract,
reference ordering, memory limits and timing definitions.

The variables are `MPLAPACK_TD_LAPACK_CUDA` (0 disables),
`MPLAPACK_TD_LAPACK_CUDA_MIN` (Cholesky default 2048),
`MPLAPACK_TD_LAPACK_FAST` (currently the same reference-order path), and
`MPLAPACK_TD_CUDA_VERBOSE` (1 prints upload/factor/download seconds).
They are read on every call. No implicit sloppy arithmetic is used.

`Rgetrf_td_cuda` / `Cgetrf_td_cuda` in the same header implement rectangular,
device-resident LU under the same variables (minimum applies to min(m,n),
LU default 1024). Three alternating rounds at 1024/2048 measured faster
resident LU for real and complex. Real/complex Cholesky, both UPLOs, also
won at 2048, so its former QD-proxy default is now backed by TD measurements.
See the dd LU notes for reference ordering, pivot semantics, two-stream
scheduling and transactional A/ipiv/info publication. The dispatchers use
the optimized CPU LU when an entry returns false.
