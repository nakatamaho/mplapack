# Device-resident Cholesky (td)

Resident R/C TRTRI/LAUUM/POTRI/GETRI (goal 09):
[API and CPU-library-speed fallback](../../../../mplapack/optimized/common/Rinverse_cuda.md).
Default is opt-in; `MPLAPACK_TD_LAPACK_INVERSE_CUDA_MIN=0` forces eligible
small calls; inverse measurements were cut short.
The [warm CPU-library retune](../../../../benchmark/results/2026/x86_64-pc-linux-gnu-Ryzen-Threadripper-3970X-TeslaA100/gate-retune/README.md)
lowers Rpotrf to 512 (both triangles win); 256 is marginal for Upper.
Cpotrf is opt-in: Upper loses through 4096 and its gate is shared with Lower.
Small-n CPU LU load anomalies are excluded, so LU remains 1024.

Resident getrs/potrs/gesv/posv (real and complex):
[API, controls, verification and timings](../../../../mplapack/optimized/common/Rsolver_cuda.md).

`Rpotrf_td_cuda` / `Cpotrf_td_cuda` are declared in `mpblas_td_cuda.h`.
See [the dd implementation notes](../../dd/cuda/README.md) for the contract,
reference ordering, memory limits and timing definitions.

The variables are `MPLAPACK_TD_LAPACK_CUDA` (0 disables),
`MPLAPACK_TD_LAPACK_CUDA_MIN` (Rpotrf default 512, Cpotrf opt-in),
`MPLAPACK_TD_LAPACK_FAST` (currently the same reference-order path), and
`MPLAPACK_TD_CUDA_VERBOSE` (1 prints upload/factor/download seconds).
They are read on every call. No implicit sloppy arithmetic is used.

`Rgetrf_td_cuda` / `Cgetrf_td_cuda` in the same header implement rectangular,
device-resident LU under the same variables (minimum applies to min(m,n),
LU default 1024). Three alternating rounds at 1024/2048 measured faster
resident LU for real and complex. Historical Cholesky crossover claims
against the old per-call GPU path are superseded by the CPU-library retune.
See the dd LU notes for reference ordering, pivot semantics, two-stream
scheduling and transactional A/ipiv/info publication. The dispatchers use
the optimized CPU LU when an entry returns false.
