# Device-resident Cholesky (qd)

Resident R/C TRTRI/LAUUM/POTRI/GETRI (goal 09):
[API and CPU-library-speed fallback](../../../../mplapack/optimized/common/Rinverse_cuda.md).
Default is opt-in; `MPLAPACK_QD_LAPACK_INVERSE_CUDA_MIN=0` forces eligible
small calls; inverse measurements were cut short.
The [warm CPU-library retune](../../../../benchmark/results/2026/x86_64-pc-linux-gnu-Ryzen-Threadripper-3970X-TeslaA100/gate-retune/README.md)
lowers Rpotrf and R/Cgetrf to 512. Cpotrf becomes opt-in because Upper
loses through 4096; the existing gate is shared with Lower.

Resident getrs/potrs/gesv/posv (real and complex):
[API, controls, verification and timings](../../../../mplapack/optimized/common/Rsolver_cuda.md).

`Rpotrf_qd_cuda` / `Cpotrf_qd_cuda` are declared in `mpblas_qd_cuda.h`.
See [the dd implementation notes](../../dd/cuda/README.md) for the contract,
reference ordering, memory limits and timing definitions.

The variables are `MPLAPACK_QD_LAPACK_CUDA` (0 disables),
`MPLAPACK_QD_LAPACK_CUDA_MIN` (Rpotrf default 512, Cpotrf opt-in),
`MPLAPACK_QD_LAPACK_FAST` (currently the same reference-order path), and
`MPLAPACK_QD_CUDA_VERBOSE` (1 prints upload/factor/download seconds).
They are read on every call. No implicit sloppy arithmetic is used.
Warm comparisons now use the CPU library, not the historical per-call GPU
baseline: Rpotrf at 512 wins 1.62x/2.06x for L/U; Lower at 256 is only
parity. Cpotrf Upper at 4096 is 1.40x slower; its Lower win cannot safely
set the shared U/L constant. Explicit minimum overrides remain available.

`Rgetrf_qd_cuda` / `Cgetrf_qd_cuda` in the same header implement rectangular,
device-resident LU under the same variables (minimum applies to min(m,n),
LU default 512). Three alternating warm rounds at 512 measured speedups
of 1.48x real and 1.60x complex against the CPU library. The Cholesky defaults above are separate.
See the dd LU notes for reference ordering, pivot semantics, two-stream
scheduling and transactional A/ipiv/info publication. The dispatchers use
the optimized CPU LU when an entry returns false.
