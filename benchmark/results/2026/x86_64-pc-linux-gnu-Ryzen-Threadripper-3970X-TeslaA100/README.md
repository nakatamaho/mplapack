# 2026 measurements: AMD Ryzen Threadripper 3970X + NVIDIA A100 80GB PCIe

Host: AMD Ryzen Threadripper 3970X (32 cores, 64 threads), NVIDIA A100 80GB
PCIe (sm_80), CUDA 13.1, GCC 15.2, Ubuntu 26.04.  The host was shared with
other jobs; each directory's README states the load and how it was measured.
Same machine as ../../2022/x86_64-pc-linux-gnu-Ryzen-Threadripper-3970X-TeslaA100.

| Directory | Contents |
|---|---|
| `rgemm/` | October 2026: Rgemm/Cgemm GPU against CPU (dd/td/qd/mpfr) and the MPFR 512-bit size sweep up to n = 9997, through the public Rgemm (transfers included) |
| `factorization-goal08/` | October 2026: GPU solvers (getrs/potrs/gesv/posv), CPU / per-call GPU BLAS / resident entry |
| `factorization-goal08b/` | October-November 2026: below-gate CPU-library speed; factor table before and after the gate changes |
| `factorization-batch-09-10/` | November 2026: CPU LU small-n follow-up, CPU/GPU inverses, Cholesky profile, batch QA receipts |
| `gate-retune/` | November 2026: factor/solver gates retuned against the CPU library (measurements cut short; inverses not finished) |
| `kernel-bench/` | November 2026: the kernels of LU / Cholesky / inverses alone (GEMM, SYRK/HERK, TRSM, panels, laswp, trtri, lauum), GPU (no transfers) against CPU for dd/td/qd/MPFR, with graphs; part of the CPU values are loaded-host reference values |
| [kernel-bench-2/](kernel-bench-2/) | November 2026: new BLAS1, triangular solves and MPFR inverses; resident CUDA events versus CPU library (32-thread setting), 7,941 points / 76 graphs, reduced matrix grids and explicit MPFR 10s cutoffs; no dispatcher retuning |
| `gemm-sweep/` | November 2026: GPU kernel-only Rgemm NN/NT/TN/TT sweep, n = multiples of 7 and 31 and powers of two |
| `roofline/` | November 2026: % of the roofline (arithmetic peaks of each type and memory bandwidth, LAWN 41 counts) for every GPU kernel |
