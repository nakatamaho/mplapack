# Kernel benchmark: CPU against GPU, per precision class

`Rkernel_bench` times, alone, the kernels from which the LU / Cholesky
factorizations, the solvers and the inverses are built (GEMM, SYRK/HERK,
TRSM, the getrf2 / potrf2 panels, laswp, trtri, lauum), on the CPU and on
the GPU, for dd, td, qd and MPFR.  The results of the A100 host are in
`results/2026/x86_64-pc-linux-gnu-Ryzen-Threadripper-3970X-TeslaA100/kernel-bench/`.

## Programs

One source, `Rkernel_bench.cpp`, three kinds of programs:

| Program | What is timed |
|---|---|
| `Rkernel_bench.<b>_opt` (b = dd, td, qd, mpfr) | the routines of `libmplapack_<b>_opt` (Rgemm, Rsyrk, Rtrsm, Rgetrf2, Rpotrf2, Rlaswp, Rtrtri, Rlauum, Cgemm, Cherk, Ctrsm; OpenMP), wall clock of the call |
| `Rkernel_bench.<b>_cuda_kernel` (b = dd, td, qd; source `Rkernel_bench_<b>_cuda.cu`) | the device kernels of the CUDA libraries, compiled into the program from their internal headers (`mpblas/optimized/common/R{gemm,blas,lapack}_cuda_common.cuh`, `Cblas_cuda_common.cuh`, `mplapack/optimized/common/R{getrf,potrf,inverse}_gpu_common.cuh`, hidden visibility) with the flags of the libraries (`-fmad=false`, no contraction), on matrices already on the device.  CUDA events: the start event is synchronized, the stop event is recorded after the kernels' own synchronization (`gemm_device`, `blas_device`, `job_device`, the stream contexts' `finish()`) |
| `Rkernel_bench.mpfr_cuda_kernel` | the entry points `<R>_mpfr_cuda` of `libmplapack_mpfr_opt_cuda` (Rgemm, Rsyrk, Rtrsm, Cgemm, Cherk, Ctrsm, Rgetrf, Rpotrf).  Their kernels are in anonymous namespaces of `mpblas/optimized/mpfr/cuda/*.cu` and cannot be called on device data from outside the library, so the program records the kernels of the call with CUPTI (activity API) and reports the span from the start of the first kernel to the end of the last one.  Packing into `cu_freal` / `cu_mpfr` and the transfers are not included; status words copied between the kernels of a factorization are |

The GPU programs call the same kernels as the libraries: `gemm_device`
(tiled, the default algorithm of Rgemm), `blas_device` (Rsyrk tiled, Rtrsm
lanes), `job_Cgemm` / `job_Cherk` / `job_Ctrsm`, the panel of the device LU
(`lu_context::recursive`), the diagonal block of the device Cholesky
(`chol_recursive`), `laswp`, and the inverse engine (`context::tri`,
`context::lauum`) with the block sizes of `iMlaenv`.  The entry points never
fall back to the CPU; a call they reject is reported as a failed point.

Usage: `Rkernel_bench.X KERNEL FORM SIZES [--label L] [--reps R] [--long S]
[--warmup W] [--dump DIR] [--dump-sizes N,N]`; output lines
`label,kernel,form,n,method,seconds,gflops,runs`.  The MPFR programs need
`MPFRXX_DEFAULT_PRECISION_BITS` (the precision of every thread, also the
OpenMP workers).

## Kernels and shapes

| kernel | form | call | FLOP count |
|---|---|---|---|
| gemm | NN, NT, TN | `Rgemm(ta,tb,n,n,n,-1,A,B,1,C)` | 2n³ |
| gemm_update | NN | `Rgemm("N","N",n,n,64,-1,A,B,1,C)` (trailing update) | 2·64·n² |
| syrk | LN, UT (k = n); LN_k64, UT_k64 | `Rsyrk(uplo,trans,n,k,-1,A,1,C)` | k·n(n+1) |
| trsm | LLNU, RLTN, LUTN, LUNN (n × n) | `Rtrsm(side,uplo,trans,diag,n,n,1,A,B)` | n³ |
| trsm | `_p64`: triangle 64 × 64, B 64 × n (Left) or n × 64 (Right) | as in the factorizations | n·64² |
| getrf2 | P64 | `Rgetrf2(n,64,A,ipiv)`: the LU panel (GPU: the panel of the device LU; MPFR: `Rgetrf_mpfr_cuda(n,64)`, nb = 64 so it is the recursive panel) | LAWN 41 getrf(m = n, 64) |
| potrf2 | L, n = 64 and 256 | `Rpotrf2("L",n,A)` (MPFR GPU: `Rpotrf_mpfr_cuda`, which at n = 256 runs the blocked algorithm with nb = 64) | n³/3 + n²/2 + n/6 |
| laswp | R64 | `Rlaswp(n,A,n,1,64,ipiv,1)`, ipiv(i) ≥ i random | none (seconds; graphs: 64·n swaps / s) |
| trtri | UN | `Rtrtri("U","N",n,A)` | n³/3 + 2n/3 |
| lauum | L | `Rlauum("L",n,A)` | n³/3 + n²/2 + n/6 |
| cgemm | NN, NC (CN) | `Cgemm(ta,tb,n,n,n,-1,A,B,1,C)` | 8n³ |
| cherk | LN, UC (k = n), `_k64` | `Cherk(uplo,trans,n,k,-1,A,1,C)` | 4·k·n(n+1) |
| ctrsm | LLNU, RLCN, LUCN, LUNN, `_p64` | `Ctrsm(...)` | 4 × TRSM |

Data are deterministic and use every bit of the type (dd/td/qd: each limb
random; MPFR: `mpfr_urandomb` at the working precision), uniform in [-1, 1).
Triangular factors have their off-diagonal entries divided by the order and
a diagonal in [1, 1.5) (well conditioned, unit and non-unit); the Cholesky
blocks are symmetric with the order on the diagonal.  The largest size is
generated once; smaller sizes use the leading blocks (CPU and MPFR GPU:
leading dimension of the largest size; dd/td/qd GPU: packed, ld = rows).
The in/out operand is restored before every run (outside the timing).

## Driver: kernel_bench.py

    kernel_bench.py probe --bindir B --out DIR      # two sizes per series, the time model
    kernel_bench.py plan  --out DIR [--max-point S] [--max-point-fine S] [--coarse]
    kernel_bench.py run   --bindir B --out DIR [--hours H] [--only gpu|cpu]
    kernel_bench.py check --bindir B --out DIR      # bit comparison at n = 64 and 200
    kernel_bench.py merge --out DIR [--cpu-sources quiet=raw/cpu-quiet,loaded=raw/cpu-loaded]
                                                    # -> DIR/kernels.csv (CPU: first source with the series; column cpu_host)
    kernel_bench.py table --out DIR --types dd,...  # markdown tables

Types: dd, td, qd, mpfr256, mpfr512, mpfr1024 (cu_freal), mpfr200rt
(cu_mpfr, run-time precision).  GPU sizes: n = 32..1024 step 32 and
1152..4096 step 128 (dd/td/qd GEMM and SYRK k = n up to 8192 step 256; MPFR
up to 4096, 1024 bits and 200 bits up to 2048), every power of two; CPU:
the powers of two 128..4096 (MPFR 2048).  `plan --coarse` keeps this grid
for the FINE series of `kernel_bench.py` and uses steps 64 / 256 (MPFR 512)
/ 512 for the others; `--max-point S` drops the GPU sizes whose predicted
run exceeds S seconds (the time model is fitted on the probe).  GPU: one
warm-up run and the best of 3 (one run when the warm-up exceeds 10 s); CPU:
the best of up to 3 (one run when a run exceeds 10 s).  `run` executes the
GPU series and the CPU series in two concurrent streams (CPU programs with
OMP_NUM_THREADS threads placed by the OS; `--pin` binds them to cpus
0..OMP_NUM_THREADS-1, which on a shared host makes the barrier-heavy panels
wait for descheduled threads; GPU programs with one OpenMP thread), records `nvidia-smi` before and after every GPU series
(`gpu-processes.log`; a series during which another process used the GPU is
listed in `gpu-busy.txt` and flagged in `kernels.csv`), re-reads the plan
before every series and starts no series after `--hours`.

`check` runs every series at n = 64 and 200 on both sides with `--dump` and
compares the results bit for bit (`kernel_bench_compare.py`: elements, info
and ipiv; a difference only in the sign of a zero is counted separately).

## Graphs: plot_kernel_bench.py

    PYTHONPATH=<matplotlib> python3 plot_kernel_bench.py DIR/kernels.csv DIR --gpu-name "NVIDIA A100 80GB PCIe"

draws, with `plot_sweep.draw()` (the code of
`results/2026/.../rgemm/mpfr512_sweep10k.png`), one graph per type and
series (`DIR/<type>/<kernel>_<form>.png`: GPU line, powers of two marked
and annotated, CPU powers of two dashed), the GPU/CPU speed ratio per type
(`ratio_<type>.png`) and the GFLOPS at n = 2048 (MPFR 1024) of every type
and kernel (`summary.png`).  CPU series of a source named `loaded` are
labelled as loaded-host reference values (legend, subtitle, ▲ in the ratio
graphs, hatched bars in the summary).

## Building

autotools: the programs are in `Makefile.{dd,td,qd,mpfr}_cuda.am` (built
with the CUDA libraries; `make -C benchmark Rkernel_bench.dd_cuda_kernel
...`).  CMake: `-DMPLAPACK_BUILD_BENCHMARKS=ON` with the CUDA libraries
enabled builds the same programs in `<build>/benchmark`.  The MPFR GPU
program links `libcupti` (CUDA toolkit).
