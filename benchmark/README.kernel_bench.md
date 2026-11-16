# Kernel benchmark: CPU against GPU, per precision class

## Goal11 extension (new routines only)

Upper-case routine names select `Rkernel_new.h` from `Rkernel_bench.cpp`.
Lower-case legacy names and their original timing paths are unchanged.

```
Rkernel_bench.dd_opt Rdot fast 1024,2048
Rkernel_bench.dd_cuda_kernel Rtrsv fast 1024,2048
MPFRXX_DEFAULT_PRECISION_BITS=512 Rkernel_bench.mpfr_cuda_resident Cgetri fast 128
```

The new path covers all11 BLAS1 entries and all6 triangular solves for
DD/TD/QD/MPFR, and all8 MPFR inverses. BF inverses are existing work and are
not included. `fast`/`default` select the CPU reduction switch; GPU reductions
accept only `fast`. Non-reductions retain reference arithmetic. CPU executables
link only the optimized CPU library. The separate MPFR resident executable
instantiates the production **checked512-bit** adapters and links the CPU
library plus libmpc_cuda; it never calls a public accelerator entry. Runtime200
is a correctness target, not a performance series. Both build systems build it.

Device operands, scalar constants, streams, events, tree/scratch, pivots, MPFR
stack and bounded arenas are allocated/packed before timing. Operands are
restored before each run. CUDA events surround only the shared device-pointer
schedules (`Rblas1_jobs_common.h`, `Rtrisolve_jobs.h`, `Rinverse_schedule.h`),
on their actual launch stream. BF AXPY/COPY use the existing lane engine;
MPFR uses the inverse-local checked element job. The stop event precedes final
status/download validation. No CUDA allocation, pack/unpack or operand transfer is
inside the interval. Inverse singular-status feedback, launch gaps and scheduler
latency between kernels **are** included: this is elapsed resident computation,
not the sum of isolated kernel durations. Arenas/stack use the production
scoped policy; memory checks cover operands, tree, workspace and arena. There
is no arithmetic optimization or gate/default change.

Data use the existing deterministic full-precision generator, independently
seeded by n. Vectors use unit strides; AXPY uses alpha=1. Solves use
Upper/NoTrans/Nonunit/inc1; TBSV has k=64, including
n<64. The shared CSV solve form `UNN_k64` records this band parameter;
it is ignored for full/packed TRSV/TPSV. Matrix storage is native
full/band/packed. Inverses use Upper/Nonunit,
well-conditioned factors, reference iMlaenv blocking; GETRI has identity pivots
and optimal workspace. Lower/T/C/unit/stride/pivoted correctness is checked by
the focused production tests, not additional performance curves.

Default repetitions: warmup1 + best3 at matrix n<=2000 / vector n<=2^21;
larger sizes one measurement (CPU without warmup, GPU with one untimed
warmup). GPU modules load eagerly before event intervals: no first-launch
module-loading latency masquerades as kernel time. A warmup exceeding10s
reduces to one measured call; a measured call exceeding10s ends further
repetitions. Options `--reps R --warmup W` override defaults for pilots/checks,
not this10s guard;
`--dump FILE` writes actual output component representations outside timing.
Failures have nonzero exits, never a fabricated zero time.

FLOPs (conventional high-precision operations, not hardware-double operations):
AXPY/DOT2n, complex DOT8n, ASUM n (complex2n), NRM2 nominal2n (complex4n,
scaled-tree implementation does more work), TRSV/TPSV n² (complex4n²),
TBSV n(2k'+1)-k'(k'+1), k'=min(64,n-1), complex×4. TRTRI n³/3+2n/3,
LAUUM n³/3+n²/2+n/6, POTRI leading2n³/3, GETRI leading4n³/3; complex×4.
COPY and IAMAX do not have meaningful FLOP counts: their CSV GFLOPS is0 and
plots use **G elements/s**, not fictitious GFLOPS. Compare times for speedup.

```
# Configure both CUDA libraries + benchmarks first; CPU32, GPU1 host thread.
make -C benchmark Rkernel_bench.{dd,td,qd,mpfr}_opt \
  Rkernel_bench.{dd,td,qd}_cuda_kernel Rkernel_bench.mpfr_cuda_resident
# CMake targets: bench_Rkernel_bench_<b>_opt / <b>_cuda_kernel,
#               bench_Rkernel_bench_mpfr_cuda_resident
python3 benchmark/test_kernel_new.py
python3 benchmark/test_kernel_new.py --bindir benchmark --type mpfr512
D=benchmark/results/2026/x86_64-pc-linux-gnu-Ryzen-Threadripper-3970X-TeslaA100/kernel-bench-2
python3 benchmark/kernel_bench.py new --out "$D" --phase pilot --methods cpu
# After the maintainer's GPU gate and focused tests/sanitizer:
python3 benchmark/test_kernel_new.py --gpu-check --bindir benchmark \
  --mpfr-bindir benchmark --out "$D/checks"
python3 benchmark/kernel_bench.py new --out "$D" --hours 5.8
# Optional: use remaining budget to replace conservative MPFR caps by actual10s stops.
python3 benchmark/kernel_bench.py new --out "$D" --types mpfr512 --methods gpu,cpu --phase extend-mpfr
PYTHONPATH=/home/docker/plotlib python3 benchmark/kernel_bench.py new --out "$D" --phase plot
python3 benchmark/kernel_bench.py new --out "$D" --phase bundle  # after all runners stop
```

`kernel_new.py` runs sequentially, never pins cores, checks the external
`/tmp/gemm-sweep-forms.done` marker (never creates it), and logs uptime before,
after and every30s, plus GPU compute applications for GPU runs. Every point's
CSV records maximum sampled load, load>=16, other GPU application, exit status,
phase, and point wall time (subprocess plus condition probes). Load includes the benchmark's own threads;
it is not an estimate of external CPU utilization. Contaminated points are ×
on graphs; no unmeasured slowdown correction is applied. CPU requests32 OpenMP
threads but retains library-internal serial/SIMD paths and cutoffs. Files in `raw/` retain
commands/output/exit receipts. `bundle` consolidates owned point logs by series,
retaining every command/output/exit, without changing CSV data. Pilot points are
not plotted as sweep points.

Two pilot sizes per series drive per-call10s MPFR/30s BF forecasts and a global
wall-budget forecast. All point wall time (including setup/pilots/probes) is
charged against a cumulative <=6h budget, retained across restarts. An optional
`budget-overhead.json` receipt can conservatively charge interruptions or older
records that did not include probe overhead; this never changes kernel times. Default
matrix grid is reduced to powers64..4096 and nearby multiples7/31 plus7/31;
`--dense` requests every multiple7/31. Vectors request all57 powers/quarter
points1024..2^24. Forecast caps, global thinning and actual MPFR>10s stops are
explicit in plan/forecast JSON and reductions.log; see the result README for
what actually ran. Restart resumes completed points, including failed receipts.
`extend-mpfr` advances existing MPFR curves toward their requested maxima,
stopping at the first reported>10s point/error or insufficient remaining budget
(twice the next point's size-scaled cost forecast, minimum60s reserve).
`extend-inverses` limits that operation to inverses. Extensions are logged;
this avoids treating a pessimistic small-size pilot as proof that1024/2048
cannot be measured.

`test_kernel_new.py --gpu-check` compares the actual benchmark operands at
n64/200 to CPU32: exact representations (including zeros) except reductions,
which additionally allow the production test envelopes (dot/asum32n²u_eff,
norm(16+8ceil(log2(2n)))u_eff|reference|, u_eff=4u). It does not replace the
production reference/oracle/error-path tests. No GPU checks run during sweeps.


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
