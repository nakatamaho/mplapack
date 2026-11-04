# Rfactor baseline benchmark

`Rfactor.<backend>`, `Rfactor.<backend>_opt`, and
`Rfactor.<backend>_cuda_total` time the linked library without dlopen or changing
legacy Rgetrf/Rpotrf output formats. Both build systems build these alongside
the existing benchmarks. CUDA factor timings use the public dispatcher:
eligible calls are resident; rejected entries use optimized CPU factors and
CPU BLAS. CUDA binaries return 77 without a device, wait for CUDA completion
inside the timer, and fail on a sticky CUDA runtime error. A CUDA-library
timing by itself is not evidence that the entry computed on the GPU.
Discard timings from failing runs.

Examples:

    OMP_NUM_THREADS=32 ./Rfactor.dd_opt --n 512 --n 1024 --n 2048 --n 4096
    OMP_NUM_THREADS=32 ./Rfactor.mpfr_cuda_total --potrf --n 1024
    MPLAPACK_LAPACK_FAST=1 ./Rfactor.qd_opt --complex --max 600

`--potrf` selects Cholesky (lower triangle); the default is LU. `--complex`
selects Cgetrf/Cpotrf. MPFR precision is fixed at 512 bits before allocating
operands. The deterministic, Hermitian strictly diagonally dominant input is
regenerated identically in every binary. Copies and matrix initialization are
outside the timer. No warmup; best of three at n <= 2000, one run above. CSV columns:
`routine,n,repeats,seconds,gflops`; leading LAPACK operation counts are 2*n^3/3
for real LU and n^3/3 for real Cholesky, times four for complex.

Without `--n`, sizes are the union of powers of two and every multiple of 13
up to `--max` (default 2048). Repeated `--n` selects only those sizes.
Environment controls (including the fast-mode and GPU BLAS algorithm/threshold
switches) are not overwritten; record them and thread count with measurements.
Reference libraries do not implement MPLAPACK_LAPACK_FAST.

For interleaved CPU-library / default CUDA / forced CUDA medians:

```sh
python3 benchmark/run_factor_dispatch_bench.py --bin-dir build/benchmark \
  --sizes 512 1024 2048 4096 --rounds 3 --threads 32 --output factors.jsonl
# Add --routines gesv posv --rhs 1 0 for driver gates.
```

The runner clears MPLAPACK environment overrides, alternates order, and
retains the benchmark's repeat policy (use `--repeat 1` to include cold-call
initialization in every sample). In particular a first timed
`cudaDeviceSynchronize` can initialize the runtime even after a CPU-only
factor; do not mistake that fixed benchmark cost for GPU BLAS activity.
The MPFR benchmark uses 512 bits. Use a fresh output file after changing a
threshold/build; completed rows are skipped on restart.

## Cold/warm calls and inverses (batch 09/10)

`--samples --repeat 3` emits **every** call as
`routine,n,nrhs,call,temperature,seconds`, with first-call `cold` and subsequent
`warm` labels. Unlike the historical best-of-N output, the cold sample includes
both the timed numerical call/synchronization and the initial CUDA runtime/device
probe; warm samples include neither that probe nor an artificial warmup. Input
copies remain outside all timers. Use one size per fresh process for cold studies
(as the runner does); with multiple `--n`, only the first size pays CUDA setup.
The labels mean first/later call for that size, not a flushed CPU/cache state.

`--trtri`, `--lauum`, `--potri`, and `--getri` time inverse/product routines on
synthetic nonsingular triangular/Cholesky/LU factors; **no preceding factorization
is timed**. GETRI uses optimal workspace, small unit-L subdiagonals (rather
than an ill-conditioned dense random unit L), and deterministic nonidentity
pivots.
TRTRI uses nonunit diagonals. `--upper` selects U instead of L for all triangular
routines (including Cholesky and POSV). The inverse leading operation counts are
n^3/3 (TRTRI/LAUUM), 2*n^3/3 (POTRI), and 4*n^3/3 (GETRI), times four for complex.
`--entry` requires acceptance of every inverse/solver call. MPFR inverse GPU
entries do not exist. Its default GETRI can still use legacy per-call GPU BLAS
in the reference outer algorithm; do not label that path a resident GPU inverse
or assume every MPFR inverse comparison is CPU-only.

```sh
python3 benchmark/run_factor_dispatch_bench.py --bin-dir build/benchmark \
  --backends dd td qd --sizes 256 512 1024 2048 4096 \
  --routines trtri lauum potri getri --samples --repeat 3 \
  --rounds 3 --threads 32 --output inverse-calls.jsonl
# Repeat triangular routines with --upper and a separate output file.
# MPFR inverse comparison: --backends mpfr --modes cpu default.
python3 benchmark/test_factor_dispatch_bench.py
```

In sample mode the runner retains all calls, reports separate median cold and
warm columns (median warm calls within a process, then median across rounds),
and requires one successful forced-entry trace **per call**, not just one per
process. Old best-of-N records and per-call records have distinct resume keys.
No correctness jobs or builds may overlap threshold measurements.

The [cut-short CPU-library gate retune](results/2026/x86_64-pc-linux-gnu-Ryzen-Threadripper-3970X-TeslaA100/gate-retune/README.md)
archives completed JSONL, reproducible cold/warm CSV, current/previous gates,
load exclusions and unmeasured inverse cases. It does not claim post-retune
performance or elimination of residual CPU-fallback overhead.
`test_factor_gate_smoke.py --bin-dir build/benchmark --mpfr-bin-dir build-mpfr/benchmark`
checks default/override/disabled dispatch traces, not performance or numerical
identity; use the dedicated fallback oracle for the latter.
