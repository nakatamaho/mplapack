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
regenerated identically in every binary. Copies and initialization are outside
the timer. No warmup; best of three at n <= 2000, one run above. CSV columns:
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
