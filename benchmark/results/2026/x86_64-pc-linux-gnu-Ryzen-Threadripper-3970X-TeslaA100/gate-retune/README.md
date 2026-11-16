# CUDA LAPACK gate retune — measurements cut short

Finalized from existing measurements only, 2026-11-04. No new performance sweep.
Numerical source HEAD: `93d1a70fb50d36227a08c19eeb894caa52badf14` (small-n LU
team-cap revert), plus the Rfactor sample/inverse extensions committed with
this report. Subsequent `c70d4812e` / `ff91620e6` only filed benchmark results.
All timings here are **before** the gate changes, not measured new-default times.

## Method and retained data

Threadripper 3970X, A100 80GB PCIe; 32 OpenMP threads, OMP_DYNAMIC=FALSE;
MPFR 512 bits. Hardware/compiler details, session metadata and observations
are in `raw/`. CPU means the optimized **CPU library**, not CUDA with a switch
disabled. Default means the CUDA-library public dispatcher. Forced uses
LAPACK_CUDA_MIN=0 and, for drivers/inverses, the bool entry; each call must
emit a successful resident-entry trace. Deterministic inputs, same binaries
per backend/mode; commands, executable hashes, timestamps, stdout/stderr and
individual calls are preserved in the JSONL. Hashes identify executables,
not all shared dependencies. The original builds were checked in both CMake
and autotools before measurement (see batch-09-10 receipts).

Each process executes two calls on reset input. Three interleaved/reversed
CPU/default/forced rounds; inverses at 2048/4096 reduced to two rounds, retaining
already completed extras. `summary.csv`: median first-call time (`*_cold`),
median later-call time (`*_warm`), min/max of per-process warm medians and entry
trace counts. With two calls there is one warm sample per process. CUDA cold
includes the initial runtime/device probe; copies/initialization are outside
the timer. GFLOPS and cold time are not used to choose gates.

The host was shared, not uniformly idle. In particular td small-n CPU LU
(256/512) is grossly inflated by host load; these points are retained but
**excluded from crossover decisions** (see optimized/common/README.md).
Lower td Cholesky is also much slower than Upper in that measurement interval;
we do not use the apparent Lower 256 win to lower its gate. Decisions require
both triangles, stable warm ranges and approximately >=1.3x improvement.
No additional timing was taken to clean up noisy points.

## Decisions (minimum n, opt-in = LONG_MAX)

Only existing gate constants changed: no arithmetic, kernels, scheduling or
new triangle/type-specific dispatch logic. Overrides retain their meanings.

| Backend | Gate | Old | New | Evidence / conservative boundary |
|---|---|---:|---:|---|
| qd | R/Cgetrf | 1024 | 512 | at 512 CPU/GPU 1.48x / 1.60x; 256 only 1.16x / 1.17x |
| td | Rpotrf | 2048 | 512 | Upper 512/1024 1.74x/2.41x; Upper 256 only 1.05x; Lower wins too |
| qd | Rpotrf | 2048 | 512 | at 512 L/U 1.62x/2.06x; Lower 256 only parity |
| dd | Cpotrf (shared U/L) | 4096 | opt-in | Upper 4096 GPU/CPU 2.25x slower |
| td | Cpotrf (shared U/L) | 2048 | opt-in | Upper 2048/4096 2.05x/1.40x slower |
| qd | Cpotrf (shared U/L) | 4096 | opt-in | Upper 4096 1.40x slower |
| mpfr | Rgetrf | 512 | 1024 | 512 CPU/GPU 0.302/0.303 s (parity), 1024 6.888/0.773 s |
| mpfr | Rgesv | 512 | 1024 | 512 CPU/default/GPU 0.415/1.437/1.437 s; 1024 CPU/GPU 7.066/5.338 s; 2048 28.954/20.737 s |
| mpfr | Cgesv | 1024 | opt-in | 1024 only 1.09x; 2048 CPU/GPU 134.737/161.626 s |

Unchanged: dd R/Cgetrf=2048, td R/Cgetrf=1024, dd Rpotrf=4096;
MPFR Cgetrf=1024, Rpotrf=2048, Cpotrf=4096. dd real LU 2048 is marginal
but not a measured regression; the shared real/complex constant cannot be
lowered to the clear complex 1024 win without slowing real LU. td small LU
is excluded rather than used to lower 1024. MPFR Cpotrf 4096 is retained from
the earlier CPU-library comparison; this cut-short sweep does not revalidate it.
BF fused drivers and all POSV remain opt-in. Standalone solve gates unchanged
(no new standalone-solve measurements). Drivers rejected by their gate may
still use the resident factor followed by CPU solve.

**Deliberate limitation:** BF complex Cholesky Lower wins clearly at large n,
but the existing constant is shared with Upper, which loses. Under this
constants-only task we prioritize no Upper regression: both remain opt-in,
with explicit LAPACK_CUDA_MIN available for Lower. A future separate
UPLO gate can recover those wins; this session does not add new dispatch logic.
Consequently this is not a claim that every profitable shape now dispatches
optimally, nor that all CUDA CPU-fallback overhead has been eliminated.

## CPU-fallback slowdown investigation (not fixed)

For dd Rpotrf 256/512/1024, Lower default/CPU warm ratios are about
1.52/1.40/1.25; Upper about 1.59/1.47/1.25. MPFR Rpotrf Lower 1024 is
1.493/1.203 = 1.24x; Upper 1.280/1.295 does not show the same slowdown.
These default points have **zero resident-entry traces**. The factor gate
rejects before allocation/upload; the warm sample excludes initial device
probe, so cold startup cannot explain the warm gap.

Source inspection: Rpotrf.cpp uses Rfactor_cpu_blas.h for CUDA CPU fallback;
GEMM is a private copy of the optimized CPU dispatcher, TRSM/SYRK go through
Rsolver_cpu_fallback.h (CPU *_omp then reference), and the diagonal panel is a
private CPU-bound reference recursion. Rlookahead_common.h explicitly restores
the CPU scheduling policy for these backends. Thus neither an intentional
GPU-BLAS fallback nor intentionally serial scheduling is indicated. Separate
compiled/private wrappers, compiler/code-layout differences and shared-host
load are plausible contributors, **not a proven root cause**. No profile or
new measurement was taken. Fixing this residual overhead is outside scope;
changing a gate to a marginal GPU path would mask rather than resolve it.

## Cut-short coverage

3078 completed successful process records, 348 grouped points, 347 complete
at the planned round count. Factors/drivers: 234 complete points (dd/td/qd:
60 each; MPFR: 54), real/complex, both triangles, n=256..4096 for BF and MPFR
real, through 2048 for MPFR complex, driver nrhs=1 only.
Inverses: dd 40 and td 40 complete points; qd 34 retained points, of which
Cgetri 4096 has only one round per mode. The next forced qd Cgetri 4096 process
was terminated along with the runner and measure.py; its unfinished output
was not invented as a JSONL record. Remaining qd large inverse points,
**all MPFR inverse comparisons**, Upper inverse performance and UNIT inverse
performance are unmeasured. BF inverse entries stay opt-in; MPFR inverse
entries do not exist. Completed raw JSONL was not deleted or rewritten.

## Reproduction

Recompute the table without any GPU or new measurements:

```sh
python3 benchmark/results/2026/x86_64-pc-linux-gnu-Ryzen-Threadripper-3970X-TeslaA100/gate-retune/summarize.py
```

For a future explicitly authorized measurement (not run in this session),
build the CPU and CUDA Rfactor executables from the same tree, then:

```sh
python3 benchmark/run_factor_dispatch_bench.py --bin-dir /tmp/factor04-cmake/benchmark \
  --mpfr-bin-dir /tmp/factor01-cmake/benchmark --backends dd td qd mpfr \
  --routines getrf potrf gesv posv --sizes 256 512 1024 2048 4096 \
  --samples --repeat 2 --rounds 3 --threads 32 --output NEW-lower.jsonl
# Upper: --routines potrf posv --upper, a NEW output file.
# Original MPFR complex maximum: --backends mpfr --kinds complex --sizes 256 512 1024 2048.
# Inverses: --backends dd td qd --routines getri potri trtri lauum;
# 3 rounds through 1024, 2 at 2048/4096, Lower only.
```

Do not resume into these archived files or mix old/new gates. Individual exact
original commands are in each raw record. No other builds/QA should overlap
future measurement.

## Limited verification

- Parser/resume/forced-trace tests: PASS 4/4.
- Rebuilt CMake CUDA libraries and Rfactor binaries; 36 default-boundary,
  explicit-enable and disable trace smoke checks: PASS. This is dispatch
  verification, not a new performance comparison.
- Rebuilt solver fallback oracle: dd/td/qd and MPFR 512 bits PASS 4/4;
  each 352 cases / 14,149,260 elements, zero errors or signed-zero errors.
- Autotools optimized libraries and benchmark executables rebuilt successfully;
  the same 36 dispatch checks PASS with external QD/MPFR/MPC/GMP library
  directories on LD_LIBRARY_PATH (the first attempt lacked libqd.so.2).
  Duplicate-symbol checks
  PASS on all eight rebuilt CUDA archives (four CMake, four autotools).
- Initial multi-target parallel CMake invocation raced archive linking; rerun
  with one target per invocation passed. First autotools invocation timed out;
  incremental retry passed. No product-code workaround was made.
- Full numerical QA, sanitizers, full precision sweep and post-change timing
  matrix intentionally NOT rerun; prior numerical receipts remain historical.

Quick checks (no full QA):

```sh
python3 benchmark/test_factor_dispatch_bench.py
python3 benchmark/test_factor_gate_smoke.py --bin-dir /tmp/factor04-cmake/benchmark \
  --mpfr-bin-dir /tmp/factor01-cmake/benchmark
ctest --test-dir /tmp/factor04-cmake -R '^(dd|td|qd)_cuda_solver_fallback$' --output-on-failure
ctest --test-dir /tmp/factor01-cmake -R '^mpfr_cuda_solver_fallback_512$' --output-on-failure
# Autotools uninstalled executable dependency paths:
LD_LIBRARY_PATH="$PWD/external/i/QD/lib:$PWD/external/i/MPFR/lib:$PWD/external/i/MPC/lib:$PWD/external/i/GMP/lib:/usr/local/cuda/lib64" \
  python3 benchmark/test_factor_gate_smoke.py --bin-dir "$PWD/benchmark" --mpfr-bin-dir "$PWD/benchmark"
```

Execution logs are copied into `verification/`; local build/work logs remain
under bench_local/results/finish-gates (not committed).
