# New BLAS1, triangular solves and MPFR inverses — resident GPU kernels

**Complete: 7,941 accepted sweep points, 180 curves, 76 graphs; no failed points.**

[Representative timings and all graphs](ratios.md) · [Full ratios](ratios.csv) ·
[Coverage](coverage.csv) · [Explicit omissions](grid-omissions.md) · [Summary](summary.json)

AMD Ryzen Threadripper3970X (32 cores/64 threads), NVIDIA A10080GB PCIe
(sm80), CUDA13.1, GCC15.2, Ubuntu26.04. November4–5,2026 (JST).
See `environment.log` for actual machine/toolchain/configure output.
Production prerequisites: `c1e4c603b`, `6caad21c7`, `fc85906fc`.
This benchmark changes no library arithmetic, public API, default or gate.

## Scope and timing

76 backend/routine pairs,180 curves:11 BLAS1 and6 triangular solves for
dd/td/qd/MPFR512, plus8 real/complex MPFR inverses. No BF inverse
rebenchmark, runtime-precision performance sweep, scalar rotg/rotm/rotmg,
or stage3 dispatcher/threshold integration.

- GPU: operands, constants, scratch, streams, events, precision/arena/stack
  state are prepared before timing. CUDA events bracket the shared production
  device schedules on their actual stream. No operand packing, upload,
  restoration, download or CUDA allocation in the interval. Inverse scheduling's
  required small scalar control feedback remains included; this is the
  resident algorithm's event span, not a sum of isolated instruction times.
  It is **not** public-entry/end-to-end latency. MPFR512 uses a separate
  checked resident-event adapter, not the old CUPTI/public-entry path.
- CPU: optimized CPU library, `OMP_NUM_THREADS=32`, `OMP_DYNAMIC=FALSE`,
  no affinity. This requests32 threads; library-internal serial/SIMD paths
  and cutoffs are retained, not every routine actually uses32 threads.
  Normal library-internal allocation and OpenMP scheduling remain in the
  CPU call time; operand restoration is outside it.
  Only the seven sum reductions explicitly enable
  `MPLAPACK_MPBLAS_FAST=1`; their CPU default/reference-order curve is also
  measured. Inverse leaves never inherit FAST; LAPACK FAST is off.
- Small sizes (matrices<=2000, vectors<=2^21): one warmup, best of3.
  Larger sizes: one measured call. GPU always has one untimed warmup and
  eager module loading, including large sizes; large CPU calls have no
  warmup. A warmup>10s reduces toone measured call; a measured call>10s
  ends further repetitions. MPFR series stop beyond a reported>10s point;
  initial conservative caps and later MPFR extensions are described below.
- Fully reproducible input seed20261104+n, full backend precision, nonzero
  low limbs. Vectors have unit strides; AXPY uses alpha=1. Solves use Upper/NoTranspose/Nonunit,
  full/native-packed/native-band storage. TBSV's k=64, lda=65 (including
  small n); the CSV's common `UNN_k64` label ignores k for TRSV/TPSV.
  Matrix diagonals are real positive and dominant. Inverses use Upper;
  GETRI receives generated LU factors and identity pivots, not an unfactored
  matrix. Untouched storage is included in operand comparisons.

FLOPs are conventional algorithm counts, not hardware instruction counts:
AXPY/DOT2n, complex DOT8n, real/complex ASUM n/2n, NRM2 2n/4n;
TRSV/TPSV n² (complex4n²); TBSV n(2h+1)-h(h+1), h=min(64,n-1)
(complex4×); TRTRI n³/3+2n/3, LAUUM n³/3+n²/2+n/6;
POTRI leading2n³/3, GETRI leading4n³/3 (complex4×).
COPY and IAMAX do not have meaningful FLOP counts: their graphs deliberately
show G elements/s rather than inventing GFLOPs; raw FLOP-count fields are0.
Ratios always use actual CPU/GPU seconds, independent of the display unit.

## Grids, budget and data provenance

Vectors request all57 powers/quarter-interval points1024..2^24. The dense
matrix request is every multiple7/31 and power2 through4096. To bound the
many expensive process/setup and matrix calls, the baseline uses23 points:
7,31, powers64..4096 and the nearest lower multiples7/31 at those powers.
Thus small powers1..32 and interior multiples are intentionally omitted.
Initial per-series pilot forecasts cap MPFR calls at10s (BF30s). The optional
post-baseline MPFR phase advances these curves using actual observations
until a reported>10s point, the requested maximum, or insufficient budget
reserve (twice the next point's size-scaled cost forecast, minimum60s). No extrapolated point is plotted.

`kernels.csv` is authoritative. Only successful `phase=sweep` rows appear in
plots/tables. Current `pilot` rows calibrate the planner; they are not plotted.
`plan-*.json`, `forecast-*.json` and `reductions.log` document initial caps,
possible thinning, actual cutoff and extensions. [Coverage](coverage.csv)
reports actual points, not merely planned points. All BF curves and all
MPFR solves complete the reduced grid through4096. All vector curves reach
2^24 except CPU-fast MPFR RCnrm2: its n=14680064 call took10.9495s, so the
final16777216 point was omitted. All remaining inverse omissions follow
reported>10s points; [the exact lists](grid-omissions.md) include every omitted
baseline size. The dense matrix grid would have712 points rather than23.
No global forecast thinning was needed.

Budget is cumulative across restarts, default5.8h (maximum6h), including
pilots, setup, discarded attempts and condition probes—not just kernel
seconds. Early rows recorded subprocess wall only; `budget-overhead.json`
conservatively charges the entire first-pilot→final-benchmark-validation
window:22:02:39 November4 to03:06:41 November5, plus2s timestamp allowance.
This is18,244s (**5.068h**, below5.8h and6h), including build/check/wait pauses,
all excluded attempts and post-sweep comparisons. Recorded point wall time
is14,071.291s; the conservative additional charge is4,172.709s. The last
performance point finished02:51:11. Plot elapsed annotations sum accepted
rows for all curves in that graph, excluding pilots/discards/shared overhead.

Two benchmark-only review corrections are retained transparently:

1.16 preliminary inverse CPU pilots inadvertently allowed FAST in leaf DOT
   calls. Corrected **before any inverse sweep**, retained as
   `pilot-invalid-fast`, and repeated under the default arithmetic policy.
2.65 unwarmed large GPU points included a material first-launch artifact
   between CUDA events. Retained as `sweep-invalid-cold-launch`, excluded,
   and repeated with eager module loading and one untimed GPU warmup.
   The already-warmed152 GPU pilots are retained as `pilot-pre-eager` and
   repeated to calibrate final startup/probe overhead. Already-warmed small
   sweep points remain valid. Example: DD AXPY n2621440 changed from the
   excluded1.23728ms to the corrected0.092352ms; no kernel math changed.

Only phase/file labels were changed; raw numbers/output were preserved.
`metadata-corrections.log` lists every relabelled record. Two controlled
interruptions produced no accepted point; partial logs are in `validation/`.
All their wall cost is still charged. Bundled `raw/` logs retain each point's
command, output and exit receipt; `checks*/runs.log` retains comparison logs.

## Acceptance and shared-host conditions

The external GPU marker was first observed present22:49:41JST November4;
it was never created or modified here. No GPU computation before it.
CPU timing was stopped before the FIRST GPU work: sequential item1–3
focused tests, then MPFR512 inverse memcheck. No GPU test overlaps a
performance sweep.

Completed acceptance (receipts in `validation/` and `checks*/`):

- Both autotools and CMake benchmark builds, sm80+sm90 cubins verified.
  An autotools device-link architecture omission was corrected; the final
  link uses both architectures. CUDA dependency files also track shared
  production headers. Final metadata-review rebuilds and checks passed.
- Focused GPU tests12/12 PASS: MPFR inverse512/200, BF BLAS1, MPFR
  BLAS1 512/200, BF triangular solves, MPFR solves512/200.
- MPFR512 inverse memcheck:132 direct calls,790 checks,66264 component
  comparisons,0 failures; `ERROR SUMMARY: 0 errors`.
- Three benchmark operand-check campaigns: **912/912 component-identical**,
  each76 pairs × n=64/200 ×2 builds. The initial304 were single-call checks;
  primed and final reviewed campaigns each requested warmup1/reps3 (the
  >10s guard applies), testing restoration. Maximum component difference0.
- Host BF9/9 and MPFR6/6 focused tests; CPU CLI smoke in all4 backends
  in both builds, including inverse FAST/default equivalence;16 reused
  CPU/CUDA archives pass unique-symbol checks. No production library source changed.
- Six policy tests pass. Historical MPFR512 CSV redraws with old/new
  `plot_sweep.py` defaults are byte-identical in light and dark styles.
- Bundling preserved the CSV SHA256:8,534 point command/output/exit receipts
  and1,824 comparison-call receipts are retained.

Before/after each point and every30s during longer points, `conditions.log`
records uptime and (for GPU points) compute applications. Other processes
are identified by PID, not executable name; own just-reaped PIDs are excluded.
CSV records before/after/max sampled load and load>=16/other-app flags.
CPU-only points do not query GPU applications (`gpu_busy=0` is then not
applicable, not a claim of GPU idleness); the initial CPU sweep overlapped
another owner's GPU GEMM sweep while respecting the GPU gate.
These flags are × on plots. Load includes our own threads; it is not an
estimate of external utilization. No slowdown correction is applied.

Accepted-point sampled one-minute load:0.979–14.908; no sample>=16.
All3,163 GPU sweep points have no observed other-GPU-application flag.
Thus no accepted plot point has a condition ×. Sampling is not proof of
continuous CPU idleness or absence of short-lived interference.

## Reproduce

From the repository root, with its QD/MPFR/GMP/MPC/CUDA dependencies built:

```sh
./gen_configure.sh && ./config.status
make -C benchmark -j8 Rkernel_bench.dd_opt Rkernel_bench.td_opt \
  Rkernel_bench.qd_opt Rkernel_bench.mpfr_opt \
  Rkernel_bench.dd_cuda_kernel Rkernel_bench.td_cuda_kernel \
  Rkernel_bench.qd_cuda_kernel Rkernel_bench.mpfr_cuda_resident
# Existing CUDA-enabled CMake configurations used here:
cmake --build /tmp/factor04-cmake -j8 --target \
  bench_Rkernel_bench_dd_opt bench_Rkernel_bench_td_opt bench_Rkernel_bench_qd_opt \
  bench_Rkernel_bench_dd_cuda_kernel bench_Rkernel_bench_td_cuda_kernel \
  bench_Rkernel_bench_qd_cuda_kernel
cmake --build /tmp/factor01-cmake -j8 --target \
  bench_Rkernel_bench_mpfr_opt bench_Rkernel_bench_mpfr_cuda_resident
python3 benchmark/test_kernel_new.py
# Only after the external marker; strictly sequential, before timing:
set -e
test -f /tmp/gemm-sweep-forms.done || exit 1
export PATH=/usr/local/cuda/bin:$PATH
export LD_LIBRARY_PATH=$PWD/external/i/QD/lib:$PWD/external/i/MPFR/lib:$PWD/external/i/GMP/lib:$PWD/external/i/MPC/lib:/usr/local/cuda/lib64
export OMP_NUM_THREADS=1
ctest --test-dir /tmp/factor01-cmake -j1 -R '^mpfr_cuda_inverse_device_(512|200)$' --output-on-failure
ctest --test-dir /tmp/factor04-cmake -j1 -R '^(dd|td|qd)_cuda_Rblas1_device$' --output-on-failure
ctest --test-dir /tmp/factor01-cmake -j1 -R '^mpfr_cuda_Rblas1_device_(512|200)$' --output-on-failure
ctest --test-dir /tmp/factor04-cmake -j1 -R '^(dd|td|qd)_cuda_Rtrisolve_device$' --output-on-failure
ctest --test-dir /tmp/factor01-cmake -j1 -R '^mpfr_cuda_Rtrisolve_device_(512|200)$' --output-on-failure
MPFRXX_DEFAULT_PRECISION_BITS=512 compute-sanitizer --tool memcheck --error-exitcode 9 /tmp/factor01-cmake/mpfr_cuda_inverse_device_test --small
python3 benchmark/test_kernel_new.py --gpu-check \
  --bindir /tmp/factor04-cmake/benchmark --mpfr-bindir /tmp/factor01-cmake/benchmark \
  --out /tmp/new-kernel-checks
D=benchmark/results/2026/x86_64-pc-linux-gnu-Ryzen-Threadripper-3970X-TeslaA100/kernel-bench-2
python3 benchmark/kernel_bench.py new --out "$D" --methods gpu,cpu --hours 5.8 \
  --bindir /tmp/factor04-cmake/benchmark --mpfr-bindir /tmp/factor01-cmake/benchmark
python3 benchmark/kernel_bench.py new --out "$D" --methods gpu,cpu --types mpfr512 \
  --phase extend-mpfr --bindir /tmp/factor04-cmake/benchmark \
  --mpfr-bindir /tmp/factor01-cmake/benchmark
PYTHONPATH=/home/docker/plotlib python3 benchmark/kernel_bench.py new --out "$D" --phase plot
python3 benchmark/kernel_bench.py new --out "$D" --phase bundle
python3 "$D/validation/derive_tables.py"
PYTHONPATH=/home/docker/plotlib python3 "$D/validation/legacy_plot_check.py"
```

Use a **new output directory** for independent measurements; an existing
CSV resumes/skips completed and failed points and retains the cumulative
budget. No numerical tests may run concurrently with these timing commands.
For autotools operand checks, use `--bindir benchmark --mpfr-bindir benchmark`
and the external-library paths in `benchmark/README.kernel_bench.md`.

## Results

The [76-row representative table](ratios.md) links every graph and gives
actual CPU/GPU milliseconds and their ratio. It prefers n=2048, then1024;
when those were not measured, it uses the largest common power of two.
For reductions the main CPU comparison is FAST; [ratios.csv](ratios.csv)
also includes CPU-default comparisons and vector n=2^21 results.
A ratio>1 means a faster GPU; no time is inferred for omitted sizes.

At vector n=2^21, CPU/GPU ratios across the11 routines range:

| Type | GPU speedup versus selected CPU curve |
|---|---:|
| DD | 9.84–107× |
| TD | 11.1–102× |
| QD | 4.91–115× |
| MPFR512 | 3.03–73.9× |

These are resident-only speedups, not transfer-inclusive application gains.
CPU FAST is not uniformly faster than CPU default; both curves are retained.

The slow end is equally important. At n=2048, DD Rtbsv's ratio is0.0300×
(about33× slower GPU); TD/QD Rtbsv ratios are0.0792×/0.1057×.
All eight MPFR inverse representative ratios are below1:

| Routine | n | CPU/GPU |
|---|---:|---:|
| Rtrtri | 512 | 0.0304× |
| Rlauum | 512 | 0.0535× |
| Rpotri | 256 | 0.0290× |
| Rgetri | 512 | 0.0840× |
| Ctrtri | 256 | 0.0220× |
| Clauum | 256 | 0.0262× |
| Cpotri | 128 | 0.0187× |
| Cgetri | 256 | 0.0447× |

Their actual cutoff points—not pessimistic pilot extrapolations—limit the
upper sweep. For example, GPU Rtrtri at1022 took33.485s, Ctrtri at496
26.655s, and Cpotri at n=248 took15.375s. CPU inverse calls also crossed10s.
These observations do **not** change any dispatcher default or threshold.

## UNVERIFIED / deliberately not done

- Full repository/LAPACK QA; exhaustive benchmark flags, strides and inputs;
  component comparison at every large performance point (benchmark checks
  explicitly use n=64/200).
- sm90 hardware execution (only sm80 A100 was available; both cubins built).
- MPFR performance at precisions other than512; correctness was checked at
  fixed512 and runtime200, not every fixed precision or FORCE_RUNTIME512.
- Dense omitted matrix points and upper sizes omitted by the10s policy.
- Large-size memcheck; the required inverse memcheck used `--small`.
- End-to-end transfer-inclusive speedups, cold public-entry latency, or
  stage3 gate/default/LAPACK-path integration. No dispatcher retuning.
