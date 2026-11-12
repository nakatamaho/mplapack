# Kernel benchmark: goal 13 continued (flag wavefront solves, one-kernel reductions, block leaves / tiled TRMM)

After-measurements of the goal 13 continuation (items A-C of
`bench_local/pi_instructions/goal13_continue.md`) on the A100 80GB PCIe,
CUDA events around the device-pointer schedules, `OMP_NUM_THREADS=1`.
"Before" for the vectors/solves = the goal 11 sweep
(`../kernel-bench-2/kernels.csv`); for the MPFR inverses = this session's
start (the goal 12 mpc_cuda v0.0.2-dev state, re-measured).  Raw rows:
`itemd_after.csv`, merged with kernel-bench-2 into
`kernels-before-after.csv`; the graphs (plot_kernel_bench.py, extended
with the goal 13 kernel names) are in `graphs/`; the roofline
(roofline.py, peaks = dd/td/qd from `../roofline/peaks.csv` + the MPFR
level-2 peaks of `../roofline/mpfr/peaks-v0.0.2-dev.csv`, level-1 variant
in `roofline-level1/`) in `roofline/`.

**CPU caveat (unchanged):** the host load average stayed at ~100/64 all
session (maintainer instruction of 11/10 18:55); CPU timings remain
suspended and every CPU target is marked unverifiable.

## GPU before/after (kernel only)

| point | before (goal 11 / goal-13 part 1) | after (this session) | speedup |
|---|---|---|---|
| dd iRamax inc1 2097152 | 0.882 ms | 0.046 ms | 19.13x |
| dd iCamax inc1 2097152 | 0.825 ms | 0.065 ms | 12.79x |
| dd iRamax inc1 16777216 | 6.316 ms | 0.204 ms | 30.99x |
| dd Rdot inc1 2097152 | 0.791 ms | 0.090 ms | 8.77x |
| dd Rdot inc1 16777216 | 5.637 ms | 0.434 ms | 12.98x |
| dd Cdotc inc1 2097152 | 0.767 ms | 0.143 ms | 5.35x |
| dd Rasum inc1 2097152 | 0.794 ms | 0.073 ms | 10.92x |
| dd Rnrm2 inc1 2097152 | 0.882 ms | 0.298 ms | 2.96x |
| dd Rtrsv UNN_k64 1024 | 6.898 ms | 1.835 ms | 3.76x |
| dd Rtrsv UNN_k64 2048 | 14.993 ms | 4.008 ms | 3.74x |
| dd Rtrsv UNN_k64 4096 | 30.230 ms | 8.273 ms | 3.65x |
| dd Rtbsv UNN_k64 2048 | 13.667 ms | 3.615 ms | 3.78x |
| dd Rtpsv UNN_k64 2048 | 14.743 ms | 4.090 ms | 3.60x |
| dd Ctrsv UNN_k64 2048 | 18.779 ms | 7.481 ms | 2.51x |
| td iRamax inc1 2097152 | 1.300 ms | 0.077 ms | 16.93x |
| td Rdot inc1 2097152 | 1.191 ms | 0.125 ms | 9.53x |
| td Rtrsv UNN_k64 2048 | 18.049 ms | 6.918 ms | 2.61x |
| qd iRamax inc1 2097152 | 1.754 ms | 0.094 ms | 18.62x |
| qd Rdot inc1 2097152 | 1.634 ms | 0.213 ms | 7.67x |
| qd Rtrsv UNN_k64 2048 | 22.634 ms | 11.111 ms | 2.04x |
| mpfr512 iRamax inc1 2097152 | 6.884 ms | 1.141 ms | 6.04x |
| mpfr512 Rdot inc1 2097152 | 10.407 ms | 4.294 ms | 2.42x |
| mpfr512 Rtrsv UNN_k64 2048 | 173.895 ms | 86.714 ms | 2.01x |
| point | before (goal-13 part 1 / goal 12 state) | after (this session) | speedup |
|---|---|---|---|
| mpfr512 Rtrtri UN 512 | 4.320 s | 2.368 s | 1.82x |
| mpfr512 Rtrtri UN 1024 | 24.332 s | 20.502 s | 1.19x |
| mpfr512 Rlauum UN 512 | 3.019 s | 2.999 s | 1.01x |
| mpfr512 Rlauum UN 1024 | 8.156 s | 8.115 s | 1.01x |
| mpfr512 Rpotri UN 512 | 7.341 s | 5.365 s | 1.37x |
| mpfr512 Rpotri UN 1024 | 32.483 s | 28.575 s | 1.14x |
| mpfr512 Rgetri UN 512 | 4.533 s | 2.579 s | 1.76x |
| mpfr512 Rgetri UN 1024 | 25.363 s | 21.484 s | 1.18x |

## Roofline (peaks-merged.csv: dd/td/qd level-2 peaks + MPFR v0.0.2-dev level-2 peaks)

| point | % of roofline | bound |
|---|---:|---|
| dd iRamax 2^21 | 52.0 | memory |
| dd iRamax 2^24 | 94.1 | memory |
| dd Rdot 2^21 | 53.2 | memory |
| dd Rdot 2^24 | 88.3 | memory |
| dd Rasum 2^21 | 33.0 | memory |
| dd Rtrsv 2048 | 0.6 | memory |
| td/qd Rdot 2^21 | see roofline.csv | memory |
| mpfr512 iRamax 2^21 | 10.5 | memory |
| mpfr512 Rdot 2^21 | 5.6 | memory |
| mpfr512 Rtrsv 2048 | 0.18 | arith |
| mpfr512 Rtrtri 512 / 1024 | 0.07 / 0.06 | arith |
| mpfr512 Rlauum 512 / 1024 | 0.05 / 0.16 | arith |
| mpfr512 Rgetri 512 / 1024 | 0.25 / 0.24 | arith |

Notes: the MPFR percentages use the v0.0.2-dev level-2 peaks (faster
arithmetic than the goal 11 peaks, so equal times give lower percentages
than the goal 11 README numbers).  The MPFR inverse percentages are far
from the arithmetic bound: the leaves are now block-parallel but the mm
column lanes still run one thread per column with a serial k sweep, the
TRSM is lane-per-column, and the complex lauu2 keeps the serial routine —
the documented remaining work of item 4.

## Targets (continued session)

| item | target | status |
|---|---|---|
| A flag wavefront | dd Rtrsv 2048 >= 10x vs goal 11 (<= 1.5 ms) | MISSED: 4.01 ms (3.7x vs goal 11). The reference chain costs one dd division (~0.53 us) plus one cross-SM flag round trip (~0.7-0.9 us) per column: floor ~1.3-1.5 us/col. Measured breakdown in the commit body. |
| A flag wavefront | faster than before (8.18 ms) | MET: 4.01 ms (2.0x); td 6.98 ms, qd 11.1 ms; MPFR512 87.4 ms (2.1x) |
| B one-kernel reductions | dd Rdot 2^21 >= 50% bandwidth | MET: 55% (91 us); 2^24: 92%. td 60%; qd 47% (just under); Rnrm2/RCnrm2 ~17% (serial scaling chain, order-fixed) |
| B one-kernel reductions | CPU >= 8x default | SUSPENDED (host load) |
| C inverses | GPU faster than CPU at n >= 1024 | UNVERIFIABLE (CPU timings suspended); GPU MPFR512 trtri/getri 1.8x faster than before at 512 |
| D benchmark/graphs/roofline | kernel-bench-3 products | DONE: itemd_after.csv, kernels.csv, kernels-before-after.csv, graphs/ (81 plots), roofline/ + roofline-level1/ (peaks-merged CSVs, roofline.csv, summary, PNGs) |
