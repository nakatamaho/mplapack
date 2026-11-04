# Kernel benchmark on the A100 host: CPU against GPU per precision class (November 2026)

The kernels from which LU, Cholesky, the solvers and the inverses are built
(GEMM, SYRK/HERK, TRSM, the getrf2 / potrf2 panels, laswp, trtri, lauum),
timed alone on the CPU and on the GPU for dd, td, qd and MPFR (cu_freal 256,
512, 1024 bits; cu_mpfr at the run-time precision 200 bits).  Programs,
methods and the driver: [`benchmark/README.kernel_bench.md`](../../../../README.kernel_bench.md).
This is the measurement the maintainer asked for before optimizing the
kernels (`chol_optimize.md`) and wiring them into the LAPACK paths.

## Host and conditions

- AMD Ryzen Threadripper 3970X (32 cores, 64 threads), NVIDIA A100 80GB PCIe
  (driver 580.178.04), CUDA 13.1, Ubuntu 26.04, gcc with ccache.
- Branch `claude/mpfr-cuda-rgemm` at 328e0dc25 (after 9ddb38e14 "gates
  retuned"); the autotools build in the source tree (configure line in
  `config.log`: `--enable-{dd,td,qd,mpfr}-cuda --with-mpfr-cuda-arch=80,90`).
  The libraries were not rebuilt (09:12-09:20 builds); the benchmark does not
  change any library code.  One source file changed after the MPFR CUDA
  library was built (`mpfr/cuda/Rsolver_bridge_body.h`, a solver gate
  constant, 09:23); no solver is called here.
- **GPU**: used only by this benchmark: `nvidia-smi` before and after every
  GPU series (`gpu-processes.log`) shows no other process; no series is
  flagged (`gpu_busy` column of `kernels.csv` empty everywhere).
- **Host load**: other jobs (not visible in this container) kept the load
  average at about 20 during the whole GPU run (11:06-14:33).  The GPU
  series ran with one host thread each and are timed by the GPU (CUDA
  events, CUPTI timestamps); only the host launches and the status reads
  between the kernels of panels and inverses can be affected by the load.
- **CPU: reference values measured on a loaded host.**  The CPU series were
  measured twice:
  - `raw/cpu-loaded/` (11:26-14:04, all 189 CPU series): concurrently with
    the GPU series and the external load of about 20;
  - `raw/cpu-quiet/` (14:41-15:07, 91 series: the main series of every type
    and all the dd/td/qd series): after the external load had stopped (load
    average 13-33 including the benchmark's own 32 threads,
    `raw/cpu-quiet/loadavg.log`).  The maintainer stopped this re-run when
    the host became busy again; the remaining 98 series (the other MPFR
    series and the complex MPFR series) only have loaded-host values.
  - `kernels.csv` and the graphs use the quiet value when the series has
    one, otherwise the loaded one (column `cpu_host`: `quiet` / `loaded`).
    **The loaded-host CPU values are reference values, about 1.5 times
    slower than on a quiet host** (on the 506 points measured both ways the
    loaded/quiet time ratio has median 1.35, quartiles 1.20 and 1.49; the
    worst small-n points are several times slower).  The graphs say so in
    the legend and the subtitle; in `ratio_*.png` these series are drawn
    with ▲ and labelled "CPU 負荷下", in `summary.png` their CPU bars are
    hatched.  GPU/CPU ratios of these series are therefore too favourable
    to the GPU by about that factor.
  - A first CPU pass pinned to cpus 0-31 with `taskset` is in
    `raw/cpu-pinned/` (not used: under the load the barrier-heavy panels
    waited for descheduled threads, Rgetrf2 up to 240x slower).
- CPU: `libmplapack_<b>_opt`, OMP_NUM_THREADS=32, threads placed by the OS.

## How it was measured

- dd/td/qd GPU: the device kernels of the CUDA libraries compiled into the
  benchmark from their internal headers, matrices already on the device,
  CUDA events.  MPFR GPU: the entry points `<R>_mpfr_cuda`, kernel span
  recorded by CUPTI (start of the first kernel to the end of the last; no
  packing, no transfers).  For the factorization entries the CUPTI span
  agrees with the entry's own factor time (getrf2 512: 22.6 ms against
  23.9 ms; potrf 256: 173 ms against 174 ms).
- GPU: one warm-up run and the best of 3 (one run when the warm-up exceeds
  10 s); CPU: the best of up to 3 (one run when a run exceeds 10 s).
  Same process for all the sizes of a series.
- FLOP counts: GEMM 2mnk (complex 8mnk), SYRK k·n(n+1) (HERK 4k·n(n+1)),
  TRSM Left n·m², Right m·n² (complex ×4), getrf2 the LAWN 41 count of
  getrf for m × 64, potrf2 and lauum n³/3 + n²/2 + n/6, trtri n³/3 + 2n/3;
  laswp: seconds (graphs: 64·n swaps per second).

## Reductions (time budget of 6 hours)

A probe (two sizes per series) predicted 51 GPU hours for the full grid of
the request.  Kept unchanged: the measuring protocol, the types, every real
kernel and form, the 8192 sweep of dd/td/qd GEMM NN and SYRK LN, and the
MPFR GEMM NN sweeps (4096; 1024 and 200 bits 2048).  Reduced (`plan.json`):

1. Fine grid of the request (32 up to 1024, 128 up to 4096, 256 up to 8192)
   only for GEMM NN, the update GEMM (k = 64), SYRK LN (k = n and 64), TRSM
   LLNU, TRSM RLTN with B n×64, and the LU panel; the other series use steps
   64 / 256 (MPFR 512) / 512, all powers of two kept.
2. GPU sizes whose predicted single run exceeds 3 s (15 s for the series of
   item 1, none for the sweeps above) are not measured: the largest n of
   each series is in the tables below (e.g. MPFR complex HERK stops at 256
   or 128, cu_mpfr complex kernels at 64).
3. Complex: GEMM NN and NC, HERK LN (k = n and 64), TRSM LLNU and RLCN with
   B n×64, the shapes of the complex LU and Cholesky (6 of the 15 forms the
   program provides).
4. MPFR GPU: no laswp (not reachable alone through the library), no trtri /
   lauum (no MPFR GPU inverse yet): CPU only.  MPFR GPU potrf2 at n = 256 is
   `Rpotrf_mpfr_cuda`, i.e. the blocked algorithm with nb = 64 (four 64×64
   diagonal blocks and their updates), not the recursive panel.
5. CPU: powers of two 128..4096 (MPFR 2048), potrf2 64 and 256.

## Bit comparison

Every series of every type, n = 64 and 200, GPU against CPU (`checks.csv`,
`kernel_bench.py check`): **354 of 354 comparisons bit-identical**, info and
ipiv included, no difference in the sign of a zero.  (The MPFR series
without a GPU version are not compared.)  The comparison script was checked
on altered dumps: a changed limb and a changed ipiv are reported as
DIFFERENT, a zero of the other sign as zero-sign.

## Results

GFLOPS at n = 2048 (MPFR 1024; potrf2 256); when the GPU was not measured at
that size (time budget), the largest size measured on both sides.  "GPU max"
is the best GPU value of the series.  CPU host: `quiet` or `loaded`
(reference value, about 1.5 times too slow, see above).  laswp: seconds.

#### dd (n = 2048 unless noted)

| kernel | form | n | GPU GFLOPS | CPU GFLOPS | GPU/CPU | GPU max GFLOPS (n) | CPU host | bit check (n = 64, 200) |
|---|---|---:|---:|---:|---:|---:|---|---|
| gemm | NN | 2048 | 320 | 44.3 | 7.23 | 322 (7936) | quiet | identical |
| gemm | NT | 2048 | 320 | 36.4 | 8.8 | 322 (7680) | quiet | identical |
| gemm | TN | 2048 | 321 | 40.9 | 7.85 | 323 (7680) | quiet | identical |
| gemm_update | NN | 2048 | 307 | 34.9 | 8.77 | 315 (4096) | quiet | identical |
| syrk | LN | 2048 | 162 | 38.4 | 4.22 | 165 (7936) | quiet | identical |
| syrk | UT | 2048 | 297 | 35.7 | 8.31 | 305 (6144) | quiet | identical |
| syrk | LN_k64 | 2048 | 158 | 42.3 | 3.73 | 162 (3968) | quiet | identical |
| syrk | UT_k64 | 2048 | 273 | 5.82 | 47 | 292 (4096) | quiet | identical |
| trsm | LLNU | 2048 | 4.36 | 34.4 | 0.127 | 5.54 (1152) | quiet | identical |
| trsm | RLTN | 2048 | 5.68 | 20.4 | 0.279 | 6.2 (2304) | quiet | identical |
| trsm | LUTN | 2048 | 16.5 | 17.7 | 0.934 | 22.4 (3072) | quiet | identical |
| trsm | LUNN | 2048 | 3.88 | 33.4 | 0.116 | 5.39 (1280) | quiet | identical |
| trsm | LLNU_p64 | 2048 | 9.95 | 25.7 | 0.388 | 19.4 (4096) | quiet | identical |
| trsm | RLTN_p64 | 2048 | 10.9 | 14.6 | 0.744 | 21.4 (4096) | quiet | identical |
| trsm | LUTN_p64 | 2048 | 18.4 | 11.1 | 1.66 | 36.8 (4096) | quiet | identical |
| trsm | LUNN_p64 | 2048 | 9.28 | 15.8 | 0.586 | 18.1 (4096) | quiet | identical |
| getrf2 | P64 | 2048 | 2.16 | 1.09 | 1.97 | 3.91 (4096) | quiet | identical |
| potrf2 | L | 256 | 0.439 | 1.12 | 0.39 | 0.439 (256) | quiet | identical |
| laswp | R64 | 2048 | 0.000109 s | 0.000196 s | 1.81 | — | quiet | identical |
| trtri | UN | 1024 | 0.322 | 2.33 | 0.138 | 0.322 (1024) | quiet | identical |
| lauum | L | 2048 | 4.76 | 7.15 | 0.666 | 6.24 (2560) | quiet | identical |
| cgemm | NN | 2048 | 195 | 31.4 | 6.19 | 199 (1536) | loaded | identical |
| cgemm | NC | 2048 | 194 | 31.9 | 6.08 | 199 (1280) | loaded | identical |
| cherk | LN | 1024 | 4.35 | 27 | 0.161 | 4.52 (992) | loaded | identical |
| cherk | LN_k64 | 2048 | 5.75 | 29.5 | 0.195 | 10.4 (3840) | loaded | identical |
| ctrsm | LLNU | 1024 | 7.45 | 23.4 | 0.318 | 9.89 (1536) | loaded | identical |
| ctrsm | RLCN_p64 | 2048 | 17.8 | 13.1 | 1.36 | 33.8 (4096) | loaded | identical |

#### td (n = 2048 unless noted)

| kernel | form | n | GPU GFLOPS | CPU GFLOPS | GPU/CPU | GPU max GFLOPS (n) | CPU host | bit check (n = 64, 200) |
|---|---|---:|---:|---:|---:|---:|---|---|
| gemm | NN | 2048 | 99.4 | 12 | 8.28 | 99.8 (7936) | quiet | identical |
| gemm | NT | 2048 | 99.5 | 9.9 | 10 | 99.8 (4608) | quiet | identical |
| gemm | TN | 2048 | 99.4 | 9.96 | 9.98 | 99.7 (4608) | quiet | identical |
| gemm_update | NN | 2048 | 98.1 | 11.8 | 8.32 | 99.1 (3456) | quiet | identical |
| syrk | LN | 2048 | 48.9 | 11.4 | 4.28 | 49.6 (7680) | quiet | identical |
| syrk | UT | 2048 | 96.5 | 11.1 | 8.7 | 97.8 (6144) | quiet | identical |
| syrk | LN_k64 | 2048 | 48.3 | 11.9 | 4.07 | 49.2 (3968) | quiet | identical |
| syrk | UT_k64 | 2048 | 92.8 | 3.06 | 30.3 | 94.8 (4096) | quiet | identical |
| trsm | LLNU | 2048 | 4.42 | 10.5 | 0.423 | 7.59 (3584) | quiet | identical |
| trsm | RLTN | 1024 | 2.23 | 4.2 | 0.531 | 3.01 (1792) | quiet | identical |
| trsm | LUTN | 2048 | 5.52 | 7.16 | 0.771 | 5.52 (2048) | quiet | identical |
| trsm | LUNN | 2048 | 4.41 | 11.3 | 0.389 | 4.41 (2048) | quiet | identical |
| trsm | LLNU_p64 | 2048 | 6.25 | 8.82 | 0.709 | 12.4 (4096) | quiet | identical |
| trsm | RLTN_p64 | 2048 | 5.18 | 5.47 | 0.948 | 10.3 (4096) | quiet | identical |
| trsm | LUTN_p64 | 2048 | 7.01 | 4.11 | 1.71 | 14 (4096) | quiet | identical |
| trsm | LUNN_p64 | 2048 | 5.77 | 5.25 | 1.1 | 11.5 (4096) | quiet | identical |
| getrf2 | P64 | 2048 | 1.74 | 0.368 | 4.72 | 3.05 (4096) | quiet | identical |
| potrf2 | L | 256 | 0.259 | 0.365 | 0.71 | 0.259 (256) | quiet | identical |
| laswp | R64 | 2048 | 0.000125 s | 0.000236 s | 1.89 | — | quiet | identical |
| trtri | UN | 1024 | 0.162 | 0.701 | 0.231 | 0.162 (1024) | quiet | identical |
| lauum | L | 1024 | 0.661 | 1.46 | 0.451 | 1.44 (1792) | quiet | identical |
| cgemm | NN | 2048 | 58.9 | 8.94 | 6.59 | 59 (2304) | loaded | identical |
| cgemm | NC | 2048 | 58.9 | 9.15 | 6.44 | 58.9 (2304) | loaded | identical |
| cherk | LN | 1024 | 2.17 | 8.52 | 0.254 | 2.17 (1024) | loaded | identical |
| cherk | LN_k64 | 2048 | 3.74 | 8.92 | 0.419 | 5.27 (2816) | loaded | identical |
| ctrsm | LLNU | 1024 | 4.02 | 7.45 | 0.539 | 4.75 (1280) | loaded | identical |
| ctrsm | RLCN_p64 | 2048 | 8.39 | 2.79 | 3.01 | 16.8 (4096) | loaded | identical |

#### qd (n = 2048 unless noted)

| kernel | form | n | GPU GFLOPS | CPU GFLOPS | GPU/CPU | GPU max GFLOPS (n) | CPU host | bit check (n = 64, 200) |
|---|---|---:|---:|---:|---:|---:|---|---|
| gemm | NN | 2048 | 45.7 | 4.82 | 9.49 | 45.9 (6912) | quiet | identical |
| gemm | NT | 2048 | 45.7 | 5.89 | 7.77 | 45.8 (2304) | quiet | identical |
| gemm | TN | 2048 | 45.8 | 5.33 | 8.59 | 45.9 (2304) | quiet | identical |
| gemm_update | NN | 2048 | 45.5 | 4.82 | 9.43 | 45.7 (4096) | quiet | identical |
| syrk | LN | 2048 | 24.4 | 5.63 | 4.32 | 24.7 (8192) | quiet | identical |
| syrk | UT | 2048 | 44.3 | 5.23 | 8.47 | 44.8 (3840) | quiet | identical |
| syrk | LN_k64 | 2048 | 24.2 | 5.75 | 4.2 | 24.6 (3968) | quiet | identical |
| syrk | UT_k64 | 2048 | 43.6 | 1.89 | 23.1 | 44.3 (3840) | quiet | identical |
| trsm | LLNU | 2048 | 3.38 | 5.14 | 0.656 | 4.8 (2944) | quiet | identical |
| trsm | RLTN | 1024 | 1.47 | 2.27 | 0.647 | 1.93 (1536) | quiet | identical |
| trsm | LUTN | 1024 | 2.14 | 1.96 | 1.09 | 3.06 (1792) | quiet | identical |
| trsm | LUNN | 1024 | 1.74 | 4.02 | 0.434 | 2.53 (1536) | quiet | identical |
| trsm | LLNU_p64 | 2048 | 3.52 | 3.83 | 0.92 | 7.01 (4096) | quiet | identical |
| trsm | RLTN_p64 | 2048 | 2.94 | 3.33 | 0.882 | 5.83 (4096) | quiet | identical |
| trsm | LUTN_p64 | 2048 | 4.26 | 2.13 | 2 | 8.48 (4096) | quiet | identical |
| trsm | LUNN_p64 | 2048 | 3.21 | 2.51 | 1.28 | 6.39 (4096) | quiet | identical |
| getrf2 | P64 | 2048 | 1.5 | 0.203 | 7.37 | 2.32 (3712) | quiet | identical |
| potrf2 | L | 256 | 0.158 | 0.191 | 0.827 | 0.158 (256) | quiet | identical |
| laswp | R64 | 2048 | 0.000141 s | 0.000354 s | 2.5 | — | quiet | identical |
| trtri | UN | 512 | 0.0709 | 0.587 | 0.121 | 0.0913 (928) | quiet | identical |
| lauum | L | 1024 | 0.371 | 0.801 | 0.462 | 0.513 (1280) | quiet | identical |
| cgemm | NN | 2048 | 26.4 | 4.38 | 6.03 | 26.4 (2048) | loaded | identical |
| cgemm | NC | 2048 | 26.4 | 4.43 | 5.96 | 26.4 (2048) | loaded | identical |
| cherk | LN | 512 | 0.688 | 3.63 | 0.19 | 1.14 (928) | loaded | identical |
| cherk | LN_k64 | 2048 | 2.42 | 4.45 | 0.544 | 2.58 (2304) | loaded | identical |
| ctrsm | LLNU | 1024 | 2.25 | 3.81 | 0.589 | 2.25 (1024) | loaded | identical |
| ctrsm | RLCN_p64 | 2048 | 4.11 | 2.65 | 1.55 | 8.21 (4096) | loaded | identical |

#### mpfr256 (n = 1024 unless noted)

| kernel | form | n | GPU GFLOPS | CPU GFLOPS | GPU/CPU | GPU max GFLOPS (n) | CPU host | bit check (n = 64, 200) |
|---|---|---:|---:|---:|---:|---:|---|---|
| gemm | NN | 1024 | 15.6 | 0.474 | 33 | 15.9 (2432) | quiet | identical |
| gemm | NT | 1024 | 15.7 | 0.344 | 45.7 | 16 (2560) | loaded | identical |
| gemm | TN | 1024 | 15.1 | 0.338 | 44.6 | 15.2 (2560) | loaded | identical |
| gemm_update | NN | 1024 | 13.3 | 0.392 | 34 | 14.4 (2816) | quiet | identical |
| syrk | LN | 1024 | 8.13 | 0.399 | 20.4 | 9.12 (4096) | quiet | identical |
| syrk | UT | 1024 | 10.4 | 0.461 | 22.7 | 11 (2560) | loaded | identical |
| syrk | LN_k64 | 1024 | 7.4 | 0.363 | 20.4 | 7.96 (4096) | quiet | identical |
| syrk | UT_k64 | 1024 | 8.79 | 0.381 | 23.1 | 9.47 (4096) | loaded | identical |
| trsm | LLNU | 1024 | 0.553 | 0.494 | 1.12 | 0.808 (1536) | quiet | identical |
| trsm | RLTN | 512 | 0.307 | 0.128 | 2.4 | 0.528 (928) | loaded | identical |
| trsm | LUTN | 512 | 0.336 | 0.461 | 0.728 | 0.597 (992) | loaded | identical |
| trsm | LUNN | 512 | 0.298 | 0.328 | 0.909 | 0.503 (928) | loaded | identical |
| trsm | LLNU_p64 | 1024 | 0.583 | 0.32 | 1.82 | 2.04 (3584) | loaded | identical |
| trsm | RLTN_p64 | 1024 | 0.541 | 0.143 | 3.77 | 2.16 (4096) | quiet | identical |
| trsm | LUTN_p64 | 1024 | 0.554 | 0.451 | 1.23 | 1.92 (3584) | loaded | identical |
| trsm | LUNN_p64 | 1024 | 0.505 | 0.303 | 1.67 | 1.74 (3584) | loaded | identical |
| getrf2 | P64 | 1024 | 0.242 | 0.0791 | 3.06 | 0.998 (4096) | quiet | identical |
| potrf2 | L | 256 | 0.0714 | 0.0959 | 0.745 | 0.0714 (256) | loaded | identical |
| laswp | R64 | 1024 | — | 0.00151 s | — | — | loaded | — (no GPU) |
| trtri | UN | 1024 | — | 0.032 | — | — | loaded | — (no GPU) |
| lauum | L | 1024 | — | 0.094 | — | — | loaded | — (no GPU) |
| cgemm | NN | 1024 | 6.75 | 0.331 | 20.4 | 7.03 (288) | loaded | identical |
| cgemm | NC | 1024 | 6.53 | 0.329 | 19.9 | 6.75 (288) | loaded | identical |
| cherk | LN | 256 | 0.0582 | 0.345 | 0.169 | 0.0805 (352) | loaded | identical |
| cherk | LN_k64 | 512 | 0.111 | 0.35 | 0.316 | 0.177 (864) | loaded | identical |
| ctrsm | LLNU | 512 | 0.251 | 0.319 | 0.787 | 0.267 (544) | loaded | identical |
| ctrsm | RLCN_p64 | 1024 | 0.0722 | 0.241 | 0.299 | 0.0722 (1024) | loaded | identical |

#### mpfr512 (n = 1024 unless noted)

| kernel | form | n | GPU GFLOPS | CPU GFLOPS | GPU/CPU | GPU max GFLOPS (n) | CPU host | bit check (n = 64, 200) |
|---|---|---:|---:|---:|---:|---:|---|---|
| gemm | NN | 1024 | 6.11 | 0.374 | 16.3 | 6.33 (4096) | quiet | identical |
| gemm | NT | 1024 | 6.1 | 0.279 | 21.9 | 6.19 (1536) | loaded | identical |
| gemm | TN | 1024 | 5.71 | 0.27 | 21.2 | 5.84 (1536) | loaded | identical |
| gemm_update | NN | 1024 | 5.7 | 0.316 | 18.1 | 5.86 (3840) | quiet | identical |
| syrk | LN | 1024 | 2.27 | 0.317 | 7.16 | 2.54 (3072) | quiet | identical |
| syrk | UT | 1024 | 3.43 | 0.35 | 9.8 | 3.74 (2048) | loaded | identical |
| syrk | LN_k64 | 1024 | 2.42 | 0.288 | 8.43 | 2.62 (4096) | quiet | identical |
| syrk | UT_k64 | 1024 | 2.9 | 0.308 | 9.44 | 3.08 (3584) | loaded | identical |
| trsm | LLNU | 1024 | 0.286 | 0.398 | 0.718 | 0.357 (1280) | quiet | identical |
| trsm | RLTN | 512 | 0.16 | 0.125 | 1.28 | 0.224 (736) | loaded | identical |
| trsm | LUTN | 512 | 0.17 | 0.361 | 0.471 | 0.236 (736) | loaded | identical |
| trsm | LUNN | 512 | 0.147 | 0.282 | 0.52 | 0.206 (736) | loaded | identical |
| trsm | LLNU_p64 | 1024 | 0.302 | 0.288 | 1.05 | 0.889 (3072) | loaded | identical |
| trsm | RLTN_p64 | 1024 | 0.241 | 0.134 | 1.8 | 0.885 (4096) | quiet | identical |
| trsm | LUTN_p64 | 1024 | 0.25 | 0.343 | 0.729 | 0.586 (2560) | loaded | identical |
| trsm | LUNN_p64 | 1024 | 0.227 | 0.275 | 0.828 | 0.552 (2560) | loaded | identical |
| getrf2 | P64 | 1024 | 0.173 | 0.0666 | 2.6 | 0.581 (4096) | quiet | identical |
| potrf2 | L | 256 | 0.0324 | 0.0827 | 0.392 | 0.0324 (256) | loaded | identical |
| laswp | R64 | 1024 | — | 0.00124 s | — | — | loaded | — (no GPU) |
| trtri | UN | 1024 | — | 0.0255 | — | — | loaded | — (no GPU) |
| lauum | L | 1024 | — | 0.0736 | — | — | loaded | — (no GPU) |
| cgemm | NN | 512 | 1.53 | 0.265 | 5.76 | 1.57 (288) | loaded | identical |
| cgemm | NC | 512 | 1.5 | 0.258 | 5.81 | 1.54 (288) | loaded | identical |
| cherk | LN | 256 | 0.0246 | 0.276 | 0.0889 | 0.0246 (256) | loaded | identical |
| cherk | LN_k64 | 512 | 0.0478 | 0.268 | 0.178 | 0.0563 (608) | loaded | identical |
| ctrsm | LLNU | 256 | 0.0505 | 0.246 | 0.205 | 0.0693 (352) | loaded | identical |
| ctrsm | RLCN_p64 | 1024 | 0.067 | 0.201 | 0.333 | 0.067 (1024) | loaded | identical |

#### mpfr1024 (n = 1024 unless noted)

| kernel | form | n | GPU GFLOPS | CPU GFLOPS | GPU/CPU | GPU max GFLOPS (n) | CPU host | bit check (n = 64, 200) |
|---|---|---:|---:|---:|---:|---:|---|---|
| gemm | NN | 1024 | 0.973 | 0.211 | 4.61 | 1.02 (2048) | quiet | identical |
| gemm | NT | 1024 | 0.973 | 0.154 | 6.31 | 0.973 (1024) | loaded | identical |
| gemm | TN | 1024 | 0.939 | 0.145 | 6.47 | 0.939 (1024) | loaded | identical |
| gemm_update | NN | 1024 | 1.05 | 0.181 | 5.82 | 1.06 (1920) | quiet | identical |
| syrk | LN | 1024 | 0.434 | 0.185 | 2.35 | 0.471 (1280) | quiet | identical |
| syrk | UT | 1024 | 0.84 | 0.195 | 4.31 | 0.84 (1024) | loaded | identical |
| syrk | LN_k64 | 1024 | 0.525 | 0.169 | 3.1 | 0.556 (512) | quiet | identical |
| syrk | UT_k64 | 1024 | 0.709 | 0.176 | 4.04 | 0.732 (2048) | loaded | identical |
| trsm | LLNU | 512 | 0.0584 | 0.217 | 0.269 | 0.111 (928) | quiet | identical |
| trsm | RLTN | 512 | 0.0542 | 0.0882 | 0.614 | 0.0576 (544) | loaded | identical |
| trsm | LUTN | 512 | 0.0608 | 0.183 | 0.332 | 0.0646 (544) | loaded | identical |
| trsm | LUNN | 512 | 0.0544 | 0.162 | 0.336 | 0.0577 (544) | loaded | identical |
| trsm | LLNU_p64 | 1024 | 0.122 | 0.161 | 0.762 | 0.245 (2048) | loaded | identical |
| trsm | RLTN_p64 | 1024 | 0.0769 | 0.102 | 0.752 | 0.153 (2048) | quiet | identical |
| trsm | LUTN_p64 | 1024 | 0.0821 | 0.182 | 0.451 | 0.112 (1536) | loaded | identical |
| trsm | LUNN_p64 | 1024 | 0.0781 | 0.147 | 0.533 | 0.111 (1536) | loaded | identical |
| getrf2 | P64 | 1024 | 0.0652 | 0.0494 | 1.32 | 0.123 (2048) | quiet | identical |
| potrf2 | L | 256 | 0.0084 | 0.0526 | 0.16 | 0.0084 (256) | loaded | identical |
| laswp | R64 | 1024 | — | 0.00134 s | — | — | loaded | — (no GPU) |
| trtri | UN | 1024 | — | 0.0155 | — | — | loaded | — (no GPU) |
| lauum | L | 1024 | — | 0.0424 | — | — | loaded | — (no GPU) |
| cgemm | NN | 512 | 0.601 | 0.153 | 3.92 | 0.656 (288) | loaded | identical |
| cgemm | NC | 512 | 0.632 | 0.149 | 4.23 | 0.64 (288) | loaded | identical |
| cherk | LN | 128 | 0.00505 | 0.151 | 0.0334 | 0.00631 (160) | loaded | identical |
| cherk | LN_k64 | 256 | 0.0099 | 0.161 | 0.0613 | 0.0159 (416) | loaded | identical |
| ctrsm | LLNU | 128 | 0.0104 | 0.147 | 0.0706 | 0.0179 (224) | loaded | identical |
| ctrsm | RLCN_p64 | 1024 | 0.0394 | 0.132 | 0.298 | 0.0394 (1024) | loaded | identical |

#### mpfr200rt (n = 1024 unless noted)

| kernel | form | n | GPU GFLOPS | CPU GFLOPS | GPU/CPU | GPU max GFLOPS (n) | CPU host | bit check (n = 64, 200) |
|---|---|---:|---:|---:|---:|---:|---|---|
| gemm | NN | 1024 | 0.545 | 0.424 | 1.29 | 0.563 (2048) | quiet | identical |
| gemm | NT | 512 | 0.526 | 0.342 | 1.54 | 0.538 (864) | loaded | identical |
| gemm | TN | 512 | 0.518 | 0.329 | 1.57 | 0.524 (864) | loaded | identical |
| gemm_update | NN | 1024 | 0.482 | 0.481 | 1 | 0.484 (2048) | quiet | identical |
| syrk | LN | 1024 | 0.192 | 0.388 | 0.493 | 0.332 (800) | quiet | identical |
| syrk | UT | 512 | 0.279 | 0.505 | 0.553 | 0.488 (800) | loaded | identical |
| syrk | LN_k64 | 1024 | 0.178 | 0.362 | 0.493 | 0.307 (1920) | quiet | identical |
| syrk | UT_k64 | 1024 | 0.253 | 0.401 | 0.63 | 0.428 (800) | loaded | identical |
| trsm | LLNU | 512 | 0.0354 | 0.487 | 0.0726 | 0.055 (800) | quiet | identical |
| trsm | RLTN | 256 | 0.0184 | 0.15 | 0.122 | 0.0301 (416) | loaded | identical |
| trsm | LUTN | 512 | 0.0457 | 0.47 | 0.0973 | 0.0457 (512) | loaded | identical |
| trsm | LUNN | 256 | 0.0174 | 0.324 | 0.0537 | 0.0284 (416) | loaded | identical |
| trsm | LLNU_p64 | 1024 | 0.0686 | 0.344 | 0.199 | 0.103 (1536) | loaded | identical |
| trsm | RLTN_p64 | 1024 | 0.0679 | 0.148 | 0.459 | 0.135 (2048) | quiet | identical |
| trsm | LUTN_p64 | 1024 | 0.0804 | 0.452 | 0.178 | 0.12 (1536) | loaded | identical |
| trsm | LUNN_p64 | 1024 | 0.0652 | 0.349 | 0.187 | 0.0979 (1536) | loaded | identical |
| getrf2 | P64 | 1024 | 0.046 | 0.0789 | 0.583 | 0.0765 (2048) | quiet | identical |
| potrf2 | L | 256 | 0.0132 | 0.0973 | 0.136 | 0.0132 (256) | loaded | identical |
| laswp | R64 | 1024 | — | 0.000906 s | — | — | loaded | — (no GPU) |
| trtri | UN | 1024 | — | 0.0317 | — | — | loaded | — (no GPU) |
| lauum | L | 1024 | — | 0.094 | — | — | loaded | — (no GPU) |
| cgemm | NN | 1024 | — | 0.36 | — | 0.00149 (64) | loaded | identical |
| cgemm | NC | 1024 | — | 0.391 | — | 0.00149 (64) | loaded | identical |
| cherk | LN | 1024 | — | 0.394 | — | 0.000185 (64) | loaded | identical |
| cherk | LN_k64 | 1024 | — | 0.408 | — | 0.000185 (64) | loaded | identical |
| ctrsm | LLNU | 1024 | — | 0.386 | — | 0.00031 (64) | loaded | identical |
| ctrsm | RLCN_p64 | 1024 | — | 0.3 | — | 0.0003 (64) | loaded | identical |


### In short (GPU kernel against CPU 32 cores; "loaded" CPU values favour the GPU)

- **GEMM** (the update of LU and Cholesky): dd 320 GFLOPS (7.2x the CPU at
  n = 2048), td 99 (8.3x), qd 46 (9.5x), MPFR 256 / 512 / 1024 bits 15.6 /
  6.1 / 0.97 (33x / 16x / 4.6x at n = 1024); MPFR 200 bits (cu_mpfr) 0.55
  (1.3x).  The update shape (k = 64) is as fast as the square one from
  n ≈ 1000.  Complex GEMM: dd 195, td 59, qd 26 GFLOPS (6x).
- **SYRK** Lower NoTrans: dd 162 (4.2x), td 49 (4.3x), qd 24 (4.3x), MPFR 256 /
  512 / 1024: 8.1 / 2.3 / 0.43 (20x / 7x / 2.3x); cu_mpfr 0.19 (0.5x, slower
  than the CPU).  Upper Trans is about twice as fast as Lower NoTrans on
  dd/td/qd (297 against 162 GFLOPS for dd).
- **TRSM**: slower than the CPU in most shapes: LLNU (LU) dd 4.4 GFLOPS
  against 34 (0.13x), td 0.42x, qd 0.66x, MPFR 256 1.1x, 512 0.72x, 1024
  0.27x, cu_mpfr 0.07x.  RLTN with B n×64 (Cholesky L panel): 0.74x dd,
  0.95x td, 0.88x qd, 3.8x / 1.8x / 0.75x MPFR 256 / 512 / 1024 (the MPFR
  CPU values of the non-main TRSM forms are loaded-host values).  On
  dd/td/qd only LUTN with B 64×n (1.7-2.0x) clearly beats the CPU.
- **Panels**: LU panel (m×64) on the GPU 2.0x (dd), 4.7x (td), 7.4x (qd),
  3.1x / 2.6x / 1.3x (MPFR 256 / 512 / 1024), 0.58x (cu_mpfr); its speed
  grows linearly with m (latency-bound: 64 pivot steps, each with launches
  and a status read).  Cholesky diagonal block (potrf2, 256×256): slower than
  the CPU for every type (0.39x dd, 0.71x td, 0.83x qd, 0.39x MPFR 512,
  0.16x MPFR 1024, 0.14x cu_mpfr).
- **laswp** (64 swaps on n×n): GPU 1.8x (dd) to 3.6x (qd) at n = 2048.
- **trtri / lauum** (dd/td/qd GPU inverse engine): trtri 0.14x (dd), 0.24x
  (td), 0.16x (qd at n = 512); lauum 0.53-0.67x.  MPFR: no GPU version.
- **Complex HERK and TRSM**: much slower than the CPU: dd Cherk 4.4 GFLOPS
  against 27 (0.16x, loaded CPU), MPFR 512 Cherk 0.025 GFLOPS (0.09x), MPFR
  1024 Cherk 0.005 (0.03x); cu_mpfr complex kernels about 0.0002-0.0015
  GFLOPS (Cgemm n = 64: 1.4 s).

### Slow GPU kernels (candidates for optimization)

1. **TRSM** (real, all types): one thread per column / row of B
   (`blas_device` lanes): with n columns there are only n threads, each
   doing O(n²) sequential work.  dd 4-6 GFLOPS against 320 for GEMM;
   blocked TRSM (diagonal blocks + GEMM updates) is the obvious remedy.
2. **Cholesky diagonal block potrf2** and **trtri / lauum**: latency-bound
   recursion with many small launches and status reads; slower than the
   CPU everywhere.
3. **Complex HERK / TRSM** (`job_Cherk`, `job_Ctrsm`: one thread per column)
   and, for MPFR, the complex kernels in general; **cu_mpfr (run-time
   precision) complex kernels** are 1000x slower than the real ones.
4. **SYRK Lower NoTrans** is half the speed of Upper Trans (dd 162 against
   297 GFLOPS): the tiled kernel's NoTrans path (alpha multiply and zero test
   per element) — the shape of the lower Cholesky update.
5. **MPFR 1024 bits and cu_mpfr**: GEMM 0.97 and 0.55 GFLOPS only; SYRK on
   cu_mpfr is slower than the CPU.
6. **LU panel**: linear in m up to 4096, i.e. dominated by per-column
   latency (launches + pivot status read), not by arithmetic.

## Graphs

Drawn with `benchmark/plot_kernel_bench.py` from `kernels.csv`, through
`plot_sweep.draw()` (the code of `../rgemm/mpfr512_sweep10k.png`; the
refactored `plot_sweep.py` redraws that PNG byte for byte, light and dark):

- `<type>/<kernel>_<form>.png` (182 graphs; type = dd, td, qd, mpfr256,
  mpfr512, mpfr1024, mpfr200rt): GFLOPS against n, GPU points joined, the
  powers of two marked and annotated, CPU powers of two dashed (legend: quiet
  host, or the loaded-host reference values); MPFR laswp / trtri / lauum are
  CPU only.  potrf2 (two sizes) is only in the tables and `summary.png`.
- `ratio_<type>.png`: CPU time / GPU time at the powers of two, log scale,
  real and complex panels, one colour per kernel, the 1x line; ▲ =
  loaded-host CPU.
- `summary.png`: GFLOPS at n = 2048 (MPFR 1024), CPU and GPU bars, every type
  and kernel; hatched = loaded-host CPU.

## Files

| File | Contents |
|---|---|
| `kernels.csv` | every point: type, kernel, form, n, method (gpu/cpu), seconds, GFLOPS, runs, bit check, gpu_busy, cpu_host |
| `checks.csv` | bit comparison at n = 64 and 200 (dumps not kept: 1.6 GB, `kernel_bench.py check` rewrites them) |
| `raw/<type>_gpu.{csv,log}` | GPU output and stderr / commands / exit codes |
| `raw/cpu-quiet/`, `raw/cpu-loaded/`, `raw/cpu-pinned/` | the three CPU passes (see above); `cpu-quiet/loadavg.log` the load during the quiet pass |
| `raw/probe_*`, `probe.csv`, `plan.json` | the probe (two sizes per series) and the plan actually run (sizes of every series, estimates, note) |
| `raw/check_*.csv` | timings printed during the bit check (one run, not used) |
| `raw/run-progress.log` | the runner's progress (series, points, seconds, estimate) |
| `gpu-processes.log` | `nvidia-smi` before / after every GPU series |
| `raw-host-at-start.txt`, `raw-host-at-end.txt`, `raw-nvidia-smi-start.txt` | load and GPU state at the start and the end |

## Rebuild, rerun, redraw

    # autotools tree configured with --enable-{dd,td,qd,mpfr}-cuda
    make -C benchmark Rkernel_bench.{dd,td,qd}_cuda_kernel Rkernel_bench.{dd,td,qd,mpfr}_opt Rkernel_bench.mpfr_cuda_kernel
    # CMake: -DMPLAPACK_BUILD_BENCHMARKS=ON with the CUDA libraries, programs in <build>/benchmark
    D=benchmark/results/2026/x86_64-pc-linux-gnu-Ryzen-Threadripper-3970X-TeslaA100/kernel-bench
    python3 benchmark/kernel_bench.py probe --bindir benchmark --out $D
    python3 benchmark/kernel_bench.py plan  --out $D --max-point 3 --max-point-fine 15 --coarse
    python3 benchmark/kernel_bench.py run   --bindir benchmark --out $D --hours 5.5
    python3 benchmark/kernel_bench.py check --bindir benchmark --out $D
    python3 benchmark/kernel_bench.py merge --out $D --cpu-sources quiet=raw/cpu-quiet,loaded=raw/cpu-loaded
    python3 benchmark/kernel_bench.py table --out $D --types dd,td,qd,mpfr256,mpfr512,mpfr1024,mpfr200rt
    PYTHONPATH=/home/docker/plotlib python3 benchmark/plot_kernel_bench.py $D/kernels.csv $D --gpu-name "NVIDIA A100 80GB PCIe"

(`run` writes the CPU results to `raw/<type>_cpu.csv`; they were moved to
`raw/cpu-loaded/` afterwards, and the quiet pass ran with `--out` in a
separate directory, moved to `raw/cpu-quiet/`.)

## Session log (summary of the local WORKLOG / REPORT)

- 10:10 start; prerequisites present (9ddb38e14), no uncommitted changes of
  other sessions (only build artifacts untracked).
- 10:20-10:50 programs written; first checks: dd 17 series and MPFR 512 10
  series bit-identical at n = 64, 200; CUPTI span agrees with the entries'
  own factor times.
- 10:45-11:05 probe; the full grid was estimated at 51 GPU hours; reductions
  above.  11:06 run started (GPU and CPU streams concurrently); the time
  model was recalibrated at 11:19 (first estimates 2-4x low) and the plan
  revised (runner restarted at a series boundary, nothing lost).
- 11:21 CPU pinning dropped (see above), CPU stream restarted at 11:26.
- 14:04 CPU stream done; 14:33 GPU stream done (no series omitted, no other
  GPU process); 14:34-14:39 bit check (354/354 identical).
- 14:41-15:07 quiet-host CPU re-run (91 of 189 series), stopped by the
  maintainer when the host became busy again.

## Not done / follow-up

- **CPU on a quiet host for the remaining 98 series** (the non-main MPFR
  series and all complex MPFR series): `kernel_bench.py run --only cpu
  --plan raw/cpu-quiet/plan-all.json --out <dir>` resumes after the series
  already in `<dir>/raw/`, then `merge` with the new directory first.
- GPU sizes above the 3 s / 15 s caps; the 9 complex forms not measured
  (NC/CN of GEMM besides NC, UC HERK, LUCN/LUNN/... TRSM).
- MPFR GPU laswp, trtri, lauum (no GPU version reachable); add them to
  `Rkernel_bench.cpp` (`KB_CALL`) and drop them from `MPFR_NO_GPU` in
  `kernel_bench.py` when the MPFR GPU inverse exists.
- MPFR GPU potrf2 at n = 256 is the blocked entry (nb = 64), not a recursive
  panel of 256.
