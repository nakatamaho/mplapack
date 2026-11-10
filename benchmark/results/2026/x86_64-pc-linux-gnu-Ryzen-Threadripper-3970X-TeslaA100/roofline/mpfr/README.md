# MPFR (cu_freal) roofline in two levels (A100 80GB PCIe, November 6, 2026)

The dd/td/qd roofline (../README.md) has one arithmetic bound because a QD_BF
operation is a fixed instruction sequence (dd: 20 FP64 instructions per add,
9 per mul). An MPFR fixed-precision operation is not: its cost depends on the
implementation (cu_fp::cu_freal of mpc_cuda) and on the operands (sign, exponent
gap). So the MPFR roofline has two levels:

- **Level 1, hardware** (any implementation of PB-bit floating point on this GPU):
  a multiplication needs at least the N^2 limb products of a full schoolbook product
  (N = PB/64; about N(N+1)/2 + N with a truncated product), an addition N funnel
  shifts and N adds with carry. Rates measured on the A100 with
  `benchmark/mpfr_path_bench.cu` (`paths.csv`, pb = 0): one 64x64->128 limb
  product with carries (mad.lo.cc/madc.hi/add.cc/addc) 9.45e11/s, add with carry
  4.25e12/s, 64-bit funnel shift 1.50e12/s.
- **Level 2, cu_freal as implemented**: the measured throughput of cu_fmul and of
  cu_fadd on each of its paths (`paths.csv`), weighted by how often the kernel's
  additions take each path (`gemm_path_mix.txt`, from `benchmark/mpfr_gemm_paths.py`:
  the reference GEMM accumulation on uniform [-1, 1) data, as the kernel benchmark):
  about 50% same-sign (cu_addmag) and 50% opposite-sign (cu_submag), exponent gaps
  almost all 1-63 bits, no cancellation beyond one limb.

GEMM (n^3 multiply-add pairs, 2 flops each) against both levels:

| PB | level 1 GFLOPS | level 2 GFLOPS | GEMM measured | GEMM / level 2 | level 2 / level 1 | GEMM / level 1 |
|---:|---:|---:|---:|---:|---:|---:|
| 256 | 97.4 | 29.1 | 15.8 | 54% | 30% | 16% |
| 512 | 26.7 | 11.7 | 6.31 | 54% | 44% | 24% |
| 1024 | 7.0 | 3.67 | 1.02 | 28% | 52% | 15% |

Per operation, cu_freal against level 1:

| PB | cu_fmul / full-product bound | cu_fadd (GEMM mix) / add bound |
|---:|---:|---:|
| 256 | 50% | 10% |
| 512 | 61% | 12% |
| 1024 | 62% | 13% |

## Reading

- **Kernel (MPLAPACK side)**: GEMM reaches 54% of what cu_freal allows at 256 and
  512 bits and 28% at 1024 bits. The earlier single-number MPFR roofline (../peaks.csv,
  madd 9.1e9/s at 512 bits, GEMM "34%") used operands converted from double
  (53 significant bits) and only same-sign additions; with full-width random
  significands and the real sign mix, cu_freal itself allows 5.9e9 pairs/s, so
  the kernel is closer to its implementation bound than 34% suggested.
- **cu_freal (mpc_cuda side)**: the multiplication is at 50-62% of the full
  schoolbook product; a truncated product (only the limbs that affect rounding)
  would lower the level 1 bound itself by up to ~45% at 1024 bits. The addition
  is the weak part, at 10-13% of its bound: `finalize` and `cu_submag` index
  their limb arrays at run-time positions, so the arrays live in local memory
  (stack frames 112/208/688 bytes at 256/512/1024 bits in this benchmark), and
  they test bits one at a time. An addition costs as much as a multiplication
  at 256 bits although it needs ~5x fewer instructions.
- The add path barely depends on the exponent gap (within ~5%); the sign matters
  more (opposite signs ~10-15% slower), heavy cancellation more (~30% at 256 bits).

## Use in the roofline of all kernels (November 6, 2026)

`../peaks.csv` now carries the level 2 MPFR peaks (madd = 1/(1/mul + 1/add),
add weighted by the GEMM path mix; `benchmark/mpfr_level_peaks.py`), and
`peaks-level1.csv` the level 1 ones; `div` stays the peak_bench value in both
(not modelled). Recomputed with them: `../roofline_summary.md`,
`../goal11/` (level 2), `level1/`, `level1-goal11/` (level 1). The path mix is
GEMM's for every kernel; TRSM, the panels and the inverses have their own mix
(not measured), so their level 2 % is approximate.

The goal 10 MPFR kernel times so far (items 1-4, from bench_local/results/gpu-opt;
`goal10/kernels-goal10-after.csv`) against the old, level 2 and level 1 bounds:
`goal10/table.md`. For example MPFR 512 n = 1024: trsm LLNU / RLTN
38.9 / 38.5% (old) -> 60.3 / 59.9% (level 2), 26.6 / 26.4% (level 1); syrk UT
24.0 -> 37.3%; Lower potrf2 (blocked, after item 4) 4.5 -> 6.9%.

## Reproduce

    nvcc -O3 -std=c++17 -fmad=false -arch=sm_80 -U_GNU_SOURCE -D_XOPEN_SOURCE=700 \
         -Iexternal/i/MPC_CUDA/include benchmark/mpfr_path_bench.cu -o mpfr_path_bench
    ./mpfr_path_bench > paths.csv
    python3 benchmark/mpfr_gemm_paths.py 1024 300 > gemm_path_mix.txt
    python3 benchmark/mpfr_level_peaks.py ../peaks-double-operands-20261104.csv paths.csv \
            gemm_path_mix.txt ../peaks.csv peaks-level1.csv
    python3 benchmark/roofline.py ../peaks.csv ../../kernel-bench/kernels.csv ../../gemm-sweep ..
    python3 benchmark/roofline.py peaks-level1.csv ../../kernel-bench/kernels.csv ../../gemm-sweep level1

## v0.0.2-dev update (goal 12, November 10, 2026)

The bundled mpc_cuda moved to upstream v0.0.2-dev (0f219b0, plus the rebased
correctly-rounded div/sqrt and canonical-zero patches). Regenerated from the
new paths: `paths-v0.0.2-dev.csv`, `peaks-v0.0.2-dev.csv`,
`peaks-level1-v0.0.2-dev.csv` (the old `paths.csv` / `peaks.csv` /
`peaks-level1.csv` remain the abd9407 "before" reference). Level 2 (old -> new):
madd 1.4533e10 -> 3.2442e10 (256), 5.8724e9 -> 1.4133e10 (512),
1.8367e9 -> 4.3469e9 (1024) pairs/s (2.23x/2.41x/2.37x); the level 1 bounds are
unchanged. Kernel-only GEMM NN sweep (best of 3): 256 n2048 15.87 -> 37.31
GFLOPS, 512 n2048 6.25 -> 16.67, 1024 n2048 1.32 -> 4.21 (3.19x). MPFR 512 GPU
factor entries (whole entry, warmup 1/best of 3, one neutral driver binary,
only the library swapped): real getrf 4096 9.26 -> 4.92 s, real potrf U 4096
9.46 -> 4.86 s, complex getrf 4096 160.2 -> 107.3 s, complex potrf L 4096
37.2 -> 19.9 s, complex potrf U 4096 58.0 -> 32.3 s; all 24 measured points
improved (raw CSVs in bench_local/results/goal12/{campaign,factor}/).
-v0.0.2-dev stack frames of the arithmetic paths (Xptxas -v): 416/224/128 bytes
at 1024/512/256 (abd9407: 688/208/112). The fixed-precision complex factor
stack budget in mpblas/optimized/mpfr/cuda/mpfr_cuda_stack.h was raised for
v0.0.2-dev's deeper cu_mpc_div frames (16/32 KiB; see the goal-12 commit).
