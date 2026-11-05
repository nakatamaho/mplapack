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

## Reproduce

    nvcc -O3 -std=c++17 -fmad=false -arch=sm_80 -U_GNU_SOURCE -D_XOPEN_SOURCE=700 \
         -Iexternal/i/MPC_CUDA/include benchmark/mpfr_path_bench.cu -o mpfr_path_bench
    ./mpfr_path_bench > paths.csv
    python3 benchmark/mpfr_gemm_paths.py 1024 300 > gemm_path_mix.txt
