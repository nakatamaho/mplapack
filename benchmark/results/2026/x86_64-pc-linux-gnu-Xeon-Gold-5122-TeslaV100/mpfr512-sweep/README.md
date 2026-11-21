# MPFR 512 Rgemm and Rgetrf/Rpotrf on the V100 (2026-11-20)

Tesla V100-PCIE-32GB (sm_70, 1380 MHz, driver 580.196.07) + Xeon Gold 5122,
CUDA 12.6 from NVIDIA HPC SDK 24.11, GCC 11.  Commit and local diff in
`commit.txt` (master 282a9955b plus the CUPTI record fix of
`benchmark/Rkernel_bench.cpp` and the GCC 11 target_clones fix of
`mpblas/optimized/common/Rgemm_simd_common.h`, committed with this
directory).  GPU otherwise idle (`gpu-processes-at-*.txt`, `run.log`).
Reproduce with `./run.sh` (`GEMM_ONLY=1` or `FACTOR_ONLY=1` runs one half).

Same method and sizes as the A100 run
`../../x86_64-pc-linux-gnu-Ryzen-Threadripper-3970X-TeslaA100/mpfr512-sweep-2`,
which the dashed lines and the A100 columns show.

## Build on this machine

The V100 needs sm_70 instead of the sm_80/sm_90 the build assumes, and the
HPC SDK keeps cuBLAS and CUPTI outside the CUDA directory, so make was given:

    ML=/opt/nvidia/hpc_sdk/Linux_x86_64/24.11/math_libs/12.6/targets/x86_64-linux
    CU=/opt/nvidia/hpc_sdk/Linux_x86_64/24.11/cuda/12.6
    CP=$CU/extras/CUPTI
    CL="-L$CU/lib64 -L$CP/lib64 -L$ML/lib -Wl,-rpath,$ML/lib -Wl,-rpath,$CU/lib64 -Wl,-rpath,$CP/lib64"
    ./configure ... --enable-mpfr-cuda=yes --with-mpfr-cuda-arch=70 --with-cudatoolkithome=$CU
    make CPPFLAGS="-I$CP/include -I$ML/include" CUDA_CFLAGS="-I$CU/include -I$ML/include" \
         CUDA_LIBS="$CL" CUDA_LDFLAGS="$CL" \
         LDFLAGS="-L$ML/lib -Wl,-rpath,$ML/lib -Wl,-rpath,$CU/lib64 -Wl,-rpath,$CP/lib64" \
         AMPERE_FLAGS="-gencode arch=compute_70,code=sm_70 -m64 -O3" HOPPER_FLAGS="" \
         KERNEL_BENCH_CUDA_ARCH_FLAGS="-gencode arch=compute_70,code=sm_70"

`-I$CP/include` matters: /usr/include has the CUPTI of CUDA 11.5 (API 16,
no CUpti_ActivityKernel9).  The bundled mpc_cuda builds and passes its tests
with CUDA 12.6 for sm_70.

## Rgemm NN/NT/TN/TT (kernel only, CUDA events)

`benchmark/run_gemm_sweep.sh`: n = multiples of 7 and 31 and powers of 2 up
to 4096, best of 3 after 1 warmup for n <= 2000, one run above.  About 82 min
per form (A100: 47 min).

![GEMM V100 vs A100](gemm_mpfr512_v100_a100.png)

GFLOPS, V100 / A100:

| n | NN | NT | TN | TT |
|---|---|---|---|---|
| 512 | 6.29 / 13.04 | 6.34 / 13.11 | 6.50 / 15.64 | 6.58 / 15.64 |
| 1024 | 7.36 / 16.37 | 7.39 / 16.35 | 7.16 / 16.32 | 7.15 / 16.32 |
| 2048 | 7.67 / 16.67 | 7.63 / 16.67 | 7.51 / 16.80 | 7.51 / 16.81 |
| 4096 | 7.92 / 16.94 | 7.92 / 16.94 | 7.85 / 17.16 | 7.85 / 17.17 |

The V100 reaches 0.46 of the A100; the four forms are within 1% of each
other at large n.

## Rgetrf / Rpotrf L / Rpotrf U (GPU entry, whole call incl. packing and transfers)

`Rfactor.mpfr_cuda_total --entry --warmup 1 --max 4096` with
`MPLAPACK_MPFR_LAPACK_CUDA_MIN=1` so every n runs on the GPU; n = powers of
2 and multiples of 13.  36, 30 and 38 min.

![Factorizations V100 vs A100](factor_mpfr512_v100_a100.png)

GFLOPS (getrf 2n³/3, potrf n³/3), V100 / A100 (ratio):

| n | Rgetrf | Rpotrf L | Rpotrf U |
|---|---|---|---|
| 512 | 0.89 / 0.98 (0.91) | 0.32 / 0.38 (0.84) | 0.28 / 0.35 (0.80) |
| 1024 | 2.03 / 2.47 (0.82) | 0.93 / 1.14 (0.82) | 0.75 / 0.97 (0.77) |
| 2048 | 3.91 / 5.63 (0.69) | 2.32 / 3.10 (0.75) | 1.69 / 2.44 (0.69) |
| 4096 | 5.82 / 9.58 (0.61) | 4.21 / 6.37 (0.66) | 3.05 / 4.87 (0.63) |

The gap to the A100 is smaller than for Rgemm and grows with n, as the
GEMM-bound share of the factorization grows; the parts outside the GEMM
kernel (panels, packing, transfers, the CPU) were not timed separately.

Plots: `plot_gemm.py`, `plot_factor.py` (matplotlib).
