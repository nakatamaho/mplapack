# Roofline of the GPU kernels (A100 80GB PCIe, November 4, 2026)

How close each GPU kernel is to the limit of the hardware for its type.

## Method

1. **Peaks** (`peaks.csv`, from `benchmark/peak_bench.cu`, measured 2026-11-04
   23:10 with no other process on the GPU): throughput of the multiple-precision
   operations with the kernels' own arithmetic (dd/td/qd: libQD3 QD_BF algorithms
   of `mplapack_gemm::{dd,td,qd}_ops`, `-fmad=false`; MPFR: `cu_fp::cu_freal<PB>`,
   division `cu_fdiv_rn`), each the best over launch shapes and independent chains:
   `madd` (c = c + c·b), `add`, `mul`, `div`, and the device-memory bandwidth
   (read + write, 1.40 TB/s). `peaks-shared-gpu-reference.csv` is an earlier
   run while another job shared the GPU (lower; kept for reference only).
2. **Counts**: multiplications and additions separately, from LAPACK Working
   Note 41 (https://netlib.org/lapack/lawnspdf/lawn41.pdf); a complex
   multiplication = 4 real multiplications + 2 additions.
3. **Bound**: t_bound = max(t_arith, t_mem), t_arith = pairs/madd +
   (muls-pairs)/mul + (adds-pairs)/add + divs/div with pairs = min(muls, adds),
   t_mem = (elements read once + written once) × element size / bandwidth.
   **% of roofline = t_bound / t_measured.** "% of GEMM plateau" relates the
   kernel's operation rate to the GEMM NN plateau of the same type
   (../gemm-sweep/).
4. Measured times: ../kernel-bench/kernels.csv (GPU kernel only, data on the
   device, CUDA events) and ../gemm-sweep/. The runtime-precision MPFR
   (cu_mpfr, mpfr200rt) has no arithmetic peak here, so no roofline value.

Check: the dd/td/qd GEMM reaches 99-101% of the madd peak, i.e. the madd
microbenchmark and GEMM agree. The dd peaks also agree with the instruction
counts: a QD_BF dd addition is about 20 FP64 instructions (A100: 4.85e12/s,
2.4e11 additions/s; measured 2.3e11), a multiplication about 9 with FMA (5.4e11;
measured 4.6e11).

## Result (% of roofline, n = 2048; MPFR 1024: n = 1024)

| kernel | dd | td | qd | MPFR 256 | MPFR 512 | MPFR 1024 |
|---|---:|---:|---:|---:|---:|---:|
| gemm NN | 100 | 101 | 99 | 46 | 34 | 14 |
| gemm update (k = 64) | 96 | 100 | 99 | 41 | 32 | 15 |
| syrk Upper/Trans | 93 | 98 | 96 | 32 | 21 | 12 |
| syrk Lower/NoTrans | 51 | 50 | 53 | 25 | 13 | 6 |
| cgemm NN | 61 | 60 | 57 | 19 | 9 | 9 |
| cherk Lower (k = 64) | 1.8 | 3.8 | 5.3 | 0.5 | 0.3 | 0.2 |
| trsm Left Lower NoTrans Unit | 1.4 | 4.5 | 7.3 | 2.3 | 2.0 | 1.5 |
| trsm Right Lower Trans (m × 64) | 3.4 | 5.3 | 6.4 | 3.1 | 2.6 | 1.1 |
| ctrsm Left Lower NoTrans Unit | 3.1 | 4.8 | 4.9 | 0.8 | 0.4 | 0.3 |
| getrf2 (LU panel n × 64) | 0.7 | 1.8 | 3.3 | 1.6 | 1.8 | 0.9 |
| potrf2 (n = 256) | 0.1 | 0.3 | 0.3 | 0.2 | 0.2 | 0.1 |
| trtri | 0.1 | 0.2 | 0.2 | - | - | - |
| lauum | 1.5 | 1.5 | 1.1 | - | - | - |
| laswp (memory bound) | 5.5 | 7.2 | 8.5 | - | - | - |

All kernels and sizes: `roofline_summary.md`, `roofline.csv`; graphs
`roofline_summary.png` (the table at n = 2048) and `roofline_<type>.png`
(% of roofline against n).

## Reproduce

    nvcc -O3 -std=c++17 -fmad=false -arch=sm_80 -U_GNU_SOURCE -D_XOPEN_SOURCE=700 -Iinclude \
         -Impblas/optimized/common -Iexternal/i/QD/include -Iexternal/i/MPC_CUDA/include \
         benchmark/peak_bench.cu -o peak_bench && ./peak_bench > peaks.csv
    PYTHONPATH=<matplotlib> python3 benchmark/roofline.py peaks.csv ../kernel-bench/kernels.csv ../gemm-sweep OUTDIR
