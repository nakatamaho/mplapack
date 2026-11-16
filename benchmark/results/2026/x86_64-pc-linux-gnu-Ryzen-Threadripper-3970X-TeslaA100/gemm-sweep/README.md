# Rgemm GPU kernel-only sweep, NN / NT / TN / TT (November 4, 2026)

Rgemm(TRANSA, TRANSB, n, n, n) on the GPU with the data already on the device,
timed with CUDA events: **no host-device transfers, no MPFR packing**
(`benchmark/Rkernel_bench.<type>_cuda_kernel gemm <FORM> ...`).

- n = multiples of 7, multiples of 31 and powers of two (712 points up to 4096).
- n ≤ 2000: one warm-up, then the best of 3 runs; n > 2000: one run, no warm-up
  (as the 2026-10 MPFR sweep in ../rgemm/).
- Upper limits: dd/td/qd all forms and MPFR 256/512 NN to 4096; MPFR 256/512
  NT/TN/TT to 2048; MPFR 1024 NN to 2048, NT/TN/TT to 1024 (time).
- Host: Threadripper 3970X, A100 80GB PCIe, CUDA 13.1. No other process was on
  the GPU at the start or end (`gpu-processes-at-*.txt` are empty). The host
  load average during the run is in `run.log` (2-15; it does not affect a
  GPU-kernel timing much). Started 2026-11-04 15:24, NN finished 19:34, all
  forms 22:46.

## Results (GFLOPS = 2n³ / time, at the largest power of two measured)

| type | NN | NT | TN | TT |
|---|---:|---:|---:|---:|
| dd | 322 (n=4096) | 322 (n=4096) | 323 (n=4096) | 321 (n=4096) |
| td | 99.8 (n=4096) | 99.7 (n=4096) | 99.7 (n=4096) | 99.7 (n=4096) |
| qd | 45.9 (n=4096) | 45.8 (n=4096) | 45.9 (n=4096) | 45.9 (n=4096) |
| MPFR 256 (cu_freal) | 15.8 (n=4096) | 15.7 (n=2048) | 15.0 (n=2048) | 15.1 (n=2048) |
| MPFR 512 (cu_freal) | 6.31 (n=4096) | 6.15 (n=2048) | 5.80 (n=2048) | 5.80 (n=2048) |
| MPFR 1024 (cu_freal) | 1.02 (n=2048) | 0.97 (n=1024) | 0.94 (n=1024) | 0.94 (n=1024) |

The four forms give the same speed for dd/td/qd; for MPFR, TN/TT (A transposed)
are about 5% slower than NN/NT. Every type reaches its plateau by n ≈ 1000
(MPFR 256 has a dip around n = 230-250).

## Files

- `<type>_<FORM>.csv`: n,seconds,gflops; `run.log`: time and load per type and form.
- `gemm_<type>.png`: the four forms of one type; `gemm_all_NN.png`: NN of every type (log scale).

## Reproduce

    benchmark/run_gemm_sweep.sh OUTDIR TYPE NMAX [FORM]      # TYPE = dd|td|qd|mpfr<bits>, FORM = NN|NT|TN|TT
    PYTHONPATH=<matplotlib> python3 benchmark/plot_gemm_sweep.py gemm_dd.png --title ... NN=dd_NN.csv NT=dd_NT.csv TN=dd_TN.csv TT=dd_TT.csv
