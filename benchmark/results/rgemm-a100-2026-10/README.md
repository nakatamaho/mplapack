# Rgemm on the A100: GPU against CPU, and the MPFR 512-bit size sweep (October 2026)

Host: AMD Ryzen Threadripper 3970X (32 cores), NVIDIA A100 80GB PCIe (sm_80),
CUDA 13.1, Ubuntu 26.04, branch claude/mpfr-cuda-rgemm (autotools build with
--enable-{dd,td,qd,mpfr}-cuda).  The host was shared with other jobs.  These
times go through the public Rgemm of the CUDA libraries, i.e. they **include**
the host-device transfers and, for MPFR, the packing into the device format.

## Files

| File | Contents |
|---|---|
| `../../Rgemm_gpu_compare.cpp` | Rgemm/Cgemm("N","N",n,n,n), best of R runs, GFLOPS; `cmp` mode compares two result dumps |
| `../../Rgemm_mpfr_sweep.cpp` | MPFR Rgemm("N","N",n,n,n) for n = powers of two and multiples of 13 up to NMAX |
| `../../plot_sweep.py` | draws a sweep CSV as in `mpfr512_sweep10k.png` |
| `build.sh` | builds both programs against the autotools libraries |
| `run_gemm_compare.sh` | runs the comparison for one type (CPU, default GPU kernel, other kernels, accuracy) |
| `gemm_compare_{dd,td,qd,mpfr}.txt` | its output (October 24-26, 2026) |
| `mpfr512_sweep.csv`, `.log`, `.png` | MPFR 512-bit sweep up to n = 2000 (best of 3) |
| `mpfr512_sweep10k.csv`, `.log`, `.png`, `_dark.png` | MPFR 512-bit sweep up to n = 9997 (best of 3 up to 2000, one run above; 783 points, about 18.3 hours) |

## Results in short

- MPFR 512 bits, GPU (cu_freal): 3.59 GFLOPS at n = 1024, 4.64 at 2048,
  5.26 at 4096, 6.03 at 8192, 6.11 at 9997.  CPU (32 cores): 0.285 GFLOPS
  (n = 512 and 1024): the GPU is 12.6 (n = 1024) to 21 (n = 9997) times faster.
- dd, n = 1024: CPU 14.9 GFLOPS, GPU 154 GFLOPS (see `gemm_compare_dd.txt`
  for td, qd, the Winograd and Ozaki variants and the accuracy columns).

## Reproduce

    benchmark/results/rgemm-a100-2026-10/build.sh ./bin
    ./bin/mpfr_sweep 512 10000 3 1 > sweep.csv 2> sweep.log      # PREC NMAX REPS REPS_LARGE [NSTART]
    benchmark/results/rgemm-a100-2026-10/run_gemm_compare.sh dd ./bin > gemm_compare_dd.txt

Plot (needs matplotlib and the Noto Sans CJK JP font):

    python3 benchmark/plot_sweep.py mpfr512_sweep10k.csv mpfr512_sweep10k.png \
      --subtitle 'Rgemm("N","N",n,n,n)、n ≤ 2000 は3回中の最良値、n > 2000 は1回。転送と詰め直しを含む。全783点、所要約18.3時間。' \
      --series "n = 13 の倍数" --cpu 0.285 --cpu-label "CPU 32コア(n=512, 1024 実測)"

The PNG files here were first drawn with an HTML/SVG page; `plot_sweep.py`
redraws the same graph from the CSV.  The environment variables in
`run_gemm_compare.sh` (`MPLAPACK_<P>_CUDA_KERNEL`, `_WINOGRAD_CUTOFF`,
`_OZAKI`, `_FORCE_RUNTIME`) are those of the CUDA libraries at the time.
