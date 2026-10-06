# Optimized Rgemm for dd, qd and binary128

Header-only code shared by `libmplapack_{dd,qd,binary128}_opt` (CPU) and
`libmplapack_{dd,qd}_opt_cuda` (GPU). It covers the matrix-multiplication
features of [BNCmatmul](https://github.com/tkouya/bncmatmul): SIMD, cache
blocking, OpenMP, Strassen/Winograd and the Ozaki scheme on the CPU, and naive,
tiled, Winograd and Ozaki kernels on the GPU. The code is an independent
implementation of the algorithms, not a copy of BNCmatmul (LGPL).

| | dd (106 bits) | qd (212 bits) | binary128 (113 bits) |
|---|---|---|---|
| CPU conventional | blocked + OpenMP + AVX2 SIMD, **bit-identical** | blocked + OpenMP, **bit-identical** | blocked + OpenMP, **bit-identical** |
| CPU branch-free SIMD | (the default is branch-free) | `MPLAPACK_QD_GEMM_BF=1` (AVX2) | no hardware SIMD for binary128 |
| CPU Winograd | `MPLAPACK_DD_GEMM_WINOGRAD_CUTOFF=N` | `MPLAPACK_QD_GEMM_WINOGRAD_CUTOFF=N` | `MPLAPACK_BINARY128_GEMM_WINOGRAD_CUTOFF=N` |
| CPU Ozaki | `MPLAPACK_DD_GEMM_OZAKI=1` | `MPLAPACK_QD_GEMM_OZAKI=1` | `MPLAPACK_BINARY128_GEMM_OZAKI=1` |
| GPU naive / tiled | bit-identical to the CPU | bit-identical to the CPU | (CPU only) |
| GPU Winograd | `MPLAPACK_DD_CUDA_WINOGRAD_CUTOFF=N` | `MPLAPACK_QD_CUDA_WINOGRAD_CUTOFF=N` | |
| GPU Ozaki (cuBLAS DGEMM) | `MPLAPACK_DD_CUDA_OZAKI=1` | `MPLAPACK_QD_CUDA_OZAKI=1` | |

"Bit-identical" means that each element of C goes through the same operations,
in the same order, as in `openmp/Rgemm_{NN,NT,TN,TT}_omp.cpp`. Only the loops
over i and j are blocked and reordered. Every combination of transposes is
supported.

Winograd, Ozaki and the qd branch-free kernel are opt-in. Their results differ
from the conventional ones in rounding, and their errors are bounded normwise,
not elementwise.

## Files

| File | Contents |
|---|---|
| `Rgemm_blocked_common.h` | blocked engine (packing, OpenMP over C blocks), generic scalar kernel |
| `Rgemm_simd_kernels_common.h` | dd kernel, qd branch-free kernel (vectorized over i) |
| `Rgemm_simd_common.h` | AVX2 / baseline function multiversioning (x86-64 ELF), `omp simd` |
| `Rgemm_qd_ops_common.h` | dd/qd arithmetic on doubles for host and device, identical to libqd |
| `Rgemm_winograd_common.h` | Winograd recursion (BNCmatmul's schedule), CPU backend |
| `Rgemm_ozaki_common.h` | Ozaki splitting, exact integer DGEMM, accumulation |
| `Rgemm_cpu_common.h` | CPU driver: environment variables, choice of algorithm |
| `Rgemm_cuda_ops_common.h`, `Rgemm_cuda_common.cuh` | GPU kernels, and their host emulation |
| `Rgemm_cuda_bridge_common.h` | GPU eligibility, choice of algorithm, CPU fallback |

Per type, `<type>/openmp/Rgemm_blocked_omp.cpp` instantiates the CPU code and
is called by `<type>/Rgemm.cpp`. `<type>/cuda/{Rgemm.cpp,
Rgemm_bridge_cuda.cpp, Rgemm_gpu_cuda.cu}` instantiate the GPU code.

## CPU

Environment variables are read at every call. `<P>` is one of `MPLAPACK_DD`,
`MPLAPACK_QD` or `MPLAPACK_BINARY128`.

| Variable | Meaning |
|---|---|
| `<P>_GEMM_BLOCKED=0` | use the plain OpenMP loops |
| `<P>_GEMM_MIN_MNK=N` | smallest `m*n*k` for the blocked engine (default 32768) |
| `<P>_GEMM_WINOGRAD_CUTOFF=N` | Winograd when `min(m,n,k) > N` (default 0: off) |
| `<P>_GEMM_OZAKI=1` | Ozaki scheme (takes precedence over Winograd) |
| `<P>_GEMM_OZAKI_SLICES=S` | number of slices (default: enough for precision + 8 bits) |
| `MPLAPACK_QD_GEMM_BF=1` | qd: branch-free SIMD kernel (Kouya's QWBFAdd/QDBFMul) |

SIMD:

- **dd:** libqd's `dd_real` IEEE addition and multiplication have no branches.
  The kernel performs them on the high and low words of 4 elements of C at a
  time (AVX2 + FMA). It uses FMA exactly where libqd does (`QD_FMS`).
- **qd:** libqd's IEEE `qd_real` addition branches on the data, so the
  conventional qd kernel is scalar. `MPLAPACK_QD_GEMM_BF=1` uses
  `qd_real::bf_add` / `bf_mul` instead, which can be vectorized; its results
  are bit-identical to libqd's `bf_*` functions, not to the default
  `qd_real` operators.
- **Instruction sets:** on x86-64 ELF systems the kernels are built twice,
  for AVX2 + FMA (x86-64-v3) and for the baseline, and the variant is chosen
  when the library is loaded. AVX-512 is not used. Both variants produce the
  same results.

The Ozaki scheme cuts the rows of op(A) and the columns of op(B), after
scaling by powers of two, into slices of
`beta = floor((53 - ceil(log2 k)) / 2)` bits. The slice products are integer
matrix products below 2^53, so they are exact in double precision. They are
summed in the target precision, slice pairs with `s + t < S` only. Calls fall
back to the other paths when an entry is not finite or lies outside
[2^-960, 2^960].

Measured on 4 cores (Xeon, AVX2 variant), seconds for an n x n x n product.
The baseline is the former `Rgemm`, i.e. `<P>_GEMM_BLOCKED=0`.

| | former | blocked | Winograd | Ozaki | qd BF |
|---|---|---|---|---|---|
| dd, n = 1000, NN | 1.42 | 0.30 | 0.52 (cutoff 128) | 0.77 | |
| dd, n = 1000, TN | 2.39 | 0.33 | | | |
| qd, n = 500, NN | 3.27 | 3.18 | 2.67 (cutoff 64) | 0.33 | 0.33 |
| binary128, n = 500, NN | 1.22 | 1.27 | 1.18 (cutoff 64) | 0.15 | |

## GPU (dd, qd)

`libmplapack_dd_opt_cuda` (`--enable-cuda`, `-DMPLAPACK_ENABLE_CUDA=ON`) and
`libmplapack_qd_opt_cuda` (`--enable-qd-cuda`, `-DMPLAPACK_ENABLE_QD_CUDA=ON`)
send eligible `Rgemm` calls to the GPU. Everything else, and every call whose
GPU run fails, uses the CPU code above.

| Variable (`<P>` = `MPLAPACK_DD` or `MPLAPACK_QD`) | Meaning |
|---|---|
| `<P>_CUDA=0` | never use the GPU |
| `<P>_CUDA_MIN_MNK=N` | smallest `m*n*k` sent to the GPU (default 32768) |
| `<P>_CUDA_KERNEL=naive` | one thread per element, no shared memory (default: `tiled`, 16 x 16 shared-memory tiles) |
| `<P>_CUDA_WINOGRAD_CUTOFF=N` | Winograd when `min(m,n,k) > N`, temporaries on the GPU, tiled leaves |
| `<P>_CUDA_OZAKI=1` | Ozaki: slices made on the host, products by cuBLAS DGEMM, sums on the GPU |
| `<P>_CUDA_OZAKI_SLICES=S` | number of slices |
| `<P>_CUDA_VERBOSE=1` | report CUDA errors that cause a CPU fallback |
| `MPLAPACK_DD_CUDA_LEGACY=1` | dd: the former 2010-2011 kernels (`Rgemm_fermi.cu`, sloppy dd arithmetic) |

The naive and tiled kernels are bit-identical to the CPU `Rgemm`. The device
functions reproduce libqd's operations, compiled with `-fmad=false`. The GPU
Ozaki result is bit-identical to the CPU Ozaki result.

The GPU is used only when libqd is configured with IEEE add and accurate
multiplication. This is how `external/qd` is built.

Limitations:

- one element of C per thread;
- qd's IEEE addition uses a 192-byte local-memory frame per thread;
- Ozaki keeps all slices on the GPU (`S (mk + kn)` doubles);
- not yet measured on a GPU.

## Benchmarks

`benchmark/go.Rgemm_algo.sh` runs `Rgemm.{dd,qd,binary128}_opt` once per
algorithm: blocked, former, Winograd, Ozaki, and qd BF. It plots MFLOPS against
the dimension, one PDF page per type (`Rgemm_algo.plt.in`).

`benchmark/go.Rgemm_cuda.sh` does the same for `Rgemm.{dd,qd}_cuda_total`:
tiled, naive, Winograd, Ozaki, and the former dd kernels. It also runs
`Rgemm.mpfr_cuda_total` at 512 and 1024 bits, using the `-PREC` option of the
MPFR `Rgemm` benchmark.

Both build systems build the programs and generate the scripts:

- autotools: `Makefile.{qd,mpfr}_cuda.am`;
- CMake: `-DMPLAPACK_BUILD_BENCHMARKS=ON`, which also builds the `*_opt` and
  `*_cuda_total` programs.

## Tests

| CMake test | autotools `make check` in | Checks |
|---|---|---|
| `{dd,qd,binary128}_opt_Rgemm` | `mpblas/optimized/<type>` (`Rgemm_opt_test`) | blocked = OpenMP loops bit for bit; Winograd (cutoffs 1, 2, 3, 8), Ozaki, qd BF exact on integers and within `2^(20-p) (k|alpha| max|A| max|B| + |beta| max|C|)` on reals; Ozaki fallback |
| `{dd,qd}_cuda_Rgemm_host` | `mpblas/optimized/{dd,qd}` (`Rgemm_gpu_host_test`) | GPU code run on the host: naive = CPU bit for bit, Ozaki = CPU Ozaki bit for bit, Winograd exact / bounded |
| `{dd,qd}_cuda_Rgemm_device` | | the same on the GPU, including the tiled kernel; skipped without a device |

autotools builds the tests from `<type>/test/*.cpp`, which include the files
of `common/test`.

The host emulation runs the naive element code. The tiled kernel, which uses
shared memory and `__syncthreads`, is covered only by the device test.
