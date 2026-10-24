# libmplapack_mpfr_opt_cuda

MPFR flavor of MPLAPACK whose `Rgemm` and Level 3 BLAS (`Rsyrk`, `Rsyr2k`,
`Rsymm`, `Rtrmm`, `Rtrsm`) run on an NVIDIA GPU, and, for debugging, Level 1
and 2 BLAS (`Rgemv`, `Rsymv`, `Rtrmv`, `Rger`, `Rsyr`, `Rsyr2`, `Rscal`,
`Rswap`, `Rrot`) with `MPLAPACK_MPFR_CUDA_L12=1`; likewise the complex BLAS
(`Cgemm`, `Csymm`, `Chemm`, `Ctrmm`, `Ctrsm`, `Csyrk`, `Cherk`, `Csyr2k`,
`Cher2k`, `Cgemmtr`, `Rgemmtr`; Level 1 and 2: `Cgemv`, `Cgeru`, `Cgerc`,
`Cher`, `Cher2`, `Caxpy`, `Ccopy`, `Cswap`, `Cscal`, `CRscal`, `CRrot`) and the
band/packed/Hermitian/triangular Level 2 routines (`Chemv`, `Chbmv`, `Chpmv`,
`Rsbmv`, `Rspmv`, `Rgbmv`, `Cgbmv`, `Ctrmv`, `Rtbmv`, `Ctbmv`, `Rtpmv`,
`Ctpmv`, `Rspr`, `Rspr2`, `Chpr`, `Chpr2`; Level 1 and 2 on request). Every
other routine, and every call that is not eligible, uses the same code as
`libmplapack_mpfr_opt`.

The GPU arithmetic comes from [mpc_cuda](https://github.com/tkouya/mpc_cuda)
(bundled as `external/mpc_cuda`):

| Precision | Kernels | Arithmetic | Result vs CPU `Rgemm` |
|---|---|---|---|
| 256, 512, 768, 1024, 2048 bits | `Rgemm_device_cuda.cu` | `cu_fp::cu_freal<PB>`: fixed precision, significand in registers | identical, except that a zero is always `+0` |
| any other | `Rgemm_rt_cuda.cu` (port of mpc_cuda `demos/matmul_mpfr.cu`) | `cu_mpfr`: MPFR 4.2.2 compiled for the device | identical, signed zeros included |

Both perform the same operations in the same order as
`openmp/Rgemm_*_omp.cpp`. The fixed precisions are listed once, in
`MPLAPACK_MPFR_CUDA_FIXED_PRECISIONS` (`Rgemm_kernel_cuda.h`). mpc_cuda
reports that `cu_freal` keeps the significand in registers up to about 1024
bits; at 2048 bits it spills to local memory and is several times slower per
operation, which may make the `cu_mpfr` kernels
(`MPLAPACK_MPFR_CUDA_FORCE_RUNTIME=1`) competitive there.

### Winograd (Strassen) variant

With `MPLAPACK_MPFR_CUDA_WINOGRAD_CUTOFF=N` (N > 0), fixed-precision calls with
`min(m, n, k) > N` use the Winograd variant of Strassen's algorithm
(`Rgemm_winograd_cuda.h`), 7 block products instead of 8 per level, recursing
until a dimension is at most N; an odd dimension is peeled at that level (the
even part recurses on views of the operands, pointers and leading dimensions,
without padded copies; the last row, column and step of k are added by the
conventional product on views). The schedule follows `mul_mpfmatrix_winograd_even()` of
[BNCmatmul](https://github.com/tkouya/bncmatmul). It is off by default: the
result differs from the conventional `Rgemm` in rounding, and its error is
bounded normwise rather than elementwise. All temporaries stay on the GPU.

## Level 1/2/3 BLAS

The kernels are in `blas3_kernel_cuda.h`, the host side in
`blas3_bridge_cuda.cpp`. They perform, for every element of the result, the
operations of the reference BLAS (`mpblas/reference/R*.cpp`, which the MPFR
libraries use for these routines) in the same order. Each operation is
correctly rounded to the working precision, so the result is the reference
result bit for bit. Zeros are always `+0` on the `cu_freal` engine.

The work is split into independent lanes:

| Routine | Lane |
|---|---|
| `Rsyrk`, `Rsyr2k`, `Rsyr`, `Rsyr2` | one element of the referenced triangle |
| `Rsymm`, `Rtrmm`, `Rtrsm` with side L | one column |
| `Rsymm`, `Rtrmm`, `Rtrsm` with side R | one row |
| `Rger` | one element of A |
| `Rgemv`, `Rsymv`, `Rtrmv`, Level 1 | one element of y (`Rtrmv`: x) |

The lanes are those of `../../common/Rblas_lanes_common.h`, shared with the
CPU (OpenMP) and dd/td/qd GPU versions; the host side
(`blas3_bridge_cuda.cpp`) takes the call descriptions of
`../../common/Rblas_calls_common.h`.

Each routine is written once against an arithmetic engine:

- 256, 512, 768, 1024 and 2048 bits use `cu_freal`;
- other precisions (or `MPLAPACK_MPFR_CUDA_FORCE_RUNTIME=1`) use `cu_mpfr`;
- `Rtrsm` always uses `cu_mpfr`, because `cu_freal` has no correctly rounded
  division. It runs only when the host MPFR exponent range is the device
  default `[1-2^30, 2^30-1]`, so that overflow and underflow happen as on the
  host.

`cuda/R*.cpp` replaces the CPU dispatcher of `../../common/blas`. Calls that
do not run on the GPU go to the CPU code of `../../common/blas/R*_omp.cpp`:
the lane-parallel OpenMP version, or the reference code, compiled unchanged
under another name. This includes invalid arguments, quick returns,
`alpha == 0` (Level 2/3), zero increments and the cases listed for `Rgemm`
below. The size threshold is `MPLAPACK_MPFR_CUDA_MIN_MNK`, compared with the
number of multiply-adds (`n*n*k` for `Rsyrk`, `m*n*n` or `m*m*n` for
`Rsymm`/`Rtrmm`/`Rtrsm`, `m*n` for `Rgemv`, ...). `MPLAPACK_MPFR_CUDA_L3=0`
disables Level 3 on the GPU; Level 1 and 2 run there only with
`MPLAPACK_MPFR_CUDA_L12=1`, since their transfers dominate.

Limitations:

- with side L or side R the parallelism is the number of columns or rows only;
- a call copies its operands to the GPU and back;
- the routines are not yet measured on a GPU.

## Complex BLAS

The complex routines (and `Rgemmtr`) run the text of the reference routines
on the GPU: `../../common/blas_dev/` holds it as templates
(`misc/gen_blas_device_reference.py`), `Cblas_mpfr_dev_types.h` the device
types on `cu_freal<PB>`, and every GPU thread runs the reference routine on
one column, row or element of the output (the jobs of
`../../common/Cblas_cuda_jobs_common.h`, shared with dd/td/qd). gmpfrxx_mkII
evaluates each operator of the reference text as one MPFR/MPC call rounded to
nearest, never fused, and the device types perform the same call:

- 256, 512, 768, 1024 and 2048 bits: `cu_freal<PB>` (`mpc_mul` is correctly
  rounded per component, as `cu_cmul`; divisions are `cu_mpc_div`); zeros
  the GPU computes are `+0`;
- other precisions up to 4096 bits, 16384 bits when built with
  `MPLAPACK_MPFR_CUDA_RT_LARGE` (CMake option of the same name; its device code
  takes long to compile), or `MPLAPACK_MPFR_CUDA_FORCE_RUNTIME=1`:
  `rreal`/`rcomplex`, the MPFR/MPC functions of libmpc_cuda themselves
  (`cu_mpfr_add`, `cu_mpc_mul`, `cu_mpc_div`, ...), signed zeros included.

So the results are the reference ones bit for bit, and elements the routine
leaves unchanged keep their sign. The host side is `Cblas_mpfr_bridge_cuda.cpp` on
`../../common/blas_cuda/Cblas_bridge_body.h`, with the dispatchers of
`../../common/blas_cuda/C*.cpp`.

Calls stay on the CPU (`../../common/blas/C*_omp.cpp`, then the reference)
when:

- the default precision (`mpfr_class`, both components of `mpc_class`) is
  above that limit, or not one of the fixed precisions with
  `MPLAPACK_MPFR_CUDA_RUNTIME=0`;
- an element of an operand or a scalar has another precision, is not finite,
  or has an exponent outside the bound that keeps every intermediate inside
  the exponent range; a result is outside the range (`Ctrsm`: outside the
  bound, since its results are factors of later operations);
- `Ctrsm` has diag = N and a zero diagonal element (the divisions of
  `Ctrsm` are `mpc_div`, on the device libmpc_cuda's `cu_mpc_div`);
- the same settings as for the real routines (`MPLAPACK_MPFR_CUDA`,
  `MPLAPACK_MPFR_CUDA_L3`, `MPLAPACK_MPFR_CUDA_L12`,
  `MPLAPACK_MPFR_CUDA_MIN_MNK`, rounding mode) exclude them.

## Building

mpc_cuda needs CUDA 13 (nvcc 12.0 rejects its early-clobber `"=&l"` asm
constraints). The bundled one is built for the GPU architectures given
(several give fat binaries, e.g. 80 and 90 for A100 and H100); an installed
libmpc_cuda.a must have been built for the same ones.

CMake (builds `external/mpc_cuda` unless `MPLAPACK_MPC_CUDA_INCLUDE_DIR` names
an installed mpc_cuda):

    cmake -B build -DMPLAPACK_ENABLE_MPFR=ON -DMPLAPACK_ENABLE_OPT=ON \
          -DMPLAPACK_ENABLE_MPFR_CUDA=ON \
          -DMPLAPACK_CUDA_ARCHITECTURES=121        # e.g. 121 for GB10, "80;90" for A100 and H100
    # installed mpc_cuda instead of the bundled one:
    #     -DMPLAPACK_MPC_CUDA_INCLUDE_DIR=/opt/mpc_cuda/include
    #     [-DMPLAPACK_MPC_CUDA_LIBRARY=/opt/mpc_cuda/lib/libmpc_cuda.a]

autotools (builds `external/mpc_cuda` unless `--with-mpc-cuda-includedir`
names an installed mpc_cuda):

    ./configure --enable-mpfr --enable-mpfr-cuda --with-mpfr-cuda-arch=121 \
                --with-cudatoolkithome=/usr/local/cuda
    # A100 and H100 (fat binaries): --with-mpfr-cuda-arch=80,90
    # installed mpc_cuda instead of the bundled one:
    #     --with-mpc-cuda-includedir=/opt/mpc_cuda/include
    #     [--with-mpc-cuda-library=/opt/mpc_cuda/lib/libmpc_cuda.a]

The CUDA sources are compiled as relocatable device code and device-linked
with `libmpc_cuda.a` inside the library, so programs link it with the C++
compiler: `-lmplapack_mpfr_opt_cuda -lcudart` (pkg-config:
`mplapack_mpfr_opt_cuda`).

## When the GPU is used

`Rgemm` sends a call to the GPU when all of the following hold; otherwise it
runs on the CPU and the result is unchanged:

- the MPFR default precision, `alpha`, `beta` and every referenced element of
  `A`, `B` and `C` have the same precision;
- the MPFR default rounding mode is `MPFR_RNDN`;
- every referenced value is finite (`C` is not read when `beta == 0`);
- exponents are small enough that no intermediate can leave the MPFR exponent
  range (host and device), and every result fits in it;
- `m*n*k >= MPLAPACK_MPFR_CUDA_MIN_MNK` (default 32768);
- a CUDA device is available and `MPLAPACK_MPFR_CUDA` is not `0`.

Environment variables:

| Variable | Meaning |
|---|---|
| `MPLAPACK_MPFR_CUDA=0` | never use the GPU |
| `MPLAPACK_MPFR_CUDA_MIN_MNK=N` | minimum `m*n*k` for the GPU path |
| `MPLAPACK_MPFR_CUDA_RUNTIME=0` | no GPU for precisions other than 256/512/768/1024/2048 |
| `MPLAPACK_MPFR_CUDA_FORCE_RUNTIME=1` | use the `cu_mpfr` kernels for 256/512/768/1024/2048 bits too |
| `MPLAPACK_MPFR_CUDA_WINOGRAD_CUTOFF=N` | Winograd/Strassen for the fixed precisions when `min(m,n,k) > N` (default 0: off) |
| `MPLAPACK_MPFR_CUDA_RT_BLOCKS`, `MPLAPACK_MPFR_CUDA_RT_THREADS` | thread pool of the `cu_mpfr` kernels (default 256 x 32) |
| `MPLAPACK_MPFR_CUDA_VERBOSE=1` | print CUDA errors that cause a CPU fallback |

The `cu_mpfr` kernels take the limbs MPFR allocates internally from a
per-thread arena of 16 KB per 1024 bits of precision, reset at every step of
the dot product (as in mpc_cuda's matmul demo); they raise the device stack
limit to 64 KB while they run.

## Tests (CMake)

`mpfr_cuda_Rgemm_{host,device}[_runtime]` compare the GPU `Rgemm` bit for bit
with the CPU `Rgemm` at 256, 512, 768, 1024, 2048 (fixed) and 64, 200, 333,
4096 (runtime) bits, for all
transpose combinations and several `alpha`/`beta`, with zeros and mixed
magnitudes, and check that ineligible calls fall back. `_runtime` runs every
precision through the `cu_mpfr` kernels.

- `host`: the kernels' element code runs on the host, with `cu_mpfr` mapped
  to the system MPFR (libmpc_cuda's own host code needs a GPU). No GPU needed.
- `device`: the kernels run on the GPU; skipped when no device is present.

`mpfr_cuda_blas3_{host,device}[_runtime]` compare the Level 3 routines bit
for bit with the reference. They cover every side, uplo, trans and diag,
`alpha` in {1, -1/2, random}, `beta` in {0, 1, 3/4, -5/4}, and odd sizes and
leading dimensions, at 256, 512, 768, 1024, 2048, 64, 200 and 333 bits.
Elements outside the referenced part must be left untouched.

`mpfr_cuda_Rblas_{host,device}` run `../../common/test/Rblas_omp_test.cpp`
through the dispatchers with every call, Level 1 and 2 included, on the GPU
(`cu_mpfr` engine), and compare with the reference routines bit for bit.

`mpfr_cuda_Cblas_{host,device}` run `../../common/test/Cblas_omp_test.cpp`
the same way for the complex routines with the MPFR/MPC engine, and
`mpfr_cuda_Cblas_{host,device}_fixed` with the `cu_freal` engine
(`MPLAPACK_MPFR_CUDA_FORCE_RUNTIME=0`), counting the zeros that differ only
in sign; the default precision is 512 bits (`MPFRXX_DEFAULT_PRECISION_BITS`
selects another).

`mpfr_cuda_Rgemm_winograd_{host,device}` check the Winograd path at
256/512/768/1024/2048 bits with cutoffs 1, 2, 3 and 8: integer matrices (no rounding) must
equal the CPU `Rgemm` exactly for odd/even sizes, all transpose combinations
and general `alpha`/`beta`; random real matrices must agree within
`2^(20-p) * (k |alpha| max|A| max|B| + |beta| max|C|)`; and the GPU backend
must equal the host backend bit for bit.

They run with `OMP_NUM_THREADS=1`: the OpenMP CPU `Rgemm` used as the
reference creates its temporaries in each worker thread at that thread's MPFR
default precision, which is not the caller's precision.

## Limitations

- Accelerated: `Rgemm`, `Rsyrk`, `Rsyr2k`, `Rsymm`, `Rtrmm`, `Rtrsm`, the
  complex Level 3 routines and `Rgemmtr` (Level 1 and 2 on request).
- One thread per element of `C`, no shared-memory tiling yet.
- The Winograd variant exists for the fixed precisions only (not for the
  `cu_mpfr` kernels) and allocates its temporaries with `cudaMalloc` at every level.
- Not yet measured on a GPU.
