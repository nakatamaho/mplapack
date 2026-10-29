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
`Ctpmv`, `Rspr`, `Rspr2`, `Chpr`, `Chpr2`; Level 1 and 2 on request).
`Rgetrf` and `Rpotrf` also have device-resident real factorization entries
(described below). Every other routine, and every call that is not eligible, uses the same code as
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
result bit for bit. Except in `Rtrsm`, zeros are always `+0` on the
`cu_freal` engine.

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
- real `Rtrsm` uses `cu_fdiv_rn` on the fixed precisions (the additive
  `cu_fround.cuh` patch in `external/mpc_cuda/patches`); runtime precisions
  and `MPLAPACK_MPFR_CUDA_FORCE_RUNTIME=1` still use `cu_mpfr`. The existing
  gate on the host MPFR exponent range, device default `[1-2^30, 2^30-1]`,
  remains. Fixed calls also reject zero divisors and any intermediate
  exponent outside that range, leaving the caller's outputs untouched.
  Only Rtrsm preserves signed zeros in the otherwise unused sign field of
  zero fixed values, so its results remain bit-identical including zeros.
  Other fixed BLAS retain their documented +0 convention. Fixed Rtrsm now
  raises the stack only to a precision-dependent limit (4 KiB at 256/512 bits,
  8 KiB at 768/1024, 16 KiB at 2048), retaining the largest limit used to avoid
  expensive per-BLAS-call resizing. It no longer requests a 64-KiB/thread
  reservation; a larger pre-existing caller limit is not lowered. Resident
  factorization entries, in contrast, raise the limit
  only for their call and restore its previous value.
  Complex divisions (`cu_mpc_div`) are unchanged.

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
- `Rtrsm` timings are in `external/mpc_cuda/patches/README.md`; the other
  routines have no timings recorded here.

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
libmpc_cuda.a must have been built for the same ones. An externally installed
mpc_cuda must also provide `mpc_cuda/cu_fround.cuh`; apply the bundled patch
if its upstream version does not yet contain it. Both bundled build systems
apply the patch after extraction and re-extract/rebuild when it changes.

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

`mpfr_cuda_cu_fround_{host,device}` compare the new correctly rounded
operations directly with host MPFR at 224 (SB=32), 256, 512, 768, 1024 and
2048 bits: 100,000 random pairs per PB on the host, 1,000,000 on the device,
plus 256 deliberately constructed pairs, domain/alias tests and integer
midpoint tests. Output reports mismatches and exact-path counts separately
for division and sqrt. The device test skips with code 77 without a GPU.
See `external/mpc_cuda/patches/README.md` for the algorithm, reproducible
benchmarks (`mpfr_cuda_cu_fround_bench`, `mpfr_cuda_Rtrsm_round_bench`) and
verification results. These benchmarks are built with the CMake CUDA tests,
but are not registered as correctness tests.

## Device-resident real LU and Cholesky (goal 07)

`Rgetrf_mpfr_cuda` and `Rpotrf_mpfr_cuda` have the LAPACK arguments and return
`bool`; they never fall back to CPU factorization. `false` leaves A, pivots
and info untouched. The `Rgetrf` / `Rpotrf` dispatchers alone fall back to the
existing optimized CPU path. Complex factorization entries are **not yet
implemented**; Cgetrf/Cpotrf retain the previous GPU-BLAS-per-call path.

The real entries pack once, upload once, keep the matrix on the device for
the whole factorization, then download/unpack once. Only numerical status
words return to the host between panels. LU pivots are one-based, first
index on ties; singular LU completes, and Cholesky stops at the reference
failure index. Both retain iMlaenv block sizes, the reference recursive
panels and BLAS arithmetic order. No Winograd/sloppy arithmetic is enabled.
The backend `MPLAPACK_MPFR_LAPACK_FAST=1` currently retains this same algorithm.

Eligibility: default precision and every matrix element's precision agree,
MPFR_RNDN, finite inputs, and the default MPFR exponent range (the same gate
as Rtrsm). Fixed precisions 256/512/768/1024/2048 use cu_freal and correctly
rounded cu_fdiv_rn / cu_sqrt_rn; other precisions use cu_mpfr. Both engines
check intermediate range/domain failures before publication. Unlike the old
fixed GEMM convention, the factorization engine retains signed zeros in
all arithmetic, including updates. Padding and the opposite Cholesky triangle
are never published. Device-memory preflight includes the matrix, pivots,
candidate reduction, runtime arena/workspace and the device stack reservation.
Allocation/launch/transfer failure also returns false without publication.

Panels and updates use two nonblocking streams with ready/done events.
Runtime launches use disjoint physical thread-index ranges and two arena /
register-workspace pools: overlapping kernels must not reset or allocate from
another stream's per-thread slab. The arena globals remain constant until
both streams finish. The existing gpu_mutex serializes the whole call against
other MPFR BLAS.
Stack limits are raised only during the call and restored before publication;
fixed limits are those above, runtime 8 KiB through 1024 bits, 16 KiB through
4096 bits, and the existing conservative 64 KiB above that.

Settings are read on every factorization call:

| Variable | Meaning |
|---|---|
| `MPLAPACK_MPFR_LAPACK_CUDA=0` | disable resident factorization |
| `MPLAPACK_MPFR_LAPACK_CUDA_MIN=N` | minimum min(m,n); defaults Rgetrf 512, Rpotrf 2048 (measured below); one override applies to both |
| `MPLAPACK_MPFR_LAPACK_FAST=1` | currently the same reference-order algorithm |
| `MPLAPACK_MPFR_CUDA=0` | disable the MPFR GPU entries, including these |
| `MPLAPACK_MPFR_CUDA_RUNTIME=0` | disable the cu_mpfr engine |
| `MPLAPACK_MPFR_CUDA_FORCE_RUNTIME=1` | select cu_mpfr at fixed precisions too |
| `MPLAPACK_MPFR_CUDA_VERBOSE=1` | print measured upload / factor / download seconds |

CMake tests `mpfr_cuda_Rfactor_{host,device}[_runtime]` use a CPU-only oracle
and call the entries directly (no fallback can hide a failure). Device
absence returns 77. The autotools host variant is
`make -C mpblas/optimized/mpfr/cuda check`. Whole padded matrices, pivot
sentinels and info are compared strictly, including signed zeros; default
coverage is six precisions, square/tall/wide matrices, both triangles,
block boundaries and odd sizes through 599, singular/tied/tiny pivots,
nonpositive diagonals, backward error <=30, rejection and output preservation.
`OMP_NUM_THREADS=1 mpfr_cuda_Rfactor_device_test PREC --quick` selects the smaller same-kernel cases
for stack/memory sanitizer sweeps; default tests are not shortened.

Benchmark targets: `mpfr_cpu_Rfactor_bench` (CPU default),
`mpfr_cuda_Rfactor_bench` (end-to-end direct entry), and
`mpfr_cuda_Rfactor_panel_bench` (512-bit panel comparison). The latter times
the device fixed / device MPFR panel phase, or actual D2H/unpack/host-MPFR/
pack/H2D for a hybrid panel. The common trailing updates are excluded from
this isolated panel comparison, not estimated. `test/run_factor_bench.py`
records raw results from separate CPU / retained-old-CUDA / new-entry binaries,
interleaving order in multiple fresh-process rounds and reporting medians.
Upload/download include conversion; allocation, initialization and cleanup
are outside the split, so the parts need not sum to the end-to-end time.

### Panel choice and measured crossover (A100 sm_80, CUDA 13.1)

Choose cu_freal panels at the five fixed precisions and cu_mpfr otherwise
(or when forced). At 512 bits the fixed GPU panel is 3.66x/4.87x faster than
the runtime GPU LU panel, and 2.21x faster for Cholesky. The host hybrid is
faster for the small Cholesky diagonal, but would require repeated panel
conversion/transfers. It is a benchmark comparator, not an entry-point path:
the default retains the required GPU panels and one conversion at each end.
No CPU fallback or hybrid panel is hidden inside a successful entry.

Panel-phase seconds, median of five interleaved/reversed rounds; OMP=32,
MPFRXX_DEFAULT_PRECISION_BITS=512 in every process. LU has 64 columns;
Cholesky is one 64x64 diagonal. Fixed/runtime phase excludes initial
allocation/upload and final download; hybrid phase measures actual
D2H/unpack/CPU-panel/pack/H2D. No trailing-update time is inferred.

| panel | shape | GPU cu_freal | GPU cu_mpfr | host hybrid |
|---|---:|---:|---:|---:|
|Rgetrf2|512x64|0.035606|0.130347|0.284334|
|Rgetrf2|2048x64|0.037426|0.182221|0.286631|
|Rpotrf2 L|64x64|0.045240|0.100045|0.009627|

End-to-end seconds, **three interleaved rounds with reversed call order**,
medians, fresh process per call, OMP=32, deterministic full-precision random
significands and strict diagonal dominance. CPU is libmplapack_mpfr_opt;
old is the retained pre-goal CUDA shared library (CPU factors, GPU BLAS per
eligible call); new calls the entry directly, with MIN=0, and must return
true/info=0. All **126** timed calls succeeded. The shared host is not
calibrated; no times below are estimates. The split excludes initialization,
allocation and cleanup, and includes matrix conversion in upload/download.
Set MPFRXX_DEFAULT_PRECISION_BITS to the requested precision also for OpenMP
workers; the runner does this and clears inherited LAPACK/MPFR tuning.

|bits|routine|n|CPU|old CUDA|new entry|upload|factor|download|old/new|
|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|
|512|Rgetrf|512|1.756656|2.341467|0.848922|0.027765|0.287701|0.010484|2.758|
|512|Rpotrf L|512|0.566250|0.889964|0.984371|0.026733|0.409643|0.007654|0.904|
|512|Rgetrf|1024|10.162399|5.611150|1.355135|0.089125|0.683307|0.044113|4.141|
|512|Rpotrf L|1024|2.701463|1.411368|1.592530|0.087411|0.941978|0.029290|0.886|
|512|Rgetrf|2048|42.035304|13.834973|3.269151|0.303249|2.228172|0.163487|4.232|
|512|Rpotrf L|2048|15.898697|4.085477|3.513640|0.306834|2.524661|0.115267|1.163|
|512|Rgetrf|4096|207.326357|66.609945|14.326270|1.317213|11.433812|0.639304|4.649|
|512|Rpotrf L|4096|108.800011|20.468114|12.022977|1.306446|9.303057|0.450669|1.702|
|256|Rgetrf|512|2.357260|1.840826|0.721356|0.017675|0.171737|0.005474|2.552|
|256|Rpotrf L|512|0.453155|0.754467|0.754167|0.018071|0.193519|0.004226|1.000|
|1024|Rgetrf|512|2.360622|2.873891|1.404849|0.041723|0.811452|0.013564|2.046|
|1024|Rpotrf L|512|0.806218|1.297241|2.124351|0.041456|1.546730|0.010665|0.611|
|200|Rgetrf|512|1.910529|12.610459|1.720275|0.013817|1.141940|0.006778|7.330|
|200|Rpotrf L|512|0.457441|6.889667|1.677536|0.014755|1.064007|0.004351|4.107|

Use routine-specific defaults: **Rgetrf 512**, **Rpotrf 2048**. LU wins at
every measured size from 512; Cholesky loses at 512/1024 but wins at both
2048 and 4096. The losing points are retained above. These are measured
eligibility points, not an extrapolated exact crossover. One explicit MIN
overrides both. Other precisions retain this conservative shared policy.
Raw records (ignored local assets): bench_local/results/goal07_real/.

### Verification (real milestone, 2026-10-29)

All values below are observed, not inferred from CPU fallback. The four main
runs are host/device, each automatic and forced runtime. Each run compares
241 LU cases / 2,832,648 matrix elements and 84 Cholesky cases / 1,361,932
matrix elements **at each precision**. Entire padded matrices, pivots and
info are checked. Small nonsingular cases also require backward error <=30.
Signed zero is strict on both paths; no tolerance was weakened.

|bits|automatic engine|LU cases (four runs)|Cholesky cases (four runs)|A/ipiv/info errors|zero-sign differences|
|---:|---|---:|---:|---:|---:|
|256|cu_freal|241 x 4|84 x 4|0|0|
|512|cu_freal|241 x 4|84 x 4|0|0|
|768|cu_freal|241 x 4|84 x 4|0|0|
|1024|cu_freal|241 x 4|84 x 4|0|0|
|2048|cu_freal|241 x 4|84 x 4|0|0|
|200|cu_mpfr|241 x 4|84 x 4|0|0|

Each main host run has 710 checks / 4,194,760 element comparisons per
precision; each device run has 1,036 / 4,194,772 (including stack restoration
and memory-budget rejection). Four main runs total 20,952 checks and
100,674,384 element comparisons. All errors and zero-sign differences are 0.

- **CMake**: `ctest -R mpfr_cuda_ -j2`: **22/22 PASS, 12,492.95 s**,
  including all 18 pre-existing MPFR CUDA tests. Sanitizer jobs ran in
  parallel, but no benchmark did. CPU/CUDA old factor default/fast tests:
  **4/4 PASS, 1,174.61 s**. GPU-hidden device tests: **2/2 SKIP (77)**;
  GPU-hidden host quick test: **506 checks, 476,131 elements, 0 errors**.
- **Autotools**: regenerated Makefile, shared/static CUDA libraries and full
  host `make check`: **1/1 PASS, 2,588 s**, all six precisions, 710 checks
  each, no errors/zero differences. A previous full run also passed.
- **compute-sanitizer memcheck**, same factor kernels, quick cases including
  the 131-sized two-stream overlap: all six precisions **730 checks each**,
  LU 165 cases / 340,907 elements; Cholesky 58 / 135,044; total 476,143
  comparisons per precision, **0 numeric/zero-sign differences, 0 sanitizer
  errors**. Fixed limits 4/4/8/8/16 KiB, runtime 200-bit limit 8 KiB;
  every factor call restores its previous limit.
- **Optional 4096-bit runtime**: final normal quick test passes 730 checks /
  476,143 comparisons, LU 165 cases / Cholesky 58, no errors or zero-sign
  differences, with a 16-KiB scoped limit. Memcheck was stopped after over
  one hour without a completed result: **UNVERIFIED**, not counted as PASS.
- **Fixed Rtrsm**: post-policy full blas3 device replay **1/1 PASS,
  1,084.25 s**; memcheck **327,920 elements, 0 failures, 0 sanitizer errors**.
  The suite's 580 allowed sign-only zeros belong to other fixed BLAS;
  Rtrsm's zero-sign checks remain strict. A separate public-entry probe
  observes retained limits 4/4/8/8/16/16 KiB for 256/512/768/1024/2048/256
  bits, exact results including zeros, and the factor entry preserves the
  inherited 16-KiB limit. The BLAS bridge retains its right-sized maximum
  to avoid repeated CUDA stack reallocations; the factor guard is scoped.
- **LAPACK real GE/PO drivers**: CPU and LD_PRELOAD CUDA both pass
  3,653 + 5,748 GE and 1,628 + 1,910 PO checks (12,939 each), including
  error exits. CPU 25.80 s; CUDA with MIN=0 1,532.51 s, **1,815 LU and
  828 Cholesky successful resident-entry traces**. MIN=0 deliberately
  exercises even tiny, latency-dominated calls; it is not a speed claim.
- **Natural default gates**: 512-sized LU / 2048-sized Cholesky eligible,
  one below each rejected with outputs untouched: 6 checks / 8,912,896
  comparisons, 0 errors or sign differences.
- **ABI/build hygiene**: `nm -D` has only the two declared entry additions,
  no removals (3,983 -> 3,989 including four generated CUDA registration/
  wrapper symbols). Seven produced/dependent archives pass
  `check_unique_symbols.sh`; both final CUDA archives rechecked. All five
  generated accelerator headers and generated backend Makefiles are fresh;
  `git diff --check` passes. Both sm_80 and sm_90 fatbin targets build.

**Deferred / UNVERIFIED**: Cgetrf/Cpotrf are intentionally not implemented
in this user-approved real-only milestone; goal 07 and hence the full goal
08 are **not** declared complete/unblocked. The optional 4096-bit memcheck,
sm_90 execution, actual injected
cudaMalloc/copy/kernel failures, an n=10000 allocation, and crossovers beyond
the measured table are UNVERIFIED. Preflight budget/overflow rejection and
unchanged outputs are tested. No common dd/td/qd header, solver or bundled
mpc_cuda patch is changed.

## Limitations

- Accelerated: `Rgemm`, `Rsyrk`, `Rsyr2k`, `Rsymm`, `Rtrmm`, `Rtrsm`, the
  complex Level 3 routines and `Rgemmtr` (Level 1 and 2 on request).
- One thread per element of `C`, no shared-memory tiling yet.
- The Winograd variant exists for the fixed precisions only (not for the
  `cu_mpfr` kernels) and allocates its temporaries with `cudaMalloc` at every level.
- `Rtrsm` has been measured on A100 (see the patch README); no timings for
  the remaining routines are recorded here.
