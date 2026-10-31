# mpc_cuda patches / proposed upstream PR

Status: **submitted upstream** as https://github.com/tkouya/mpc_cuda/pull/1
(2026-10-30; the upstream branch also adds tools/test_fround.cu to
`make check-fixed`, and documents the API in doc/manual.md and manual.ja.md).
Base: upstream https://github.com/tkouya/mpc_cuda commit **abd9407**, bundled
as `mpc_cuda-0.0.1-abd9407.tar.xz`. The distfile and its recorded SHA256 are
unchanged. Remove the patch when the bundled upstream revision contains it.
This branch and 3.1 use the same patch/commit; no separate backport branch.

## Proposed PR: correctly rounded finite cu_freal division and square root

### Change

`patch-cu_freal-correctly-rounded-div-sqrt` adds the host/device header
`include/mpc_cuda/cu_fround.cuh`, includes it from `mpc_cuda.cuh`, and adds
it to the explicit `Makefile.in` header install list (the old rule only
wildcards `.h`, not `.cuh`). It adds, in namespace `cu_fp`:

```cpp
template<int PB> __host__ __device__ int
cu_fdiv_rn(cu_freal<PB>& q, const cu_freal<PB>& a, const cu_freal<PB>& b);
template<int PB> __host__ __device__ int
cu_sqrt_rn(cu_freal<PB>& r, const cu_freal<PB>& a);
```

Return 0 on success, 1 for a zero divisor or negative square-root input;
output on error is unspecified. Output may alias either input. These are
round-to-nearest, ties-to-even operations for finite normalized values.
Like cu_freal, they have no NaN, infinity or signed zero: zero is **+0**.
There is no MPFR emin/emax restriction, but exponents/results must fit in
`long` and avoid cu_freal's reserved zero exponent, as with its other
operations. The existing `cu_fdiv`, `cu_sqrt`, `operator/`, complex division,
and all other existing arithmetic are unchanged. No external ABI symbols.

### Algorithm and error bound

Normalize operands before Newton (including exponent parity for sqrt),
work at **W = PB + 64** independently of `CU_FGUARD`, and restore the
exponent after rounding. Reuse `cu_recip_w` / `cu_sqrt_w` and their rounded
multiplication. The header derives conservative relative-error recurrences
`e_next <= e^2 + 4u` and `e_next <= 2e^2 + 8u`, `u=2^-W`, from a double
seed error at most `2^-50`. The existing iteration count has two surplus
iterations beyond quadratic convergence. Including the last multiply the
error is below 22 W-ulps; use a conservative **64 W-ulp** ambiguity filter.

Explicit convergence argument (including, not discarding, the rounding
floor): both recurrences are bounded by `e_next <= 2e^2 + 8u`. Set
`B_0=2^-50`, `B_next=4B^2`, hence `B_j=2^(-2-48*2^j)`. Since
`36B_j + 162u < 1`, induction gives `e_j <= B_j + 9u` (expand
`2(B_j+9u)^2+8u <= 4B_j^2+9u`). With
`K=ceil(log2(W/52))`, at iteration `K+1`, `B_j/u < 2^-80` for W>=96.
The final, `(K+2)`th iteration therefore gives `<5u` for reciprocal and
`<9u` for rsqrt using their respective recurrences; the last multiply
adds at most another relative u (the header's 7u/11u bounds leave margin).

The low 64 significant bits of the approximation measure distance from a
PB midpoint, including **SB=32** layouts. Farther than 64 W-ulps, direct PB
rounding is safe. Otherwise construct the exact PB+1-bit midpoint `t` and
compare **a with t*b** or **a with t*t** using full `cu_umul` limb products,
integer exponent alignment and comparison. This is the sign of the exact
residual, without a rounded floating-point subtraction. Choose the correct
neighbor; exact midpoint ties use the PB low-bit parity. No unbounded Ziv
loop, allocator, global state or MPFR dependency. Instrumentable internal
helpers return per-call exact-path flags for tests, not global counters.

### Tests supplied by MPLAPACK

The upstream patch is deliberately independent of MPLAPACK. The oracle
program is `mpblas/optimized/mpfr/cuda/test/cu_fround_test.cpp`, also
included by its device `.cu` wrapper. It compares against system MPFR
`mpfr_div` / `mpfr_sqrt`, MPFR_RNDN, at PB=224,256,512,768,1024,2048.
It covers wide random exponents, signs, zero/domain errors, aliasing,
powers of two, all-one significands, exact quotients/squares, exponent
edges, and deliberately chosen near-midpoint values on both sides.

**Specification clarification:** with equal-P finite binary inputs and an
unbounded exponent range, an exact division midpoint would require an odd
(P+1)-bit numerator times an integer, which cannot fit a P-bit numerator.
Likewise a square-root midpoint squared has more than P significant bits.
Thus those public-operation ties cannot be constructed without first
rounding away the tie. Instead test the exact midpoint/tie resolver directly
with wider exact integers (both parity choices, midpoint and +/- one
integer unit for square comparisons), and separately force the actual
public-operation exact paths. For division choose `b=floor/ceil(2^P/3)` and
`a=b+1`: the quotient is on opposite sides of `1+3*2^-P` by O(2^-2P), within
one W-ulp. Additional random odd denominators use modular inverses modulo
`2^P` to construct midpoint residuals +/-1. `sqrt(1+2^(1-P))` is just below
`1+2^-P` by O(2^-2P); Hensel-lifting `k^2 == -7 (mod 2^(P+2))` constructs
PB-exact `a=(k^2+7)*2^(-2P)` whose root is just above the midpoint `k*2^-P`.
Unscaled hard cases assert that the public operation actually took the
exact path, for division and for sqrt on both sides of the boundary.

Full CMake host/device oracle runs (A100 80GB PCIe, CUDA 13.1):

| PB | host pairs / mismatches | device pairs / mismatches | exact div / sqrt, each run |
|---:|---:|---:|---:|
|224 (SB=32)|100256 / 0|1000256 / 0|80 / 48|
|256|100256 / 0|1000256 / 0|80 / 48|
|512|100256 / 0|1000256 / 0|80 / 48|
|768|100256 / 0|1000256 / 0|80 / 48|
|1024|100256 / 0|1000256 / 0|80 / 48|
|2048|100256 / 0|1000256 / 0|80 / 48|

Each pair checks a division (except 16 zero divisors) and an absolute-value
square root: 200496 host / 2000496 device MPFR result comparisons per PB.
The exact-path counts exclude the independent midpoint-resolver tests.
`compute-sanitizer --tool memcheck --error-exitcode 1 ...cu_fround_device_test
1000` reports zero errors at all six precisions.

## MPLAPACK integration / reproducibility

The mpc_cuda ExternalProject actually lives in the top-level
`CMakeLists.txt`, not under `cmake/`. Both autotools and CMake apply this
patch immediately after extraction;
patch failure stops the build. Autotools' extract stamp depends on the
patch; CMake has an independent `patch_input` step before download/extract,
so a changed/touched patch invalidates extraction, patching, configure,
build and install, rather than trying to patch an already patched tree.
The C locale and sm_80/sm_90 / glibc feature flags are retained.

Fixed real Rtrsm now calls `cu_fdiv_rn`. Runtime precisions and
`MPLAPACK_MPFR_CUDA_FORCE_RUNTIME=1` retain the old cu_mpfr engine. Keep the
existing default-exponent-range gate; detect every intermediate range or
domain error and return false without committing output. For Rtrsm alone,
its arithmetic engine retains MPFR signed-zero semantics in the unused sign
field of zero fixed values (the new upstream functions still produce +0).
Other real fixed BLAS, complex BLAS and dd/td/qd arithmetic are unchanged.
No `mpblas/optimized/common/` header is changed. Only fixed Rtrsm launches
request a 64 KB stack for Newton temporaries.

CMake tests (device absence skips with 77):

```sh
ctest --test-dir BUILD --output-on-failure \
  -R 'mpfr_cuda_(blas3|Rblas|Cblas|Rgemm|cu_fround)_'
compute-sanitizer --tool memcheck --error-exitcode 1 \
  BUILD/mpfr_cuda_cu_fround_device_test 1000
MPLAPACK_MPFR_CUDA_FORCE_RUNTIME=0 compute-sanitizer --tool memcheck --error-exitcode 1 \
  BUILD/mpfr_cuda_blas3_device_test --fixed-only
# Memory/stack checking across both engines (see the legacy API note below):
compute-sanitizer --tool memcheck --report-api-errors no --error-exitcode 1 \
  BUILD/mpfr_cuda_blas3_device_test
```

Standalone host oracle:

```sh
c++ -O2 -std=c++17 -I/path/to/patched/mpc_cuda/include \
  mpblas/optimized/mpfr/cuda/test/cu_fround_test.cpp -lmpfr -lgmp -o cu_fround_test
./cu_fround_test 100000
```

Benchmarks are built with the CMake MPFR CUDA tests, not run by ctest:

```sh
BUILD/mpfr_cuda_cu_fround_bench 65536
MPLAPACK_MPFR_CUDA_FORCE_RUNTIME=0 BUILD/mpfr_cuda_Rtrsm_round_bench 512 1024 L U N
MPLAPACK_MPFR_CUDA_FORCE_RUNTIME=1 BUILD/mpfr_cuda_Rtrsm_round_bench 512 1024 L U N
```

The operation benchmark reports **kernel-throughput ns/operation**, not
single-thread latency, for approximate fixed / correctly rounded fixed /
cu_mpfr div and sqrt; warm up then five rounds alternating operation order.
Rtrsm includes packing, upload, kernel, download, unpack; initialize full-PB
off-diagonal/RHS significands outside timing, one warmup, then one measured
call per process (A has diagonal exactly 2, keeping the system well conditioned).
The runtime selection is cached, so compare in separate processes: five
rounds, alternating new/old then old/new, and take medians. Host sharing
makes the absolute times non-calibrated; never report estimated timings.

## Build / patch rebuild verification

PASS: fresh source extraction and missing install prefixes, both build
systems, CUDA 13.1, **sm_80 and sm_90**, followed by touching the final patch
and rebuilding both a second time. Final installed headers match the patch
byte for byte; both libmpc_cuda.a archives pass `check_unique_symbols.sh`.
The previous install was retained only as a `/tmp` measurement baseline,
outside the build paths. Fresh `patch -p1 --dry-run` passes against an
independent pristine tar extraction. SHA256 remains
`bcac28655a62d381ffe38fc09704aca39f45610d1f697c351b4fc0acaa2278e7`;
`git diff -- external/mpc_cuda/download external/distfiles.sha256` is empty.

Actual patch-touch rebuild output (abridged, both completed):

```text
# autotools
rm -rf work ; mkdir -p work
patching file Makefile.in
patching file include/mpc_cuda.cuh
patching file include/mpc_cuda/cu_fround.cuh
env MAKEFLAGS= LC_ALL=C CUDA_LIBDIR=/usr/local/cuda/lib64 make -j1
>> compile mini-gmp
>> compile MPFR objects
>> compile MPC objects
mpc_cuda 0.0.1 installed under .../external/mpc_cuda/../i/MPC_CUDA

# CMake ExternalProject
No patch_input step for 'mpc_cuda_bundled'
Performing download step (verify and extract) for 'mpc_cuda_bundled'
-- extracting... [tar xf]
-- extracting... [clean up]
Performing patch step for 'mpc_cuda_bundled'
patching file Makefile.in
patching file include/mpc_cuda.cuh
patching file include/mpc_cuda/cu_fround.cuh
Performing build step for 'mpc_cuda_bundled'
>> compile mini-gmp
>> compile MPFR objects
>> compile MPC objects
mpc_cuda 0.0.1 installed under .../external/i/MPC_CUDA
Completed 'mpc_cuda_bundled'
```

## Library and ABI checks

PASS: autotools shared/static MPLAPACK library and CMake static library,
both sm_80/sm_90. All six archives checked (both libmpc_cuda.a, both
libmplapack_mpfr_opt_cuda.a, CMake libmplapack_mpfr_opt.a and
libmplapack_mpfr.a) have no duplicate strong global definitions.

`nm -D --defined-only` on the autotools shared library: **3983 symbols
before and after**. Excluding the 340 `__cudaRegisterLinkedBinary_*` and
340 `__fatbinwrap_*` compiler-generated identifiers (whose build hashes/PIDs
naturally change), the **3303 symbol names match exactly**. There are no
new public MPLAPACK symbols or entry points. The common headers are unchanged.

Full CMake regression run: **18/18 passed, zero skipped**, 8602.90 seconds.
This includes both engines of Rgemm and blas3, host/device Rblas and Cblas,
fixed Cblas, Winograd and the new cu_fround oracle tests. Rtrsm comparisons
remain strict even for signed zeros; no tolerance was weakened.

Sanitizer qualification: the full, unfiltered blas3 run numerically passed
all 524672 comparisons, but `--report-api-errors explicit` exits 1 with
**two handled API diagnostics**: `cudaDeviceSetLimit` and the subsequent
`cudaGetLastError` in the unchanged `limit_guard` destructor of
`Rgemm_rt_cuda.cu`. CUDA disallows restoring the smaller malloc heap limit
after device malloc was used; that guard already documents and clears this
error. The retained **pre-change shared library reproduces the same two
diagnostics** with a runtime-only Rtrsm call. This is not an invalid-memory
or RN-stack report, and the old runtime implementation is deliberately not
changed by 07a.

`blas3_cuda_test --fixed-only` selects just the five fixed precisions
(default/no-argument behavior is unchanged), so the new RN kernels can be
checked with all API reporting enabled without executing that legacy heap
cleanup. The full-program memory/stack check uses `--report-api-errors no`
only to avoid counting the handled legacy API return as a memory failure.
Final focused results (PASS):

- cu_fround device memcheck, full API reporting, all six PB: **0 errors**
  (1256 pairs / 2496 MPFR comparisons each, exact div/sqrt counts 80/48).
- fixed-only blas3 memcheck, full API reporting: **327920 elements,
  0 failures, 0 sanitizer errors**, all five fixed precisions and every
  side/uplo/trans/diag; signed-zero comparisons for Rtrsm stay strict.
- full blas3 memory/stack check, API-return reporting disabled as explained
  above: **524672 elements, 0 failures, 0 sanitizer errors**.
- `CUDA_VISIBLE_DEVICES='' ...cu_fround_device_test`: exit **77** as required.
- After adding the test-only selector, host default/runtime tests were
  rechecked (2/2 PASS); the default device path was also checked by the full
  memory/stack run. No library numerical code changed after the 18-test run.

## Rtrsm timings (measured, not estimates)

A100 80GB PCIe / CUDA 13.1; shared Threadripper host. Every cell is the
median of **five interleaved rounds**, reversing old/new order each round.
These are end-to-end entry-point seconds, with one warmup per process.
Flags are side/uplo/trans, diag=N, alpha=-1/2, m=n. Nonzero off-diagonal A
and RHS entries have full-PB significands; every A diagonal is exactly 2.
Old is cu_mpfr, new is cu_freal + cu_fdiv_rn. No CPU fallback: every measured
entry point returned true. Absolute times are not calibrated for host load.

### Self-contained CMake library (acceptance cases)

| PB | n | side/uplo/trans | old cu_mpfr (s) | new cu_freal (s) | old/new |
|---:|---:|:---:|---:|---:|---:|
|512|512|LUN|6.684604|0.964564|6.93|
|512|512|RLT|6.197510|0.892115|6.95|
|512|1024|LUN|26.616132|4.039389|6.59|
|512|1024|RLT|24.934635|3.795860|6.57|

Reproduce with `python3 mpblas/optimized/mpfr/cuda/test/run_round_bench.py
BUILD 5 acceptance`.

### Wider precision/size sweep

This sweep isolates the new bridge and fixed kernel translation unit
(sm_80, -O3, -fmad=false) and links the retained pre-change shared library
for the old cu_mpfr engine. It uses exactly the same Rtrsm entry-point
source and inputs, not a different solve algorithm. The CMake-library
measurements above independently confirm the required 512-bit cases.

| PB | n | side/uplo/trans | old cu_mpfr (s) | new cu_freal (s) | old/new |
|---:|---:|:---:|---:|---:|---:|
|256|256|LUN|1.007908|0.126803|7.95|
|256|256|RLT|0.953827|0.118358|8.06|
|256|1024|LUN|15.488797|2.129618|7.27|
|256|1024|RLT|14.817153|1.994794|7.43|
|512|256|LUN|1.712389|0.249919|6.85|
|512|256|RLT|1.594273|0.235774|6.76|
|512|512|LUN|6.684362|0.962640|6.94|
|512|512|RLT|6.182458|0.876348|7.05|
|512|1024|LUN|26.661375|4.029639|6.62|
|512|1024|RLT|24.935813|3.843811|6.49|
|512|2048|LUN|106.908493|15.886198|6.73|
|512|2048|RLT|99.843389|15.047063|6.64|
|1024|256|LUN|3.737435|0.661448|5.65|
|1024|256|RLT|3.405729|0.668427|5.10|
|1024|1024|LUN|59.661911|9.422702|6.33|
|1024|1024|RLT|53.752054|9.507042|5.65|

Raw local records: `bench_local/results/goal07a_trsm_interleaved.csv` and
`goal07a_trsm_cmake_interleaved.csv` (ignored working-tree assets).

## Standalone operations (measured kernel throughput, ns/operation)

CMake device benchmark, 65536 random positive input pairs in [1,2), all
significand bits filled; median of five rounds, alternating operation order.
Excludes allocation/transfers; includes each kernel's input views and output
store. **Not single-thread latency.** Approximate cu_fmath uses its unchanged
128 guard bits; the new RN API uses 64 plus the exact-boundary filter.

| PB | cu_fdiv | cu_fdiv_rn | cu_mpfr_div | cu_sqrt | cu_sqrt_rn | cu_mpfr_sqrt |
|---:|---:|---:|---:|---:|---:|---:|
|256|2.922|2.391|5.672|3.969|4.047|40.328|
|512|7.406|6.625|11.672|14.312|11.984|79.313|
|1024|41.031|37.984|32.406|58.438|53.469|209.703|

At 1024 bits, standalone cu_mpfr division is faster than the RN fixed
division in this throughput test; the whole fixed Rtrsm still wins because
its much more numerous multiply/subtract updates use cu_freal too. No
claim is made that every fixed primitive is faster at every precision.
Raw local output: `bench_local/results/goal07a_operations.txt`.

UNVERIFIED: execution on **sm_90 hardware** (only an A100/sm_80 is available;
both build systems successfully generated sm_80 and sm_90 code). All
requested A100 numerical/stack checks are verified, with the handled legacy
runtime API diagnostic explicitly qualified above. Goal 07 can now measure
cu_freal panels against cu_mpfr and host panels.
