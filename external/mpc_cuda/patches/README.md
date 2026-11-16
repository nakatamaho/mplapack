# mpc_cuda patches / proposed upstream PR

## Canonical cu_freal zero results (goal 12)

Status: **not yet submitted** (PR in preparation by the maintainer). Verified
against v0.0.2-dev: the cu_fzero oracle passes at every PB on host and A100
(all combinations of +/-0 inputs, poisoned zero limbs, random operands,
+,-,*, neg/abs/scale2, conversions: 298424 MPFR comparisons and 287098
frozen-original comparisons per execution mode, 0 failures). No upstream PR
has been made. This section does not change the submission status of the
division/sqrt patch below.

`patch-cu_freal-zero-sign` applies to **v0.0.2-dev 0f219b0** (rebased from the
abd9407 version; the same five minimal edits, since v0.0.2-dev kept the old
zero handling in cu_fadd/cu_fsub/cu_neg/cu_abs/cu_scale2), independently of
PR #1; MPLAPACK applies it after the division/sqrt patch. Both extraction
stamps depend on both patches. The tarball changes only with the goal-12
bundled update. The standalone tests are carried by MPLAPACK
(`test/cu_fzero_test.cpp` with the frozen abd9407 oracle
`test/cu_freal_before_zero.cuh`, host and device); the CUDA test returns 77
without a GPU.

Reproducer: `cu_freal<256> z(0.0); auto r=z-z;` formerly returned a zero with
`sign=-1`. `cu_fsub` negated the zero subtrahend, and `cu_fadd` returned it
verbatim. The patch handles both-zero addition before returning an operand,
avoids negating a zero subtrahend, and canonicalizes zero results of
`cu_neg`, `cu_abs` and `cu_scale2`. It retains all signatures and nonzero paths.
Existing multiplication, numeric/precision conversions and the RN div/sqrt
APIs already canonicalize zero. Copying the struct is still a bitwise copy;
this patch does not redefine copy semantics or special-function algorithms.

The oracle compares nonzero results with MPFR RNDN and the frozen abd9407
implementation, field by field. For either MPFR zero sign, the required result
is canonical +0: zero sentinel, sign=1, all limbs zero. Inputs include both
zero signs, poisoned zero limbs, powers of two, both signs of full-width random
significands, and every pair in the test corpus. Helpers, double conversion and
precision round-trips are also checked. PB=224,256,512,768,992,1024,2048:
**42,632 MPFR comparisons and 41,014 nonzero old-result comparisons per PB,
zero failures on host and A100**, for each execution mode. CMake focused tests
2/2 PASS; autotools default tests 6/6 PASS. Full CMake MPFR CUDA regression:
**84/84 PASS**, no skips, 67963.31 seconds, including complete forced-runtime
768/1024/2048 solver tests. The initial failures were independently reproduced
with pre-patch headers and repaired separately in a3c422820 (runtime HERK and
its right-looking column caller). Touching this patch re-extracts/repatches/
rebuilds/installs mpc_cuda in both build systems. Both installed headers match.
Six relevant archives have unique symbols; stable dynamic exports remain 3346.

MPLAPACK's fixed Rtrsm and resident real/complex factor policies save operand
zero signs before invoking arithmetic and restore signed zero afterward. Their
strict signed-zero tests are not weakened. Generic public fixed BLAS packing
canonicalizes zeros. The goal10 public fixed-complex corpus remains at 45996
checks / 1096 existing zero-sign exceptions, and quick TRSM at 504 checks /
22194 existing exceptions, zero failures. CMake Rgemm has 6 existing sign-only
zero differences; fixed Cblas has 1096 and fixed/runtime blas3 has 580/0.
No new exception or tolerance was added. Runtime768 solver --quick memcheck
passes 2761 checks / 1779332 comparisons with zero signed-zero differences and
zero sanitizer errors; this focused pass does not lift the general upstream
complex memcheck qualification. Raw logs: bench_local/results/goal12/.

### Proposed upstream PR text (independent of PR #1)

**Title:** Canonicalize zero results of cu_freal arithmetic and helpers

cu_freal documents a single zero representation. Currently `(+0)-(+0)` creates
a flagged zero with sign=-1, and neg/abs/scale2 can propagate noncanonical zero
inputs. Handle the both-zero add case explicitly, skip negation of a zero
subtrahend, and use set_zero in the unary helpers when the input is zero.
Nonzero arithmetic and API signatures are unchanged. Multiplication and
conversions already use the canonical representation. This patch is independent
of the correctly-rounded division/square-root addition.

The included standalone host/CUDA test compares the sign, exponent and every
limb against MPFR RNDN, mapping MPFR's signed zeros to cu_freal's documented +0,
and checks nonzero results against a frozen original implementation. It covers
seven precisions including two with 32 unused bits, zero signs and poisoned
zero limbs, full-width operands, binary operators and unary/conversion helpers.
Applications that require MPFR signed-zero semantics must retain a separate
sign policy; cu_freal itself has no signed-zero contract.

## Correctly rounded division/sqrt patch

Status: **submitted upstream** as https://github.com/tkouya/mpc_cuda/pull/1
(2026-10-30; the upstream branch also adds tools/test_fround.cu to
`make check-fixed`, and documents the API in doc/manual.md and manual.ja.md).
Base: upstream https://github.com/tkouya/mpc_cuda branch `v0.0.2-dev`,
commit **0f219b0** (goal 12 update; the bundled distfile is
`mpc_cuda-0.0.2-dev-0f219b0.tar.xz`, `git archive` of that commit plus the
generated `configure`). The old abd9407 0.0.1 distfile and its SHA256 are
kept for reference. Both patches were rebased to the v0.0.2-dev tree; remove
a patch when the bundled upstream revision contains it.
This branch and 3.1 use the same patch/commit; no separate backport branch.
Upstream checks in the patched tree (goal 12, sm_80): `make check-fixed`,
`check-fixed-gpu`, `check-fmath` all pass.

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
