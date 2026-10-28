# GOAL 07a: correctly rounded division and square root for cu_freal (mpc_cuda patch)

Prerequisites: the MPFR CUDA fixes of b626121f (mutex around the MPFR GPU
entry points, device stack limit for the complex kernels, Rswap/Rrot);
goals 01-06 are done. Goal 07 depends on this goal. Read 00-overview.md first, and
mpblas/optimized/mpfr/cuda/README.md.

## Background (verified on the bundled mpc_cuda-0.0.1-abd9407)

- cu_fp::cu_freal<PB> (include/mpc_cuda/cu_freal.cuh): fixed precision PB
  (multiple of 32), significand in registers, value = sign · M ·
  2^(exp − N·64), N = ceil(PB/64) limbs, M left-justified, zero is a flag
  (exp == CU_FREAL_ZERO_EXP). +, −, × are bit-exact with MPFR RNDN.
- include/mpc_cuda/cu_fmath.cuh says in its header: "high-accuracy, NOT
  correctly-rounded". cu_fdiv is a Newton reciprocal followed by a
  multiplication, cu_sqrt a Newton reciprocal square root at a wider
  working precision W, then converted to PB: near a rounding boundary they
  can round the wrong way, so they are not mpfr_div / mpfr_sqrt bit for
  bit. sqrt of a negative number returns 0; there is no Inf or NaN.
- The precisions MPLAPACK uses on cu_freal (256, 512, 768, 1024, 2048,
  MPLAPACK_MPFR_CUDA_FIXED_PRECISIONS in Rgemm_kernel_cuda.h) are multiples
  of 64 (no 32 unused low bits); the functions must still handle SB = 32.
- Consequence today: MPLAPACK's MPFR Rtrsm always runs on cu_mpfr (the
  device port of MPFR, slow), and goal 07 cannot use cu_freal for the
  divisions and square roots of the getrf / potrf panels.

## Policy: a patch kept in MPLAPACK, offered upstream

Upstream (https://github.com/tkouya/mpc_cuda) may adopt this late or never.
The change is therefore kept as a patch file in MPLAPACK, applied to the
bundled tarball, and the same file is the upstream submission:

- `external/mpc_cuda/patches/patch-cu_freal-correctly-rounded-div-sqrt`:
  a unified diff (git format-patch style, with a header that explains the
  change, the algorithm and the tests) against the tree of
  mpc_cuda-0.0.1-abd9407 as extracted from
  external/mpc_cuda/download/mpc_cuda-0.0.1-abd9407.tar.xz. Do not change
  the tarball or its SHA256 in external/distfiles.sha256.
- Applied at extract time by both build systems, as external/lapack does
  with its patches/<version>/patch-* files:
  external/mpc_cuda/Makefile.am (after the tar extraction, before
  configure; with a stamp) and the CMake ExternalProject mpc_cuda_bundled
  (PATCH_COMMAND). A failing patch must stop the build.
- The patch only ADDS functions; the existing cu_fdiv / cu_sqrt and
  operator/ keep their behaviour, so the patch applies and builds for any
  user of mpc_cuda. (If upstream prefers to make the operators correctly
  rounded, that is their choice; MPLAPACK calls the new names.)
- Upstream submission material: `external/mpc_cuda/patches/README.md`
  (what each patch does, the upstream commit it is made against, its
  status: "not yet submitted / submitted <url> / merged in <commit>";
  remove the patch once the bundled tarball contains it) and the PR
  description text in the commit body. Opening the upstream PR needs the
  maintainer's account: prepare the text, do not attempt it.
- Branch 3.1 now equals claude/mpfr-cuda-rgemm (fast-forwarded; it builds
  external/mpc_cuda and the MPFR CUDA library), so the patch is carried on
  3.1 simply by the commit of this goal; no separate 3.1 branch is needed.
- Rebuilds: the stamp of the extract/patch step must depend on the patch
  file(s), so that editing a patch re-extracts, re-patches and rebuilds
  mpc_cuda in both build systems.

## Tasks

1. In the patch, new header include/mpc_cuda/cu_fround.cuh (included by
   mpc_cuda.cuh and installed with the other headers — check mpc_cuda's
   install rule), namespace cu_fp:

       template<int PB> __host__ __device__ int
       cu_fdiv_rn (cu_freal<PB> &q, const cu_freal<PB> &a, const cu_freal<PB> &b);
       template<int PB> __host__ __device__ int
       cu_sqrt_rn (cu_freal<PB> &r, const cu_freal<PB> &a);

   Result: q = RNDN(a / b), r = RNDN(sqrt(a)) exactly as mpfr_div /
   mpfr_sqrt with MPFR_RNDN at precision PB (ties to even), with an
   unbounded exponent range. Return value: 0 on success; nonzero, output
   unspecified, for b = 0 or a < 0 (the caller handles these: MPLAPACK
   falls back to the CPU). a = 0 gives +0 (sign rules of MPFR for −0 do not
   apply: cu_freal has no signed zero; document it). The returned int may
   also carry the MPFR ternary value (sign of rounded − exact) if cheap.
2. Algorithm (Ziv-style, documented in the header):
   - fast path: Newton reciprocal / reciprocal square root at W = PB + 64
     bits (reuse cu_fmath's iterations at the wider type), giving an
     approximation y with a proven error bound below a few ulps of W;
   - if y rounded to PB is not ambiguous (y is farther than the error bound
     from every rounding boundary of PB), return it;
   - otherwise decide exactly: let t be the PB+1-bit midpoint candidate;
     division: sign of a − t·b computed exactly (t·b is a 2N+1-limb product,
     compare with a shifted, no rounding); sqrt: compare t² with a exactly.
     Exact results (remainder 0) and ties (to even) handled here.
   - only integer limb operations in the exact step (cu_umul, adds, shifts
     over the limb arrays); no floating point there.
   Write the error bound of the Newton step into the header with its
   derivation; the tests below must hit the exact path.
3. Tests, in MPLAPACK (not in the patch): mpblas/optimized/mpfr/cuda/test/
   cu_fround_test.cpp, built as a device test (CMake, SKIP_RETURN_CODE 77)
   and a host build (the header is host-usable C++), comparing with host
   MPFR (mpfr_div / mpfr_sqrt, RNDN, same precision) bit for bit for PB in
   {256, 512, 768, 1024, 2048} and one PB with SB = 32 (e.g. 224 or 480):
   - random operands over a wide exponent range (≥ 10^6 per PB on the
     device, ≥ 10^5 on the host);
   - hard cases constructed on purpose: exact quotients and squares
     (a = q·b, a = r² with PB-bit q, r), quotients and roots exactly at a
     midpoint (ties to even), and within 1 ulp(W) of a midpoint on both
     sides (build them from a chosen midpoint t: b random, a = t·b ± tiny);
   - powers of two, operands with all significand bits set, b = 1, a = b;
   - report how many cases took the exact path (must be > 0 per PB).
4. Use in MPLAPACK (same goal, real Rtrsm only; the complex routines divide
   with libmpc_cuda's cu_mpc_div and are unchanged): the MPFR blas3 Rtrsm on
   the fixed precisions moves from cu_mpfr to cu_freal with
   cu_fdiv_rn (mpblas/optimized/mpfr/cuda/blas3_kernel_cuda.h /
   blas3_bridge_cuda.cpp), keeping the result bit-identical to the
   reference Rtrsm (existing mpfr_cuda_blas3 tests); update the README
   ("Rtrsm always uses cu_mpfr" no longer holds for fixed precisions) and
   keep the cu_mpfr path for the runtime precisions. Keep the exponent
   range gate Rtrsm has today. The Newton steps at W = PB + 64 bits may need
   more device stack than the real fixed kernels have today (the default
   1 KB; the complex kernels already raise it, see mpfr_cblas_stack in
   Cblas_mpfr_device_cuda.cu): check with compute-sanitizer and raise the
   limit for the kernels that call cu_fdiv_rn / cu_sqrt_rn if needed.
5. 00-overview.md and goal 07: state that cu_fdiv_rn / cu_sqrt_rn exist, so
   the getrf / potrf panels can stay on cu_freal (goal 07 still measures).

## Acceptance criteria

- Clean tree check: delete external/i/MPC_CUDA and external/mpc_cuda/work,
  reconfig, build — the patch applies (autotools and CMake), mpc_cuda builds
  (C locale, sm_80 and sm_90), `git diff` of the tarball and
  distfiles.sha256 is empty.
- cu_fround_test: 0 mismatches against MPFR for every PB, on the A100 and
  on the host; exact-path counts reported per PB.
- mpfr_cuda_blas3_device and _runtime, mpfr_cuda_Rblas_device,
  mpfr_cuda_Rgemm_device pass (bit identity unchanged); Rtrsm timings on
  the fixed path before / after for 512 bits, m = n = 512 and 1024, side L
  and R, in the commit body.
- `patch -p1 --dry-run` of the patch file applies to a fresh extraction of
  the tarball outside the build (the upstream-submission check).
- misc/check_unique_symbols.sh on libmplapack_mpfr_opt_cuda.
- Editing the patch file (e.g. touching it) re-applies it and rebuilds
  mpc_cuda (show the make output of both build systems).
- One commit on claude/mpfr-cuda-rgemm: "mpc_cuda: correctly rounded cu_freal division and sqrt (factorization goal 07a)". No push.
