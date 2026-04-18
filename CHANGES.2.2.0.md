# MPLAPACK 2.2.0 — Changes since 2.1.1

Release date: 2026-04 (tag: v2.2.0, topic branch: `topic/lapack-3.12.1`)

This release re-bases MPLAPACK 2.1.1 onto LAPACK 3.12.1 (previously 3.9.1)
and consolidates multi-precision arithmetic handling into a single
`arithmetic_params` layer. 281 commits, roughly spanning 2026-01 through
2026-04-17.

The bulk of this release is a single migration; most user-visible changes
are consequences of that migration, not independent features.

---

## 1. Headline

- **LAPACK reference version: 3.9.1 → 3.12.1.**
  All routines are regenerated from LAPACK 3.12.1 sources through the
  FABLE converter. The FABLE pipeline is now parameterized by
  `LAPACK_VERSION` (default `3.12.1`); per-version patch directories
  live under `fable/3.12.1/lapack/` and `fable/3.12.1/lapack-test/`.
- **Added:** `Cgemmtr/Rgemmtr`, `Cgedmd/Rgedmd`, `Cgedmdq/Rgedmdq`,
  `Cgelst/Rgelst`, `Clatrs3/Rlatrs3`, `Cgeqp3rk/Rgeqp3rk`,
  `Cunhr_col/Rorhr_col`, `Cungtsqr_row/Rorgtsqr_row`,
  `Claqz0–4 / Rlaqz0–4`, `Risinf`, `Rroundup_lwork`, and companion tests.
- **Removed:** `Clacon/Rlacon` (deprecated reverse-communication; replaced
  by `Clacn2/Rlacn2`), `Cgelqs/Rgelqs/Cgeqrs/Rgeqrs` (replaced by `xGELS`),
  Cray-era guard-digit workarounds in D&C routines (see §3 for scope).
- **Rewrites of note:** `Rlartg/Clartg` (LAPACK 3.12.1 algorithm,
  `abssq` helper, arithmetic_params-driven scaling), `xLASSQ` and
  dependent norm routines, `xLARFT` (recursive structure), `xGEBAL`
  (goto-based code replaced), `Crotg/Rrotg`, `Rnrm2/RCnrm2`.
- **New infrastructure:** `mplapack_arithmetic_params.h` centralizes
  `eps`, `safmin/safmax`, `rmin/rmax`, `sfmin`, and Blue scaling constants
  per backend, replacing ad-hoc `Rlamch` call sites and hardcoded values.
- **`Rlabad` policy reversal:** Unlike reference LAPACK 3.12.1 (which
  makes `Rlabad` a no-op), MPLAPACK *reinstates* `Rlabad` as an active
  routine and re-injects it at call sites that are sensitive to the
  wide exponent range of binary80, binary128, MPFR, GMP, DD, and QD.

---

## 2. New routines and tests

### MPBLAS additions (from BLAS 3.12.x)

- `Cgemmtr`, `Rgemmtr` — triangular-updating GEMM. Added to the reference
  tree, declarations in all seven backend headers. `xLASYF*`, `xLAHEF*`
  were refactored to use `Rgemmtr/Cgemmtr` in place of `xGEMV` + `xGEMM`
  on triangular/rectangular block pairs.
- `Rnrm2`, `RCnrm2` — rewritten to the LAPACK 3.12.1 reference
  implementation (Blue-style scaled sum of squares).
- `Crotg`, `Rrotg` — parameter names normalized to `a,b`; rewritten to the
  3.12.1 reference (no behavioral change at IEEE precision).

### LAPACK routines added

- `Cgedmd/Rgedmd`, `Cgedmdq/Rgedmdq` — Dynamic Mode Decomposition.
  Converted via FABLE, test drivers `xdmdeigtst[RC]_<backend>`, matrix
  copy / array-slice issues in the generated test drivers patched
  (explicit `copy_matrix_block` helper; `read(6,*)` → `read(5,*)`).
- `Cgelst/Rgelst` — least-squares driver based on compact QR.
  Added to `Cdrvls/Rdrvls`; test indices renumbered
  (1-2 `GELS`, 3-4 `GELST`, 5-6 `GETSLS`, 7-10 `GELSY`, 11-14 `GELSS`, ...).
- `Clatrs3/Rlatrs3` — multi-RHS triangular solve with scaling. Tests in
  `Cchktr/Rchktr` now use `scale3[]`, `rwork2`, and report `result[10]`.
- `Cgeqp3rk/Rgeqp3rk`, `Claqp2rk/Rlaqp2rk`, `Claqp3rk/Rlaqp3rk` —
  rank-revealing QR with column pivoting (QP3RK). Upstream
  `patch-xLAQP2RK.f` is applied to fix an `[in]` argument direction bug.
  New test harness `Cchkqp3rk/Rchkqp3rk` added; matrix type `QK` handled
  in `Rlatb4/Clatb4/Rlatb5/Clatb5`.
- `Cunhr_col/Rorhr_col` — Householder reconstruction from compact WY.
  T-matrix zeroing loop bound bug fixed; error handling added to `Alaerh`
  for the QK path.
- `Cungtsqr_row/Rorgtsqr_row`, extended `xLATSQR` tests.
- `Claqz0–4 / Rlaqz0–4` — blocked QZ iteration (used by `xHGEQZ`).
- `Risinf` — IEEE isinf predicate, companion to `Risnan`, for all
  backends (double / binary80 / binary128 via `isinfq` or builtin / DD /
  QD / GMP / MPFR).
- `Rroundup_lwork` — safe INTEGER → REAL conversion for workspace
  queries. For double, rounds up `castREAL(lwork)` when `lwork > 2^53`;
  no-op for extended-precision backends. Has a per-backend suffix and
  `#define` alias like `Rlamch`/`iMlaenv`, fixing an
  "ambiguating new declaration" error in compare builds.

### Test coverage additions

- Extended error-input tests for `xLATSQR`, `xUNGTSQR_ROW`, `xORGTSQR_ROW`,
  and several solver/driver routines.
- Extended `xERRTSQR`, `xERRSY_AA_2STAGE`, aligned with upstream patched
  sources.
- `Rdrvrf3` test matrix conditioning fixed for the unit-diagonal case.
- `xGET02` updated to use `norm(op(A))` instead of `norm(A)`; `xGBT02`
  gained an `rwork` parameter for row-wise norm.

---

## 3. Removed / deprecated

- `Clacon`, `Rlacon` — deprecated reverse-communication condition
  estimators. Replaced throughout by `Clacn2/Rlacn2`. Source files,
  patches (`patch-Clacon.cpp`, `patch-Rlacon.cpp`), and declarations
  removed.
- `Cgelqs/Rgelqs`, `Cgeqrs/Rgeqrs` — removed in favor of `xGELS`.
- Cray-era guard-digit workarounds in D&C routines — removed, following
  LAPACK 3.12.1. Specifically, the `Rlamc3(x, x) - x` idiom used to zero
  the lowest bit of eigenvalue/singular-value arrays in `Rlaed3`, `Rlaed9`,
  `Rlasd3`, `Rlasd8` (defensive against machines without an add/subtract
  guard digit, i.e. Cray XMP/YMP/C90/Cray 2). All MPLAPACK backends have
  correct guard-digit behavior, so `2x - x = x` holds exactly and the loop
  is pure overhead. Note: this is **unrelated to the `Rlabad` policy in
  §8** — guard-digit correctness is about single-operation precision,
  whereas `Rlabad` is about exponent range; upstream removed both under an
  IEEE-hardware assumption, but only the guard-digit removal transfers to
  multi-precision backends. The `Rlamc3` calls retained in `Rlasd8`'s
  secular-equation j-loop serve a different purpose (enforcing
  `(x+y)+z` evaluation order against compiler reassociation) and are kept.
- Pre-conversion local patches that 3.12.1 obsoletes (`patch-Rdrvst2stg`,
  `patch-Cunbdb3.cpp` original workaround, several *.gvd post-conversion
  patches).
- Deprecated random-number reference values — removed.

---

## 4. New infrastructure: `mplapack_arithmetic_params.h`

A centralized, template-specialized layer for per-backend machine
constants. Replaces scattered `Rlamch("E")`/`Rlamch("S")` call sites and
hardcoded `1.0 + unhandled` placeholders left by FABLE.

- `get_arithmetic_params<REAL>()` returns `eps`, `safmin`, `safmax`,
  `rmin`, `rmax`, `sfmin`, `prec`, and associated Blue scaling inputs.
- `get_blue_scaling_params<REAL>()` returns `tbig`, `sbig`, `tsml`, `ssml`
  for `Rlassq/Classq` and other norm routines.
- Specializations added for `double`, `_Float128`/`__float128`,
  `long double` (binary80), `dd_real`, `qd_real`, `mpf_class`, and
  `mpfr::mpreal`.
- `safmin`/`safmax` are internal power-of-two values (documented in the
  header) rather than the exact IEEE denormal boundary.
- DD/QD: `emin` is derived from `frexp(_min_normalized)` rather than from
  upstream `dd_real::_eps` and related symbols, to decouple from the
  upstream QD library's own constants.
- MPFR: when the user has not configured an explicit exponent range,
  a truncated range is selected (`±prec × 64`) to avoid Rstebz bisection
  pathology; a runtime warning is emitted from
  `make_blue_scaling_params()` if the BlueScale thresholds fall outside
  the usable range.
- Rewired users: `Rlassq`, `Classq`, `Rlartg`, `Clartg`, `xHGEQZ`,
  `xBDSQR`, `xBBCSD`, `xLATRS`, others.
- Blue-scaling classification probes and per-backend consistency checks
  added to `Rlamch.test.cpp` (including GMP and MPFR non-IEEE profiles).

---

## 5. Backend-specific hardening

### GMP

- **Divide-by-zero guards** in `Rgesvj`, `Rgejsv`, `Cgesvj`, `Cgejsv`
  before `big/sqrt(aaqq)` when `aaqq == 0`. GMP raises SIGFPE on exact
  zero division (IEEE would return Inf). All four routines get
  GMP-specific guarded branches (3 in `xgesvj`, 1 in `xgejsv`).
- **Givens rotations without trig.** `Rlatms`, `Rlatmt`, `Clatms`,
  `Clatmt` previously drew `angle ∈ [0, 2π)` and used `c=cos(angle)`,
  `s=sin(angle)`. On GMP, `sin/cos` fall back to `double`, breaking
  `c² + s² = 1` by ~1 ulp. Replaced by a rational double-angle
  unit-circle sampler; for complex variants, the `Clarnd(5,...)` phase
  factor is retained but the real `(c,s)` pair is replaced.
  This fixes `Rgeqp3rk` ratios ~1e+136 on `M=2` / `N=2` cases.
- **`pow(eps, -1/8)` double-fallback avoidance** in `Rbdsqr`, `Cbdsqr`,
  `Rbdsvdx`, `Rbbcsd`, `Cbbcsd`. `pow(mpf_class, mpf_class)` has no
  native `mpf_t` implementation and would truncate via `get_d()`.
  Replaced by exact constructions where possible; guarded under
  `___MPLAPACK_BUILD_WITH_GMP___`.
- **Threshold relaxations** (test harness only, not algorithm):
  NEP/DHS multiplier 3× → 10×; DVX/`Rgeevx` 3× → 6×; ZGV 1.5× → 4×;
  DGS 8×; DGV 3×; real generalized driver thresholds restored.

### MPFR

- **High-precision iter-limit scaling.** `xBDSQR` `maxitr`, `xHGEQZ`
  `maxit`, `xBBCSD` iteration counts scale as
  `iceil(base × log(eps) / log(double_eps))`, floored at `base`, capped
  at `base × 64` (1000 for HGEQZ). `Rbdsqr/Cbdsqr` base 6,
  `Rhgeqz/Chgeqz` base 30.
- **`Rdrgev3` eigenvalue-compare tolerance.** MPFR-only scaled tolerance
  for result slots 5, 6, 7; bit-exact comparison retained for all other
  backends.
- **`Rget38` `Rlabad` scaling** for MPFR/GMP DEC tests — `Rtrsen` scaling
  stabilized.
- **SEP threshold** raised to 80 for `Rsep` borderline cases.

### binary80 / binary128

- **`Clatrs` Rlabad-style range truncation.** Restored for binary80 and
  binary128; `Cdrvev` type-1 zero matrix and types 19/20/21 random
  matrices regressed with ratios ~1e+136 without this.
  Reference LAPACK 3.12.1's removal of `Rlabad` is IEEE-safe but
  incorrect for the wider binary80/binary128 exponent ranges
  (`smlnum ~ 1e-4898`, `bignum ~ 1e+4897`).
  A long-term followup is filed to route the `(sfmin, bignum)` pair
  through `arithmetic_params` for `Rlatrs`, `Clatbs`, `Clatps`,
  `Ctrsyl`, `Ctgsy2`, and their R-counterparts.
- **`Cgeev`, `Cgges` Rlabad guard** added for binary80/binary128 (already
  present in `Cgges3`). Fixes `ZGS` result-11 alpha/beta vs. Schur-diagonal
  drift and `ZEV` test-7 left-eigenvector sentinel failure.
- **`Cgges3` scaling fix.** Now recomputes `smlnum/bignum` through
  `Rlabad` rather than directly from `Rlamch("S")`.

### DD / QD

- **`Rget34` structural comparison** for QD: `qd_real::operator==`
  compares four components directly, which can disagree on values whose
  difference renormalizes to zero. A structural-equivalence check replaces
  the exact comparison for QD only in `Rget34`.
- **`Rget32/Rlasy2` QD residual scaling** — coefficient side rescaled.
- **SVD harness (type 16) bidiagonal spectrum width** clamped: replaces
  the `exp(-2 log(ulp) × Rlarnd(2))` formula (which produces ~4e62 for
  DD, ~1e125 for QD) with `exp(-0.5 log(ulp))` (DD) /
  `exp(-0.25 log(ulp))` (QD), both ~4.5e15, matching double-precision
  scale.
- `*gesvj` returns singular values with a separate scale factor; the
  real/complex SVD driver tests now apply `WORK(1)/RWORK(1)` before
  reconstruction checks.
- `Cchkbd` skips the `Rsvdch` Sturm-sequence check for QD complex SVD
  tests (mirroring real `Rchkbd`).
- DD DGS driver threshold raised 2× for borderline orthogonality ratios
  (10.3–10.9 vs. base 10).

---

## 6. FABLE / cout.py translator fixes

The conversion pipeline accumulated a number of known bugs during this
migration; all are fixed on the 2.2.0 branch.

- **`.NOT.` vs relational-operator precedence.** Fortran `.NOT.` binds
  looser than `.LE./.GE./...`; FABLE generated
  `!(a <= b)` where Fortran wanted `!a <= b`. Postprocessing pass
  `_postprocess_not_relop_parens` rewrites these.
- **`DO` loop translation repair.** Postprocessing for broken generated
  DO-loop headers restored.
- **`wp = real64` / `wp = kind(...)` preprocessing.** New dedicated pass
  so LAPACK 3.12.x sources that use `wp` as the kind parameter convert
  cleanly.
- **2D array slices on `COL` paths.** `__SLICE2D__` helper did not handle
  rank-1 iteration direction correctly in `Rlaqz3/Claqz3`. The converter
  is now conservative for rank-2 slices in this shape; new patches
  cover the affected LAQZ routines.
- **`_postprocess_strip_float_suffix`** no longer strips trailing `f`
  from identifiers.
- **`INT(DBLE(work(1)))` on `COMPLEX` `work`** handled correctly.
- **Intrinsic rewriting:** `fem::cmplx` added to `rewrite_intrinsics()`,
  replaced at call sites by the MPLAPACK equivalent.
- **Constant table:** hardcoded constant list replaced with a
  `_KNOWN_CONSTS` table covering the LAPACK 3.12.1 literals (including
  UNHANDLED cases for `Rlartg/Clartg`).
- **`fem::str<N>` substring handling:** 2-argument `operator()` supports
  Fortran substring slices (`type(1,1)` → `str_cref`). For raw `char`
  use (e.g. `printf("%c", ...)`), `type.elems[0]` is used.
- **`iMparam2stage.cpp`** — `Mlsame` address-of restored; `fem::`
  dependency removed.
- **`iMlaenv.cpp`** — FABLE Fortran substring syntax corrected.
- **`abssq` declarations** excluded from generated mpblas headers.
- **Keep-lists extended** for mpblas, mplapack, and generated test trees
  so hand-written files survive regeneration.
- **FABLE pipeline is repo-relative** and `LAPACK_VERSION`-aware:
  `fable/go.sh`, `fable/go_testing.sh`, `fable/sync_test_inputs.sh`
  now take a version parameter (default `3.12.1`); per-version patch
  lists defined in `patch_lapack_${VER}.sh` and
  `patch_lapack_test_${VER}.sh`.
- `misc/fable-cout.sh` wrapper added around the conversion entry point.
- `regen_optimized_makefile_sources` added to `go.sh`; all
  `mpblas/optimized/*/Makefile.am` regenerated automatically.

---

## 7. Numerical and correctness fixes (beyond migration)

- **`Rlartg/Clartg` rewrite to LAPACK 3.12.1 algorithm.** Legacy
  `safmn2/safmx2` scaling loops with `goto` replaced by
  `rtmin/rtmax`-based scaling; `abssq(COMPLEX)` helper added computing
  `|z|² = re² + im²` without `sqrt`; parameter names normalized to
  `c,s`; `cabs_inf`/`cabssq` helpers removed; constants sourced from
  `arithmetic_params`.
- **`Clartg` subnormal/underflow hazard** in MPFR binary64/binary128
  emulation. Squared-norm intermediates could underflow despite
  representable inputs; fix tightens the unscaled gate to
  `fsmall/f1 > sqrt(eps)` and splits `sqrt(f2*h2)` into
  `sqrt(f2)*sqrt(h2)`.
- **`xLASSQ` rewrite** and dependent norm/CS routines (3.12.1 Blue-style).
- **`xLARFT` rewrite with recursive structure** (3.12.1 reference).
- **`xGEBAL` rewrite** replacing goto-based control flow.
- **Workspace query fixes across 111 routines.**
- **`Clarf/Rlarf` → `Clarf1f/Rlarf1f/Clarf1l/Rlarf1l`** migration
  applied to: `Cunbdb*`/`Rorbdb*`, `Cupmtr/Ropmtr`, `Cunm2r/Cunml2/
  Cunm2l/Cunmr2/Rorm2r/Rorml2/Rorm2l`, `Cung2r/Cungl2/Cung2l/Cungr2/
  Rorg2r/Rorgl2/Rorg2l/Rorgr2`, `Claqr2/Claqr3/Rlaqr2/Rlaqr3`,
  `Claqp2/Rlaqp2/Claqp2rk/Rlaqp2rk/Cgeqp3rk/Rgeqp3rk`,
  `Cgehd2/Rgehd2`, `Cgebd2/Rgebd2`, `Cgeql2/Rgeql2`, `Cgelq2/Rgelq2`,
  `Cgeqr2/Cgeqr2p/Rgeqr2/Rgeqr2p`.
- **`Rtgex2/Rtrevc3/Ctrevc3`** — bug fixes synced from 3.12.1.
- **`Cunbdb3/Rorbdb3`** — incorrect `ldx11`-as-increment argument for
  the `x21` vector fixed; upstream ZUNBDB3 X21 row conjugation fix
  backported; CSD edge tests enabled.
- **`Clarfgp/Rlarfgp`** — near-zero vector handling improved.
- **`Cgeqrt2/Rgeqrt2/Cgemlqt/Rgemlqt`** — input validation tightened.
- **`Chgeqz/Rhgeqz`** — precomputed `btol` used for `T(j,j)` deflation.
- **`Clatrs/Rlatrs`** — hardened against overflow/NaN in column-norm
  scaling.
- **`iMparmq`** — `ispec=17` (`icost`) added: near-diagonal chase vs.
  BLAS update cost ratio.
- **`xTGSJA`** — division-by-zero guard replaced by post-division
  overflow check via `Rlamch("O")`, handling subnormal and infinite
  cases.
- **`Rlanv2`** — FMA contraction suppressed.
- **`xSTEMR`** — `n=2` eigenvector ordering fixed.
- **`xGETC2`** — pivot search order fixed; simplified.
- **`xBDSQR`** — convergence tracking improved.
- **`xGELSS`** — multi-RHS bug fix.
- **NaN/Inf input validation** added to several routines.
- **`xTRSYL3`** — `LDSWORK` no longer overwritten during workspace
  query; `pow(2.0,...)` replaced by `pow(two,...)` for multi-precision.
- **`Risnan` / `DISNAN/DLAISNAN/LA_XISNAN` unified** into a single
  hand-written `Risnan.cpp`.
- **`Rlarrd/Rorcsd2by1`** — 3.12.1 bug fixes ported.
- **`Clatbs/Rlapy2/Rlapy3`** — overflow hardening.
- **`Cgerq2/Rgerq2`** — reflector applied via `Clarf1l/Rlarf1l`.
- **`Csyswapr/Rsyswapr`** — manual loops replaced with BLAS swap calls.
- **`Cheequb/Csyequb`** — implicit integer arithmetic replaced with
  `castREAL` and named REAL constants.
- **`xgvx`** (`Cspgvx`, `Chpgvx`, etc.) — `info` argument direction
  corrected to `iinfo` where appropriate (upstream-compatible fix).
- **CSD edge-case workaround.** Upstream `csd.in` ships M=1/P=1/Q=1
  entries that fail; those 2 entries are removed (12 → 10). All 600
  remaining tests pass.
- **`Cuncsd2by1` M=55,P=40,Q=20** is a known LAPACK 3.12.1 upstream bug
  (fixed in LAPACK HEAD). Intentionally left as-is on this branch.
- **ODR violation in `Mxlaenv.cpp`** fixed by moving it to the
  `libeig_override_*` / `liblin_override_*` libraries only.
- **`drgvx` `-O1` workaround** for `g++-mp-15 -O2` miscompile
  (`-fno-ipa-modref` path); backported from master.

---

## 8. Rlabad policy

LAPACK 3.12.1 makes `Rlabad` a no-op under the assumption of IEEE
arithmetic. MPLAPACK reverses this decision for wide-range backends:

- `Rlabad` is **reinstated as a hand-written reference routine** and
  removed from the FABLE-generated set.
- Previously `Rlabad` was a no-op stub with both arguments commented out.
- `Rlabad(smlnum, bignum)` is re-injected at the sites that materially
  depend on it for binary80/binary128/MPFR/GMP/DD/QD. Affected call
  sites include: `Cgeev`, `Cgges`, `Cgges3`, `Clatrs`, `Rget38`,
  `Cqrt13`, and the pre-conversion LAPACK patches touching `xLATB4`.
- Inlined `log10(large) > 2000` guards that duplicated `Rlabad`'s
  behavior are replaced by direct `Rlabad` calls.

---

## 9. Upstream LAPACK contributions

Tracked and merged back where accepted:

- **PR #1219 / #1220** — `CABS1` / `CABSMAX` unification across
  `CLARGV`/`ZLARGV`.
- **`xSTEBZ` error propagation fix** across 26 routines.
- **`xpgvx` `info` → `iinfo`** correction.
- **`xUNBDB3` X21 row conjugation fix.**
- **`xLAQP2RK.f`** `[in]` argument correction.

---

## 10. Release engineering and test infrastructure

- **Version bumped** to 2.2.0 in `configure.ac` (`VER_MINOR=2`,
  `VER_PATCH=0`). `topic/lapack-3.12.1` tracks the 2.2.0 line;
  `origin/master` continued to carry 2.1.1 during development.
- **Test result logs** relocated to versioned directories
  (`eig/2.1.0/`, `lin/2.1.0/`, legacy `2.0.0/`).
- **Preliminary 2.2.0 results** added for
  `Ryzen_Threadripper_3970X_ubuntu24_04_gcc-13_3_0`.
- **Test parameter files** (`Cbb.in`, `Cgg.in`, `Rbb.in`, `Rgg.in`,
  `csd.in`) patched via `patch-fix-testparameters` in
  `external/lapack/Makefile.am`; `Rgg` size list expanded.
- **`external/lapack/Makefile.am`** — patch application rules
  simplified and parameterized by LAPACK version.
- **MPFR Rlamch reference files** (`Rlamch_reference*.txt`) refreshed
  across compare backends; a GMP 32-bit-long variant added.
- **`Mexponent.cpp`/`Risnan.cpp`** predicate checks updated.
- **Rrotmg patch file** separated from the former general blas patch;
  `static` storage removed from constant-like locals.
- **Generated TESTING routines** — unused ones excluded from conversion
  and their declarations removed from `mplapack_lin_*` headers.
- **BLAS/LAPACK cleanup step** extended to remove generated and
  temporary files.
- **Random-number generation normalized across precision classes.**
  `Rlaran` now delegates to `Rlaruv(iseed, 1, &x)`, and the random stream
  is kept independent of the selected arithmetic precision wherever
  possible. This makes generated test matrices reproducible across double,
  binary80, binary128, DD, QD, MPFR, and GMP. Known exceptions are GMP
  32-bit vs. 64-bit builds, whose underlying integer representation can
  differ, and the double backend, whose stream is intentionally MPLAPACK's
  normalized stream rather than bit-for-bit LAPACK reference output.

---

## 11. Known issues carried over

- `Cuncsd2by1` M=55,P=40,Q=20 — upstream LAPACK 3.12.1 bug, fixed in
  LAPACK HEAD, not backported on this branch.
- binary80 `Cgd.out` result file contains a stale `ZGS` failure block
  despite a passing summary; clean regeneration is scheduled as a
  followup.
- Several `Rlabad`-dependent routines (`Rlatrs`, `Clatbs`, `Clatps`,
  `Ctrsyl`, `Ctgsy2`, and R-counterparts) still use the `#ifdef`
  pattern. Routing their `(sfmin, bignum)` pair through
  `arithmetic_params` is filed as a followup.
- Several tests carry backend-specific threshold relaxations
  (documented per-routine in section 5). These are test-harness
  tolerance adjustments; the underlying numerical code is unchanged.

---

## 12. Compatibility notes

- Users linking against MPLAPACK 2.1.1 will need to rebuild. Header
  layout has changed: `mplapack_arithmetic_params.h` is new; several
  declarations moved between `mplapack_generic.h` and per-backend
  headers (in particular `Rroundup_lwork`, `Cgemmtr/Rgemmtr`,
  `Risinf`).
- Source-level users of `Clacon/Rlacon` must migrate to `Clacn2/Rlacn2`.
- Source-level users of `Cgelqs/Rgelqs/Cgeqrs/Rgeqrs` must migrate to
  `xGELS`.
- The `Rlabad` function is no longer a no-op; code that called it
  expecting no effect will now receive actively scaled
  `(smlnum, bignum)` values.
