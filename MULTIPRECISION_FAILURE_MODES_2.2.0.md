# Multi-Precision Backend Failure Patterns Exposed During LAPACK 3.12.1 Migration

Companion note to `CHANGES.2.2.0.md`. The MPLAPACK 2.1.1 → 2.2.0 re-base onto
LAPACK 3.12.1 surfaced a recurring set of failure patterns that are not
visible at IEEE binary64 precision. Each pattern is backend-specific,
structural (not a one-off bug), and likely to recur in any future LAPACK sync.

Classification: **A** = GMP arithmetic gaps, **B** = wide-exponent backends,
**C** = DD/QD format semantics, **D** = high-precision stress on test-harness
assumptions.

Commit short-hashes link to
[github.com/nakatamaho/mplapack](https://github.com/nakatamaho/mplapack).

---

## ⚠ Warning — Multi-precision result reproducibility is platform-dependent

**An identical integer seed and an identical pseudo-random integer
stream do not guarantee bit-identical numerical results across
different (compiler, architecture, operating system, libm)
combinations.** This applies not only to MPLAPACK's own test harnesses
but to *any* downstream computation built on top of MPLAPACK that
exercises a high-precision backend. Expect last-digit differences
between platforms as the default, not the exception.

**What is deterministic across platforms:**

- The integer seed `iseed[4]` passed to `Rlaruv`.
- The raw integer output stream of `Rlaruv` itself — a pure linear
  congruential generator over 32-bit integer operations.
- The `REAL`-valued output of `Rlarnd` is also expected to match across
  platforms for a fixed backend and precision. MPLAPACK checks this
  explicitly by comparing the generated stream against
  `mplapack/test/compare/*/Rlaruv_reference.txt`. The known exception is
  GMP/MPFR when crossing ABI families such as LP64 and LLP64; their
  exponent and limb-related widths can differ even when the nominal
  precision setting is the same.
- MPFR transcendental functions are deterministic for a fixed MPFR
  precision/exponent configuration, MPFR library version, and ABI; MPFR
  should not be grouped with the GMP `double` fallback case.

**What is *not* deterministic across platforms:**

- Any accumulated quantity that passes through FMA contraction,
  instruction scheduling, or inlining choices the compiler is free to
  make.
- DD / QD intermediate results, whose renormalization order depends on
  the evaluation order of their internal `double` flops (not fully
  controlled by `-ffp-contract=off`).
- binary80 computations affected by `FLT_EVAL_METHOD`, FPU control-word
  state, or internal-register extended-precision rounding.
- binary128 computations on toolchains that disagree on
  `__float128` / `_Float128` ABI, libquadmath vs. compiler-builtin
  transcendentals.
- GMP transcendentals, which fall back through platform-specific
  `double` libm (see pattern **A1**).

**What this has been observed to affect so far:**

- **Test harness**: `Cget23` test(7) `VL`-consistency, failing
  divergently on x86-64 + icx, ARM64 + Apple clang (M4), and x86-64 +
  mingw64 gcc. A targeted tolerance fix covers the specific harness
  (see D6); the underlying cause is not addressed.
- **Test residuals and histograms** more generally: last-digit
  divergence in reported ratios across the same three environments plus
  others. Stays within documented thresholds but defeats bit-identical
  cross-platform comparison. Under investigation (see O1).

**Practical guidance:**

1. If you need bit-identical numerical results for backends that route
   through native floating-point or libm behavior (DD/QD, binary80,
   binary128, and GMP fallback paths), pin the complete toolchain
   (compiler + version + flags + libm + OS + architecture). Matching
   only the seed and algorithm is not sufficient for those paths. MPFR
   is different: for a fixed MPFR library version, precision, exponent
   range, ABI, and code path, MPFR arithmetic and transcendentals are
   expected to be deterministic.
2. Last-digit differences in test ratios between your machine and a
   reference run can be expected behavior for non-`double` backends that
   use native floating-point/libm paths. MPFR mismatches should be
   treated more narrowly: first verify the `Rlarnd` reference stream and
   then look for a non-MPFR operation or a genuine regression.
3. When filing a suspected numerical regression, always include the
   full (compiler, version, flags, libm, OS, architecture) tuple and
   the full backend (GMP precision bits, MPFR exponent range, etc.);
   without these, the failure may simply not reproduce.
4. For DD/QD, pay special attention to FMA contraction in the internal
   `double` operations. If the compiler rewrites the error-free
   transforms into fused operations, the compensation terms can be
   destroyed and the effective result may collapse toward ordinary
   double precision (~16 decimal digits). `-ffp-contract=off` is a
   necessary guard, but not a complete proof that every relevant internal
   operation was preserved as intended.

---

## Category A — GMP arithmetic gaps

GMP `mpf_class` provides the four basic operations and `sqrt` natively, but
lacks transcendentals and several C++ utility overloads. The shortfall is
papered over by implicit conversion, which silently collapses results to
double.

| Pattern | Mechanism | Diagnostic fingerprint | Mitigation applied in 2.2.0 | Commits |
|---|---|---|---|---|
| **A1. Transcendental functions fall back to double** | `pow(mpf_class, mpf_class)`, `sin`, `cos`, `log`, `exp` have no native `mpf_t` implementation. MPLAPACK's utils route them through `get_d()`, producing a 53-bit intermediate that is then promoted back to `mpf_class`. | Routine visibly uses high precision at input and output, but an intermediate (e.g. `pow(eps, -1/8)` in BDSQR tolmul) carries 53-bit noise. Error ratio capped at ~2⁻⁵³ relative (e.g. ~1e+136 when `eps ~ 1e-155`). | Replace transcendentals with exact or algebraic reformulations: `pow(eps, -1/8)` → explicit root; `sin/cos` in Givens generation → rational double-angle unit-circle sampler enforcing `c² + s² = 1`. Guarded under `___MPLAPACK_BUILD_WITH_GMP___`. | [`fbb8a4f`](https://github.com/nakatamaho/mplapack/commit/fbb8a4f6c5577ba101a72d4fa972a5acf29456af) Rbdsqr pow fix · [`a9510a3`](https://github.com/nakatamaho/mplapack/commit/a9510a32e71a019177560eff2667eafda0180295) extend to Cbdsqr/Rbdsvdx/Rbbcsd/Cbbcsd · [`9e8abb4`](https://github.com/nakatamaho/mplapack/commit/9e8abb44cb7c8a1b413c811679a9de69eb3205a2) Rlatms c²+s²=1 · [`c671869`](https://github.com/nakatamaho/mplapack/commit/c671869376b5ca92308a768705da396aa2071172) Rlatms rational rotation · [`90fd854`](https://github.com/nakatamaho/mplapack/commit/90fd854f8b91bb0da3b4ceb4a51498f6d3e1dfeb) extend to Rlatmt/Clatms/Clatmt |
| **A2. Exact-zero division raises SIGFPE** | Unlike IEEE 754 which returns `±Inf` or `NaN`, GMP's `mpf_div` aborts the process on exact zero divisor. | `SIGFPE`/`abort` at call sites where upstream LAPACK relies on Inf propagation (e.g. `big/sqrt(aaqq)` after `Classq` returns `aaqq = 0` for a zero column). | Explicit `if (aaqq == 0)` guard before `sqrt(aaqq)` / division, per-branch. Applied in `Rgesvj` (3 branches: upper/lower/general), `Rgejsv`, `Cgesvj` (3 branches), `Cgejsv`. Guard is GMP-specific; IEEE behavior is preserved for other backends. | [`433d057`](https://github.com/nakatamaho/mplapack/commit/433d0577a38e) Rgesvj/Rgejsv/Cgesvj/Cgejsv guard |
| **A3. Missing utility overloads cascade to `std::` via `operator double()`** | `mpf_class` defines an implicit `operator double()`. When a utils header lacks an `mpf_class` overload for `abs/floor/ceil/trunc/fmod/copysign/...`, ADL resolves to `std::` through the implicit conversion, silently truncating to 53 bits. | Same ~2⁻⁵³ relative fingerprint as A1. Distinguishing signal: symptom reproduces on routines that do not call any transcendental — the leak is in a utility call. | Diff GMP utils header against MPFR utils header first during debugging; add the missing overload. Defensive addition of `abs(mpf_class)` and `sqrt(mpf_class)` overloads backed by `mpf_abs` / `mpf_sqrt` directly. | [`fe2839a`](https://github.com/nakatamaho/mplapack/commit/fe2839a76586b1eec41c9066686c7dcae8d48af6) add abs/sqrt mpf_class overloads |
| **A4. IEEE feature probes are not portable to GMP/DD/QD arithmetic** | `iMieeeck` deliberately forms infinities, NaNs, and signed zeros by evaluating expressions such as `one / zero`. GMP aborts on exact zero division, and DD/QD do not implement a complete IEEE 754 exception/signaling/signed-zero model. Running the probe therefore tests backend semantics rather than hardware IEEE behavior. | GMP can abort during the probe; DD/QD can produce false IEEE capability conclusions. Downstream LAPACK code may then select IEEE-dependent paths that the backend does not actually support. | Return `0` immediately from `iMieeeck` for GMP, DD, and QD. This advertises "do not rely on IEEE exceptional arithmetic" and avoids executing the unsafe probe. | [`96802ee`](https://github.com/nakatamaho/mplapack/commit/96802ee691c05a00baf8cf1103a6e3f88949d946) initial FABLE patch set |
| **A5. Some CSD bidiagonal-block reductions remain unsupported for GMP** | `Cunbdb` and `Rorbdb` require `atan2` when extracting CSD angles. GMP `mpf_class` has no native `atan2`, and falling back through `double` would destroy the intended precision. Allowing these routines to run would therefore give misleading failures or backend aborts rather than a valid GMP result. | GMP CSD-related tests can terminate inside the reduction or silently pass through a double-precision angle path. | Fail explicitly at entry for GMP with a clear MPLAPACK error message. This is a known limitation, not a numerical tolerance adjustment. | [`96802ee`](https://github.com/nakatamaho/mplapack/commit/96802ee691c05a00baf8cf1103a6e3f88949d946) initial FABLE patch set |

---

## Category B — Wide-exponent backends (binary80, binary128, MPFR, GMP)

LAPACK 3.12.1 removed `Rlabad` under the assumption that IEEE binary64's
`smlnum ~ 2.2e-308` and `bignum ~ 1.8e308` are moderate enough that
`smlnum/eps`, `sqrt(smlnum)/eps`, etc. stay representable. This assumption
breaks at binary80 (`smlnum ~ 1e-4932`), binary128 (`smlnum ~ 1e-4898`), and
arbitrary MPFR/GMP configurations.

| Pattern | Mechanism | Diagnostic fingerprint | Mitigation applied in 2.2.0 | Commits |
|---|---|---|---|---|
| **B1. `Rlabad` removal is not safe outside IEEE binary64** | `Rlabad` clamps `smlnum`/`bignum` toward `sqrt(smlnum)`/`sqrt(bignum)`. Upstream's removal is correct for IEEE; for wide-exponent backends, `smlnum/eps` can evaluate to a value far outside the safe range, poisoning downstream scaling. | `Cdrvev` type-1 zero matrix and types 19/20/21 random matrices exhibit ratios ~1e+136. `Cgges` result-11 (alpha/beta vs. Schur diagonal) drifts on scaling-sensitive cases. `Cgeev` ZEV test-7 sentinel (left eigenvectors from combined vs. left-only path) diverges. | Reinstate `Rlabad` as an active hand-written reference routine (not a no-op stub). Re-inject `Rlabad(smlnum, bignum)` before scaling at: `Cgeev`, `Cgges`, `Cgges3`, `Clatrs`, `Rget38`, `Cqrt13`, `xLATB4`. Guard per-backend as needed. Long-term plan: route `(sfmin, bignum)` through `arithmetic_params` to eliminate the `#ifdef` pattern. | [`e6823f8`](https://github.com/nakatamaho/mplapack/commit/e6823f8f0d44) revert upstream Rlabad removal · [`21687206`](https://github.com/nakatamaho/mplapack/commit/21687206b8d32923140bb1e031c6012b450f2e85) reinstate as active routine · [`af29559`](https://github.com/nakatamaho/mplapack/commit/af29559cdfbb36611ad1961243d7eec9b551728e) follow-up · [`d3e12a8`](https://github.com/nakatamaho/mplapack/commit/d3e12a8e76a5d72a0d027c7ffa4d665e3ff390a5) refactor wide-exp range reduction · [`936e97b`](https://github.com/nakatamaho/mplapack/commit/936e97bff2d59966e39bbef0b11dd61d6fb3a714) Clatrs · [`9fce663`](https://github.com/nakatamaho/mplapack/commit/9fce663736bae1046fba7b15bfa510a87bd6aefd) Cgeev/Cgges · [`4c05b55`](https://github.com/nakatamaho/mplapack/commit/4c05b55468f4161ce99f826f51b26f11b6f05bc4) Cgges3 · [`f465c3d`](https://github.com/nakatamaho/mplapack/commit/f465c3d05d96ee6c2e3c736026375c964c93ff6a) Rget38 · [`394ca36`](https://github.com/nakatamaho/mplapack/commit/394ca3680552) Cqrt13/Cqrt15 · [`004992b`](https://github.com/nakatamaho/mplapack/commit/004992b6408b) xLATB4 scope narrowing · [`710a134`](https://github.com/nakatamaho/mplapack/commit/710a13400769) xLATB4 restore |
| **B2. `xGEJSV/xGESVJ` underflow guards need a reciprocal safe bound** | `Rgejsv/Cgejsv/Rgesvj/Cgesvj` use `sfmin`, `sqrt(sfmin)`, `big`, and `sqrt(big)` together when deciding whether scaled column norms can be manipulated without underflow. For DD/QD/MPFR/GMP, `Rlamch("Overflow")` may not be the reciprocal partner of the `sfmin` used by the routine, especially after wide-range safe-minimum adjustments. Mixing those values can make the lower-end guard region inconsistent. | SVD/Jacobi-SVD paths become sensitive to extreme small scaling: columns near `sfmin` can be classified against mismatched lower and upper bounds, causing avoidable underflow-guard failures or backend-specific test drift. | For DD, QD, MPFR, and GMP builds, set `big = one / sfmin` in `Rgejsv/Cgejsv/Rgesvj/Cgesvj`; keep upstream `big = Rlamch("Overflow")` for binary IEEE-style backends. Mirror the guarded change in the FABLE patches (`patch-Cgejsv.cpp`, `patch-Cgesvj.cpp`, `patch-Rgejsv.cpp`, `patch-Rgesvj.cpp`). Long-term plan: route this through the shared arithmetic-params `(sfmin,bignum)` pair. | [`1e2d306`](https://github.com/nakatamaho/mplapack/commit/1e2d3068d) regenerate LAPACK 3.12.1 reference patches |
| **B3. MPFR binary64/binary128 emulation underflows faithfully** | MPFR with restricted exponent range reproduces IEEE denormal behavior exactly. Squared-norm intermediates `\|z\|²` can underflow even when the individual inputs `\|re\|, \|im\|` are representable. The scaled-vs-unscaled dispatch based on componentwise maxima is not a sufficient gate. | `Clartg` Z-only unitarity loss (~1e-04) with Q clean, firing only on high `jtype` with mixed large/small scaling. Asymmetric: unitarity of Q stays intact. | Tighten the unscaled-path gate from componentwise max to a ratio test: `fsmall/f1 > sqrt(eps)`. Split `sqrt(f2*h2)` into `sqrt(f2) * sqrt(h2)` so that neither product under/overflows. Applied in `Clartg`. | [`e6536f6`](https://github.com/nakatamaho/mplapack/commit/e6536f65adfa98cd2aaca1a27a07ee1394272723) Clartg subnormal hazard |
| **B4. MPFR default exponent range is effectively unbounded** | If the user does not configure MPFR's exponent range, the defaults are large enough that `Rlamch` reports `rmin/rmax` as values LAPACK's algorithms never anticipated. Rstebz bisection, Sturm counts, and Blue-scaling constants become pathological. | Rstebz enters bisection loops; Blue-scaling thresholds fall outside the usable range; subnormal-region classification is ambiguous. | Truncate MPFR exponent range at initialization. Emit a runtime warning from `make_blue_scaling_params()` if BlueScale thresholds are still outside range. | [`780f4c9`](https://github.com/nakatamaho/mplapack/commit/780f4c906e59) mpfr fallback lower exponent · [`1fbf252`](https://github.com/nakatamaho/mplapack/commit/1fbf252eebfe) Blue-scaling guard warnings |
| **B5. Iteration limits were calibrated for double precision** | `xBDSQR` (`maxitr = 6`), `xHGEQZ` (`maxit = 30`), `xBBCSD` inherit their limits from LAPACK's double-precision tuning. High-precision backends use tighter O(eps) stopping thresholds and legitimately need more sweeps to reduce off-diagonal entries below tolerance. | Non-convergence reported (`INFO != 0`) on matrices that converge cleanly at double. No numerical wrongness — simply premature abort. | Scale iteration limit as `iceil(base × log(eps) / log(double_eps))` in REAL arithmetic, floored at `base`, capped at `base × 64` (or 1000 for HGEQZ). Apply symmetrically to `Rbdsqr/Cbdsqr`, `Rhgeqz/Chgeqz`, `Rbbcsd/Cbbcsd`. Mirror in FABLE patches. | [`0082767`](https://github.com/nakatamaho/mplapack/commit/0082767f8be355861b96bfe83a4106566e3489f0) BDSQR scale for xGEDMD · [`5602c2f`](https://github.com/nakatamaho/mplapack/commit/5602c2fb6a8a5d21796bf8ccd6f04de76c1d083a) HGEQZ scale |
| **B6. Eigenvalue-interval refinement counts scale with exponent span** | `Rstebz` and the `Rlarr*` refinement routines compute bisection/refinement iteration counts from `log((range + pivmin) / pivmin)`. With MPFR/GMP safe minima and wide exponent spans, this can become O(1e8) even for a test that is not intended to spend that much time in interval refinement. | SEP/STEBZ-style tests appear to hang or run for impractical time in MPFR/GMP builds; the failure is dominated by exponent-span bookkeeping rather than eigenvalue accuracy. | Cap `itmax`/`maxitr` at 100000 for MPFR/GMP in `Rstebz`, `Rlarrb`, `Rlarrd`, and `Rlarrk`. This keeps the test finite while preserving a much larger budget than binary64 ever uses. | [`1e2d306`](https://github.com/nakatamaho/mplapack/commit/1e2d3068dbaf66d3e51a48f7977bd6b7085d2af6) regenerate LAPACK 3.12.1 reference patches |

---

## Category C — DD / QD format semantics

DD (double-double) and QD (quad-double) represent numbers as an unevaluated
sum of 2 or 4 IEEE doubles. This gives extended precision but not an IEEE
format, and the upstream library's conventions leak into LAPACK's
assumptions.

| Pattern | Mechanism | Diagnostic fingerprint | Mitigation applied in 2.2.0 | Commits |
|---|---|---|---|---|
| **C1. `operator==` is componentwise, not structural** | `qd_real::operator==` compares the four internal `double` components directly. Two `qd_real` values whose mathematical difference renormalizes to exactly zero can still disagree in the last component pair. | `Rlaexc` canonical real Schur form check (`Rget34`) reports spurious inequality for entries that are numerically identical. | Replace the raw `==` comparison with a structural-equivalence check in `Rget34` for QD only. Preserve exact comparison for all other backends. Scope the change to this single test; do not touch `operator==` globally. | [`f882e14`](https://github.com/nakatamaho/mplapack/commit/f882e14b69e14a0f7445ae3a7548b4b0e1743df3) QD structural comparison in Rget34 |
| **C2. Upstream-defined `_eps` and `_min_normalized` inconsistent with LAPACK's safe-range pair** | The QD library exposes `dd_real::_eps`, `qd_real::_eps`, `_min_normalized`, etc. Using these directly inside `arithmetic_params` conflates the format's advertised precision with the `(sfmin, bignum)` pair LAPACK needs. The format's precision is also not always a clean `2^(-n)`. | Blue-scaling thresholds land on the wrong power-of-two; test tolerances drift backend to backend in a way that does not track `eps`. | Derive `emin` for `dd_real`/`qd_real` from `frexp(_min_normalized)` rather than from the upstream constants. Stop deriving `eps`/precision from `dd_real::_eps` / `qd_real::_eps`. Treat `arithmetic_params::safmin`/`safmax` as MPLAPACK-internal power-of-two bounds, documented as such. | [`231b5ac`](https://github.com/nakatamaho/mplapack/commit/231b5ac4f76d) emin from frexp · [`6ec9ca4`](https://github.com/nakatamaho/mplapack/commit/6ec9ca435226) drop DD_EMIN constant · [`d887678`](https://github.com/nakatamaho/mplapack/commit/d887678102fc) pin canonical DD/QD E and P · [`1fa244c`](https://github.com/nakatamaho/mplapack/commit/1fa244ca3e2a) align dd/qd params with shared layout |
| **C3. QD residual harnesses can amplify unevaluated-product artifacts** | In QD, forming tiny products first and only scaling the completed residual later can preserve component-level cancellation artifacts that are irrelevant to the mathematical residual. `Rget32`'s Rlasy2 check divides that residual by a denominator involving `tnrm * eps * xnrm`, so a harmless artifact becomes a huge ratio. | `Rget32/Rlasy2` fails at QD with enormous ratios while the solution entries are the expected small integer-like values and the error is localized to cancellation in the test residual construction. | For QD only, scale the coefficient side and solution/right-hand-side side before forming the four residual rows. Keep the original residual expression for all other backends. | [`93ecec0`](https://github.com/nakatamaho/mplapack/commit/93ecec0b1cb149c63834ff2bf91096a5ab3b8286) QD Rlasy2 residual check |

---

## Category D — High-precision stress exposes test-harness assumptions

LAPACK's test harnesses were written for IEEE binary64. Many of their
implicit assumptions — about generator magnitudes, scale conventions, and
comparison granularity — stop holding once `eps` is many orders of magnitude
smaller and the exponent range is wider.

| Pattern | Mechanism | Diagnostic fingerprint | Mitigation applied in 2.2.0 | Commits |
|---|---|---|---|---|
| **D1. Random-matrix generators explode at high precision** | Generators such as `bd[j] = exp(temp1 × Rlarnd(2))` with `temp1 = -2 × log(ulp)` were sized for double (`temp1 ≈ 73`, `exp(temp1) ≈ 8e31`). Applying the same rule gives `temp1 ≈ 144` for DD (`exp ≈ 4e62`) and `temp1 ≈ 289` for QD (`exp ≈ 1e125`). The test then measures underflow-boundary and Sturm-count edge behavior, not the algorithm under test. | SVD type-16 tests fail with ratios dominated by backend-specific scaling behavior rather than by the SVD itself. Failures concentrate near boundary types; the complex QD path can also trip the `Rsvdch` Sturm precheck before the SVD residuals are the limiting issue. | Clamp `temp1` per backend: `-0.5 × log(ulp)` for DD, `-0.25 × log(ulp)` for QD — both give `temp1 ≈ 36`, `exp ≈ 4.5e15`, matching double's `1/eps`. Keep IEEE backends unchanged. Skip the `Rsvdch` precheck in the QD complex bidiagonal harness so the test measures the SVD path rather than the auxiliary Sturm checker. | [`43cc5fc`](https://github.com/nakatamaho/mplapack/commit/43cc5fc74f50d77e013fb45594da0de08f4348a7) DD/QD SVD harness edge cases |
| **D2. Bit-exact comparison promotes absolute-scale noise to `ulpinv`** | Test harnesses using `if (a == b) pass else fail = ulpinv` treat any deviation as catastrophic. At double precision, deviations between two mathematically-equal computations are either exactly zero or visibly large. At high precision, O(eps × \|value\|) absolute noise on correctly-computed quantities fails bit-exact equality while remaining numerically identical. Root cause: **precision**. | `Rdrgev3` result 5/6/7 (W consistency) reports `ulpinv` on MPFR/GMP while the numerical answer is identical to within a tiny fraction of the eigenvalue magnitude. | Use a scaled tolerance for MPFR and GMP: compare relative to the eigenvalue scale. Keep bit-exact comparison for other backends. Patch mirrored in FABLE test-patch set. | [`b617627`](https://github.com/nakatamaho/mplapack/commit/b61762713e373179d1caeb856f1af2dc7c53cfcf) MPFR Rdrgev3 W-compare tolerance · [`d581a013`](https://github.com/nakatamaho/mplapack/commit/d581a0136cae9f6c4b9d2a98d26a392272901c95) extend tolerance to GMP |
| **D3. Jacobi SVD's separate scale-factor convention unhandled in tests** | `*gesvj` returns singular values paired with a scale factor in `WORK(1)` (real) or `RWORK(1)` (complex). The usual SVD driver tests, written for the standard convention, reconstruct `A = U Σ V^T` without applying the scale and declare a large residual. | At double, the scale is typically `1.0` and the bug is latent. At high precision, generators produce scale factors that are non-trivial, and the test fails. | Apply `WORK(1)` / `RWORK(1)` before the reconstruction residual check in both real and complex SVD driver tests. | [`43cc5fc`](https://github.com/nakatamaho/mplapack/commit/43cc5fc74f50d77e013fb45594da0de08f4348a7) part of DD/QD SVD harness commit |
| **D4. Upstream fatal-exit codes hide recoverable outcomes** | `Rgges`/`Rgges3` return `INFO = N+3` when eigenvalue reordering is imperfect but the Schur form is still valid. The harness treats `INFO ≠ 0` as fatal and emits `ulpinv` before running residual, Schur-form, or SDIM checks. | DD DGS type-17, GMP DGS cases report `ulpinv` while all downstream checks would pass. | Treat `INFO = N+3` like `INFO = N+2` in `Rdrges`/`Rdrges3`: non-fatal, continue to downstream checks. Mirror in patched test sources. | [`13fe16a`](https://github.com/nakatamaho/mplapack/commit/13fe16a34376d8a1ad10ce00dfa5edd4947eabba) Rdrges3 for DD · [`e366886`](https://github.com/nakatamaho/mplapack/commit/e3668865a3129ada7c9712b514b01a38b7758184) Rdrges extension |
| **D5. `CNDNUM = 0.1/EPS` produces meaningless stress cases at high precision** | GB driver type-6 tests construct band matrices with condition number `0.1/EPS`. At double, this is ~`1e15`, a reasonable stress. At DD, `~1e32`; at binary128 / MPFR, far beyond what iterative refinement can usefully distinguish from random noise, and the expected backward error is dominated by the artificial condition number rather than by the algorithm. | `Cgbsvx`/`Rgbsvx` backward-error checks fail on type-6 while all other GB types pass. | Cap type-6 `CNDNUM` in `Rlatb4/Clatb4`: `1e24` for DD, binary128, MPFR; `1e30` for QD, GMP. Leaves the stress test meaningful without driving the harness into noise. | [`3781953`](https://github.com/nakatamaho/mplapack/commit/378195313a8eb88abe6797c8969ccf446b73947b) cap type-6 CNDNUM |
| **D6. binary80/binary128 bit-exact checks can fail on numerically-identical eigenvectors** | `Cget23` compares left eigenvectors from two nominally equivalent driver paths with bit-exact equality. On binary80 and binary128, the two paths can differ by a few ulps while all residual, orthogonality, and normalization checks remain clean. The precise root cause is not yet known. It has been observed across multiple toolchains and architectures, so this should not be documented as an FMA-only or single-compiler effect. | `Cget23` test(7) `VL`-consistency fails on binary80 / binary128 with icx on x86-64 and on binary128 with Apple clang on Apple M4; a similar-looking `VL`-consistency failure was also seen on mingw64. test(1)–(6) pass cleanly throughout. The mismatching entries are numerically identical within a small scaled tolerance, but the bit-exact sentinel reports `ulpinv`. | For binary80 and binary128 only, compare each left-eigenvector column by max-norm against a 100-ulp tolerance in `Cget23`. Keep bit-exact comparison for all other backends. Matching FABLE test patch added. **The underlying reason for the bit differences is not resolved in 2.2.0**; only this over-strict harness check is adjusted. | [`27ce68e`](https://github.com/nakatamaho/mplapack/commit/27ce68e267b92bee50c50775ca5ffb068657034e) binary80/binary128 Cget23 VL tolerance |

---

## Related infrastructure (Category B, arithmetic_params layer)

The `mplapack_arithmetic_params.h` layer introduced during 2.2.0 is the
intended long-term home for the per-backend constants that B1–B6
currently handle ad hoc. The refactor series:

- [`614456b`](https://github.com/nakatamaho/mplapack/commit/614456b38c7d)
  introduce integer Rlamch metadata and unified Blue scaling builder (WIP)
- [`808b32f`](https://github.com/nakatamaho/mplapack/commit/808b32f64711)
  migrate double and binary128 params
- [`deebc40`](https://github.com/nakatamaho/mplapack/commit/deebc40472a7)
  migrate binary80 params
- [`1fa244c`](https://github.com/nakatamaho/mplapack/commit/1fa244ca3e2a)
  align dd/qd params
- [`f3f7178`](https://github.com/nakatamaho/mplapack/commit/f3f7178d5f17)
  align MPFR params
- [`515f2de`](https://github.com/nakatamaho/mplapack/commit/515f2deb09f2)
  align GMP params
- [`e0069b6`](https://github.com/nakatamaho/mplapack/commit/e0069b661315)
  move Rlamch backend logic to arithmetic params layer
- [`0a153fc`](https://github.com/nakatamaho/mplapack/commit/0a153fcdafe9)
  revert compute_safmin exponent selection to Fortran reference (`max`, not `min`)
- [`640ea36`](https://github.com/nakatamaho/mplapack/commit/640ea3635f17)
  clarify internal safmin/safmax helpers
- [`575ffc5`](https://github.com/nakatamaho/mplapack/commit/575ffc5f4f24)
  make safmax finite; relax MPFR Blue checks in stress environments
- [`aa731f0`](https://github.com/nakatamaho/mplapack/commit/aa731f0fbb23)
  strengthen binary80 Rlamch / BlueScale tests
- [`b005349`](https://github.com/nakatamaho/mplapack/commit/b005349a5f9c)
  strengthen binary128 Rlamch / BlueScale tests
- [`62e8b63`](https://github.com/nakatamaho/mplapack/commit/62e8b63c8a71)
  GMP arithmetic-params consistency checks and Blue invariants
- [`c388114`](https://github.com/nakatamaho/mplapack/commit/c388114666b3)
  BlueScale near-threshold classification checks
- [`0dfb8d3`](https://github.com/nakatamaho/mplapack/commit/0dfb8d366b08)
  strict boundary tests for Blue scaling constants
- [`85211ce`](https://github.com/nakatamaho/mplapack/commit/85211cea0dd5)
  fix Rlassq Blue scaling via arithmetic_params
- [`6a710ea`](https://github.com/nakatamaho/mplapack/commit/6a710ea91a70)
  fix Classq Blue scaling; add abssq helper to Clartg
- [`b79005a`](https://github.com/nakatamaho/mplapack/commit/b79005a4588d)
  add arithmetic parameters headers

---

## Open investigations

Failure modes that are observed but not yet resolved in 2.2.0; documented
here to prevent repeated rediscovery and to mark them as known unknowns
rather than silent anomalies.

### O1. Cross-platform non-determinism of LAPACK RNG-seeded test output

See the prominent warning at the top of this document for the
user-facing statement. This entry is the engineering-side investigation
record.

**Observation.** Despite a bit-identical integer seed and a
bit-identical raw integer stream from `Rlaruv` (a pure linear
congruential generator), downstream test residuals and histograms
differ in the last several digits across (compiler, OS, CPU, libm)
combinations with identical backend and identical MPLAPACK revision.
The divergence is real, reproducible, and not reducible to a single
toolchain or platform.

**Hypotheses under consideration**, in rough order of suspicion:

- The first divergence point is not expected to be `Rlarnd` itself:
  both the integer stream and the real-valued `Rlarnd` outputs are
  checked against reference files and should match across platforms for
  a fixed backend, precision, and ABI. The known exception is crossing
  ABI families such as LP64 and LLP64 for GMP/MPFR, where exponent and
  limb-related widths can differ and therefore change generated
  high-precision values.
- For DD/QD, `dd_real`/`qd_real` renormalization depends on the
  evaluation order of internal `double` flops. This is distinct from
  FMA contraction and is not fully suppressed by `-ffp-contract=off`.
- For binary80, `FLT_EVAL_METHOD` and FPU control-word state (80-bit
  vs. 64-bit internal rounding) differ across `gcc`, `icx`, `clang`,
  and even across build modes of the same compiler.
- For binary128, `__float128` vs. `_Float128` ABI handling differs
  between toolchains; transcendentals go through `libquadmath` on some
  and through compiler builtins on others.
- For GMP, the transcendental fallback path (pattern **A1**) routes
  through platform-specific `double` libm, so any code that calls
  `pow`, `sin`, `cos`, `log`, `exp` on a high-precision value is
  effectively running libm-specific code at that step. MPFR does not
  share this failure mode; its transcendentals are computed by MPFR.

**Status.** Under active investigation. No mitigation in 2.2.0.
Observed divergences still fall within the documented test thresholds
on all supported platforms, so pass/fail outcomes remain stable — but
bit-identical cross-platform reproduction is not currently achievable
for any backend other than `double`, and cannot be assumed even there
in all circumstances.

**Relationship to D6.** D6 documents the same class of phenomenon at
harness granularity, confirmed in three independent environments
(`icx` on x86-64, Apple clang on Apple M4, and — pending file —
mingw64 gcc on x86-64 Windows). That cross-toolchain + cross-architecture
+ cross-OS evidence upgrades the link to O1 from "plausible conjecture"
to **strong circumstantial support**: if compiler/libm/ABI-driven
reordering visibly splits one harness result across three unrelated
environments, the same mechanism is the most economical explanation for
the broader divergence in O1. Under this reading, O1 and D6 are two
manifestations of one root cause. D6's targeted fix (a tolerance in
one harness) does not generalize; the proper fix for O1, if the link
holds, is to reduce compiler-visible evaluation-order freedom in the
pipeline (strictly-ordered transcendentals, pinned libm, explicit
renormalization points in DD/QD) rather than to accumulate more
tolerances.

**Proximate investigation goal.** Pin down the first divergence point
on a controlled two-platform comparison: run the same seed through
`Rlaruv` (expect bit-identical integer stream), then through `Rlarnd`
(expect divergence to start), and bisect by replacing individual libm
calls with hand-coded alternatives until the integer-stream-to-REAL
translation is deterministic. Without this, any "fix" is only validated
on the machine it was developed on.

### O2. DD `Claror` is optimization-sensitive with GCC + libqd

`Claror` carries a DD-only `#pragma GCC optimize("O0")` workaround in
the FABLE patch set. The immediate symptom was localized to the
Householder-vector update around `x[kbeg-1] += xnorms` under GCC +
libqd, but the exact optimizer transformation has not yet been
isolated. The patch therefore documents an observed compiler/backend
interaction rather than a proven arithmetic algorithm defect.

**Status.** Workaround only. Keep this separate from the threshold and
test-harness entries above: if it reappears, investigate compiler
optimization of DD renormalization and complex Householder updates
before changing numerical tolerances.

---

## Recurring diagnostic heuristics

These shortcuts paid off repeatedly during 2.2.0 debugging and are worth
keeping in mind for any future LAPACK sync:

1. **Relative error ratio capped at ~2⁻⁵³** (e.g. ratio ~`1e+136` when
   `eps ~ 1e-155`) is a decisive fingerprint of double-precision fallback
   inside a multi-precision routine. Look for A1/A3, in that order.
2. **Asymmetric unitarity loss** — one of `Q`/`Z` clean, the other degraded
   — under mixed large/small scaling points at B3 (subnormal hazard), not
   at the matrix algorithm itself.
3. **Failures concentrated at high `jtype`** typically indicate a
   test-harness / generator issue (D1, D5), not a numerical bug in the
   routine under test.
4. **Symptom vanishes at double precision, reappears at DD/QD/binary128**
   without any pattern to the matrix content — suspect B1 (`Rlabad`
   removal) or B5 (iteration limit) before digging into the algorithm.
5. **When in doubt, diff the GMP utils header against the MPFR utils
   header first** (A3). Missing overloads are the single most common cause
   of silent precision loss on GMP.
6. **A bit-exact variant-consistency check fails but residual/normalization
   pass** — suspect D6 (compiler reordering), not an algorithm bug. Swap
   `-ffp-contract=off` or switch compiler to confirm before changing any
   source. If confirmed, also check whether the same computation feeds
   into O1 (cross-platform non-determinism).

---

## What is *not* in this taxonomy

For completeness, failure modes deliberately excluded from the above:

- **Test threshold relaxations** — these are test-harness tolerance
  adjustments, not backend-specific failure patterns. See `CHANGES.2.2.0.md`
  §5 for the list.
- **Upstream LAPACK 3.12.1 bugs** (e.g. `Cuncsd2by1` M=55,P=40,Q=20) —
  these fail identically at double precision and are not multi-precision
  specific.
- **FABLE/cout.py translator bugs** — these affect code generation, not
  arithmetic behavior. They happen to surface during multi-precision
  migration because that is when the code is regenerated, but the defect
  class is different.
- **Build-system issues** (ODR in `Mxlaenv.cpp`, `g++-mp-15` IPA-modref
  miscompile of `xdrgvx`) — unrelated to numeric precision.
