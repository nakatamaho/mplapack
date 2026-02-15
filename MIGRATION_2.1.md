# MIGRATION_2.1.md

This document describes how to migrate **MPLAPACK users and packagers** from older `_Float128/_Float64x`-based naming to the **IEEE 754 naming scheme** used in the `2.1` line:
- `_Float128` → **binary128**
- `_Float64x` → **binary80**

The goal is to remove compiler-/libc-specific type names from user-facing glue code and build targets, while keeping the numerical API stable.

---

## 0. Scope (who should read this)

Read this if you:
- include MPLAPACK headers that previously exposed `_Float128` or `_Float64x`
- have downstream code that used those types directly (including `std::complex<_Float128>` etc.)
- rely on example Makefiles / helper scripts / target names that contained `__Float128` / `__Float64x`
- maintain a package that must keep ABI/SOVERSION and link targets consistent

If you only use `double`, `dd`, `qd`, `mpfr`, or `gmp`, you may not need changes.

---

## 1. Breaking changes (practical)

### 1.1 C++ standard
MPLAPACK 2.1 requires **C++17** (GNU extensions).

**Action**
- CMake: `set(CMAKE_CXX_STANDARD 17)` (or higher)
- Autotools projects: add `-std=gnu++17` (or an equivalent toolchain setting)

### 1.2 Naming standardization
The naming visible in code, example targets, and some build glue was standardized:
- `_Float128/_Float64x` **labels** → `binary128/binary80`
- compiler-specific types → MPLAPACK **portable type aliases**

---

## 2. Quick replacement table (do this first)

### 2.1 In C/C++ source code (types)
| Old (compiler-specific) | New (portable alias) |
|---|---|
| `_Float128` | `mplapack_binary128_t` |
| `_Float64x` | `mplapack_binary80_t` |

Typical patterns:
- `std::complex<_Float128>` → `std::complex<mplapack_binary128_t>`
- `std::complex<_Float64x>` → `std::complex<mplapack_binary80_t>`

### 2.2 In file/dir/target *names* (strings)
Use these when you have filenames, directories, or target names containing `_Float*`:
| Old label | New label |
|---|---|
| `_Float128` | `binary128` |
| `_Float64x` | `binary80` |

**Important:** code types use `mplapack_binary{80,128}_t`, while labels in names use `binary{80,128}`.

---

## 3. Header-level changes you will actually see

### 3.1 Convenience headers now typedef to portable aliases
The umbrella headers (examples: `include/mplapack_eig.h`, `include/mplapack_lin.h`, `include/mplapack_matgen.h`) used to typedef `REAL`/`COMPLEX` using `_Float*` in their binary80/128 sections.
They now typedef using:
- `mplapack_binary80_t`
- `mplapack_binary128_t`

**Action**
- If your downstream code relied on those typedefs being `_Float*`, update your code to accept the alias type.
- If you wrote overloads/templates specialized on `_Float128`/`_Float64x`, update them.

---

## 4. Macro / build-configuration deltas (binary128 I/O selection)

Some code previously selected binary128 I/O via:
- `___MPLAPACK_WANT_LIBQUADMATH___`
- `___MPLAPACK_LONGDOUBLE_IS_BINARY128___`

This moved to an enum-style selector:
- `MPLAPACK_BINARY128_IO == MPLAPACK_BINARY128_IO_QUADMATH_SNPRINTF`
- `MPLAPACK_BINARY128_IO == MPLAPACK_BINARY128_IO_SNPRINTF_LDBL`

**Action**
- Replace direct `#if defined ___MPLAPACK_WANT_LIBQUADMATH___`-style checks with comparisons against `MPLAPACK_BINARY128_IO`.
- If you injected these macros from your build system, migrate to the new selector mechanism.

---

## 5. Backend selection: binary80/binary128 “mode” (advanced / packagers)

The 2.1 series introduces mode-based selection for the binary80/binary128 backends (compile-time):
- `MPLAPACK_BINARY80_MODE` (binary80 backend selection)
- `MPLAPACK_BINARY128_MODE` (binary128 backend selection)

The implementation supports multiple modes (conceptually):
- `_Float64x`-based
- 80-bit `long double`-based
- `_Float128`-based
- binary128 `long double`-based
- libquadmath-based
- explicit disabled mode

**Action**
- If you hard-coded assumptions about what “binary80” or “binary128” means on a platform/compiler,
  audit those assumptions and switch to mode-aware configuration.
- To see the exact mode constants used by the branch, grep the source:
  - `grep -R "MPLAPACK_BINARY80_MODE" -n`
  - `grep -R "MPLAPACK_BINARY128_MODE" -n`

---

## 6. Old → new *target name* mapping (examples + build glue)

This section is mainly for people who:
- build/ship the example programs
- copied the example Makefiles into external projects
- refer to “opt” targets by name

### 6.1 Example program targets (observed)
| Old target name | New target name |
|---|---|
| `Rgemm__Float64x_opt` | `Rgemm_mplapack_binary80_t_opt` |
| `Rgemm__Float128_opt` | `Rgemm_mplapack_binary128_t_opt` |

(Windows `.exe` variants follow the same rename.)

### 6.2 Example Makefile variables (observed)
| Old variable name | New variable name |
|---|---|
| `_Float64x_cxxflags` | `mplapack_binary80_t_cxxflags` |
| `_Float64x_libdepends` | `mplapack_binary80_t_libdepends` |
| `_Float64x_opt_libdepends` | `mplapack_binary80_t_opt_libdepends` |

### 6.3 Example link library names (observed)
| Old library (link flag) | New library (link flag) |
|---|---|
| `-lmpblas__Float128_opt` | `-lmpblas_mplapack_binary128_t_opt` |

**Action**
- If your downstream build references these targets/variables/libs, rename them.
- If you do **not** build examples, you can ignore this section.

---

## 7. Mechanical migration recipe (recommended order)

### Step 1: bump compiler standard
- Ensure C++17 is used.

### Step 2: update type usage in your code
Search for the old types:
```sh
git grep -nE '\b_Float128\b|\b_Float64x\b'
````

Then replace:

* `_Float128` → `mplapack_binary128_t`
* `_Float64x` → `mplapack_binary80_t`

### Step 3: update any `_Float*`-labeled filenames/targets

Search:

```sh
git grep -nE '_Float128|_Float64x|__Float128|__Float64x'
```

Then apply label changes where appropriate:

* `_Float128` → `binary128`
* `_Float64x` → `binary80`

And apply the known example target renames (if relevant):

* `Rgemm__Float128_opt` → `Rgemm_mplapack_binary128_t_opt`
* `Rgemm__Float64x_opt` → `Rgemm_mplapack_binary80_t_opt`

### Step 4: rebuild and let the compiler complain productively

Most remaining issues will be:

* missing includes for the alias types
* stale macro checks for binary128 I/O
* stale link names (especially example-derived glue)

---

## 8. Troubleshooting heuristics (fast)

* **`mplapack_binary128_t` unknown**
  You are missing the header that defines it, or binary128 support is disabled in your build configuration.

* **Link error for `-lmpblas__Float128_opt`**
  Your build is using the pre-2.1 name. Rename to the new link target (or stop linking the “opt” variant if you do not need it).

* **You expected `_Float64x` but got `long double` behavior**
  binary80 is a *format label*; the backend may differ by compiler/platform and is selected via mode logic. Treat `mplapack_binary80_t` as the truth.

---

## 9. Known problems (2.1.0)

### 9.1 MPFR default profile: Csg (Hermitian generalized eigenvalue) test failures

**Symptom**

When running the `Csg` eigenvalue test suite under the MPFR default profile (512-bit mantissa, exponent range ≈ ±323 228 496), two classes of failures are observed:

| NB | Failures | Description |
|-----|----------|-------------|
| 1 | 9 / 11172 | Test ratios of order 10^153 (type 8, 9, 11, 18) |
| 3 | 6 / 11172 | Test ratios of order 10^153 (type 17, 19) |
| 20 | fatal | `Cstein` parameter 6 error → test abort |

The same test suite passes completely under the MPFR `binary128` profile (precision=113, emin=−16381, emax=16384).

**Root cause 1 — bisection iteration cap (`itmax` clamping in `Rstebz`)**

`Rstebz` computes the maximum number of bisection iterations as:

```
itmax = ⌊(log(gu − gl + pivmin) − log(pivmin)) / log(2)⌋ + 2
```

For the MPFR default profile this yields itmax on the order of 10^8 or more, which is clamped to 20 000:

```cpp
#if defined ___MPLAPACK_BUILD_WITH_MPFR___
    if (itmax >= 20000)
        itmax = 20000; // XXX itmax can be too large for MPFR
#endif
```

For test matrices with large eigenvalues (kmagn=2, `anorm = rtovfl × ulp × aninv ≈ 10^{161 614 093}`), 20 000 bisection steps reduce the interval by a factor of 2^20000 ≈ 10^6021, leaving the interval width at approximately 10^{161 608 072} — far above the required tolerance of `ulp × |eigenvalue|`. The bisection therefore does not converge, producing meaningless eigenvalue approximations and the observed enormous test ratios.

Under the `binary128` profile the required `itmax` is approximately 32 750, so 20 000 iterations are more than sufficient for the much smaller exponent range.

**Root cause 2 — missing error propagation from `Rstebz` to `Cstein` in `Cheevx`**

When `Rstebz` fails to converge, it marks unconverged eigenvalues with **negative** block indices (`iblock[j] = −jb`). `Cheevx` does not check `Rstebz`'s return code before calling `Cstein`:

```cpp
// Cheevx.cpp — no info check between these two calls
Rstebz(range, &order, n, vll, vuu, il, iu, abstll, ...rwork..., info);
if (wantz) {
    Cstein(n, ...iwork[indibl]..., info);   // receives negative iblock values
```

`Cstein` validates that `iblock` is monotonically non-decreasing; negative values violate this, triggering `info = −6` (`parameter 6 had an illegal value`). This is the direct cause of the NB=20 fatal error. The same latent issue exists for NB=1 and NB=3 but is masked by the test driver's skip logic.

**Scope and impact**

- Affects **only** the MPFR default profile where the exponent range is extremely large relative to the mantissa precision.
- Does **not** affect the `binary128` or `binary80` profiles, nor `dd`, `qd`, or `gmp` backends.
- Does **not** affect user code that calls `Cheev`/`Cheevd` (which use QR iteration, not bisection) or code that uses `Cheevx`/`Chegvx` with `abstol = 0` on matrices whose eigenvalues are within a "normal" floating-point range.
- The test driver (`Cdrvsg2stg`) sets `abstol = 2 × underflow`, which forces the `Rstebz`→`Cstein` code path in `Cheevx`. Typical user code with `abstol = 0` would take the `Rsterf`/`Csteqr` path first, avoiding this issue.

**Planned resolution**

- Root cause 1 (`itmax` clamping): deferred to a future release. A proper fix requires either adapting the bisection iteration limit to the actual precision/exponent-range ratio, or restructuring the test matrices to avoid eigenvalue magnitudes that exceed the feasible bisection range.
- Root cause 2 (error propagation): a defensive check will be added to `Cheevx` to skip the `Cstein` call when `Rstebz` reports non-convergence (negative `iblock` entries), preventing the spurious parameter error.

### 9.2 MPFR binary128 profile: Cgesvj (Jacobi SVD) test failures on type 4 matrices

**Symptom**

When running the complex SVD test suite (`Cdrvbd`) under the MPFR `binary128` profile (precision=113, emin=−16381, emax=16384), `Cgesvj` fails test 15 for all type 4 matrices ("evenly spaced singular values near underflow"). The failure is consistent across all matrix sizes and workspace variants:

| NB settings | Failures | Test ratio |
|-------------|----------|------------|
| All 5 NB/NBMIN/NX combinations | 28 / 14340 each | +5.1922968585348276e+33 |

The test ratio equals exactly 1/ulp = 2^112, which is the LAPACK sentinel for total decomposition failure. All other SVD drivers (`Cgesvd`, `Cgesdd`, `Cgesvdx`, `Cgesvdq`) and the `ZGESJV` driver pass the same tests. The `ZBD` bidiagonalization routines also pass completely.

**Root cause — MPFR lacks IEEE 754 subnormal (gradual underflow) emulation**

MPFR internally maintains a fully normalized mantissa at all times. Although the binary128 profile sets `emin=−16381` (corresponding to IEEE 754 binary128's smallest normal number ≈ 3.36×10^−4932), MPFR does not perform gradual underflow. Without explicit calls to `mpfr_subnormalize()` after every arithmetic operation, values below the smallest normal are flushed to zero instead of being represented as subnormal numbers with reduced precision.

IEEE 754 binary128 has a subnormal range spanning exponents [−16494, −16382), which covers 112 bits (approximately 34 decimal orders of magnitude). Any intermediate result landing in this range is:

- **IEEE 754 binary128**: a subnormal number (non-zero, with reduced but non-zero precision)
- **MPFR without `mpfr_subnormalize()`**: exactly zero

This divergence is catastrophic for `Cgesvj`, which uses Jacobi rotations with column-norm-based convergence criteria. Specifically:

1. **`Rlassq` (column norm computation)**: For type 4 matrices, some column entries fall below `smlnum`. IEEE 754 preserves these as subnormals, yielding small but non-zero column norms. MPFR flushes them to zero, producing zero norms for columns that should have non-zero norms.

2. **Jacobi rotation decisions**: `Cgesvj` compares products of the form `aapp * smlnum` against column norms `aaqq`. When `aapp < 1` (which is virtually always for near-underflow matrices), `aapp * smlnum` falls in the subnormal range — non-zero in IEEE 754 but zero in MPFR. This causes different branch paths in the rotation logic.

3. **Scaling via `Rlascl`**: When `Rlascl` is called with `cfrom` values that are themselves subnormal (valid in IEEE 754), MPFR sees `cfrom = 0` and returns `info = −4` (invalid input).

The fact that the test ratio is identical (exactly 1/ulp) for all matrix sizes and workspace settings confirms that the failure is deterministic and originates from a systematic flush-to-zero effect, not from a numerical precision issue.

**Why other SVD drivers are unaffected**

`Cgesvd` and `Cgesdd` use bidiagonalization followed by QR iteration or divide-and-conquer, which do not rely on column-norm comparisons in the subnormal range. `Cgesvdx` uses bisection on the bidiagonal matrix. These algorithms are structurally less sensitive to the subnormal/zero distinction because their intermediate quantities remain in the normal floating-point range after initial scaling.

**Scope and impact**

- Affects **only** the MPFR `binary128` profile. The MPFR default profile (512-bit mantissa) has a much smaller `smlnum` relative to its precision, so subnormal-sensitive code paths are not triggered by the test matrices.
- Does **not** affect `dd`, `qd`, `gmp`, native `binary128`, or native `binary80` backends.
- Does **not** affect user code that calls other SVD drivers (`Cgesvd`, `Cgesdd`, `Cgesvdx`, `Cgesvdq`), which all pass the same tests.
- The analogous real-valued routine `Rgesvj` is expected to exhibit the same failure under identical conditions.

**Planned resolution**

Two approaches are under consideration:

- **Option A — subnormal emulation in the MPFR binary128 wrapper**: Introduce `mpfr_subnormalize()` calls after each arithmetic operation in the MPFR wrapper layer when operating under the binary128 profile. This would produce IEEE 754–compliant behavior at the cost of significant performance degradation.
- **Option B — adjust `Cgesvj` scaling thresholds for MPFR**: Modify the underflow-sensitive comparisons in `Cgesvj` (and `Rgesvj`) to account for MPFR's flush-to-zero behavior when subnormal emulation is not active. This is less invasive but makes the algorithm MPFR-aware.

Neither fix is included in the 2.1.0 release. Users requiring Jacobi SVD on near-underflow matrices under the MPFR binary128 profile should use `Cgesvd` or `Cgesdd` as a workaround.

---

## 10. References (implementation context)

Primary commits that introduced the renames and mode selection:

* Rename floating-point types + I/O macro selection:
  [https://github.com/nakatamaho/mplapack/commit/3d70c98d324b8d302c4e4e179724cd1dfc1404a8](https://github.com/nakatamaho/mplapack/commit/3d70c98d324b8d302c4e4e179724cd1dfc1404a8)

* Type aliasing cleanup in public umbrella headers:
  [https://github.com/nakatamaho/mplapack/commit/fa64d200dd1a9f705767b372d07bc44171515a7d](https://github.com/nakatamaho/mplapack/commit/fa64d200dd1a9f705767b372d07bc44171515a7d)

* Mode-based architecture for binary80/binary128 and Rlamch-related refactor:
  [https://github.com/nakatamaho/mplapack/commit/f92d11d085e78055b6918aa1c863c4ad4008f7e3](https://github.com/nakatamaho/mplapack/commit/f92d11d085e78055b6918aa1c863c4ad4008f7e3)
