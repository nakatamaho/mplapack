# MPLAPACK Migration Guide

## MPLAPACK 2.3.0: Unified libraries

MPLAPACK 2.3.0 replaces the separate MPBLAS and MPLAPACK products with one
self-contained library for each backend/flavor. This is an ABI-breaking change;
the shared-library ABI is now 2.

### Link-line translation

| Old link line | New link line |
|---|---|
| `-lmplapack_dd -lmpblas_dd` | `-lmplapack_dd` |
| `-lmpblas_dd_opt` | `-lmplapack_dd_opt` |
| `-lmpblas_dd_opt -lmpblas_dd_cuda` | `-lmplapack_dd_opt_cuda` |
| `-lmpblas_binary128_opt_opencl` (+ opt library) | `-lmplapack_binary128_opt_opencl` |

Link exactly one `mplapack` library per backend/flavor. The selected library
already contains both the BLAS and LAPACK routines and all lower-priority
implementations needed by that flavor.

```sh
# DD reference flavor
c++ app.cpp -lmplapack_dd -lqd

# DD optimized CUDA flavor
c++ app.cpp -lmplapack_dd_opt_cuda -lqd -lcudart -lcublas
```

### CMake target translation

| Before 2.3.0 | MPLAPACK 2.3.0 |
|---|---|
| `mplapack::mpblas_<b>` | Gone; use `mplapack::mplapack_<b>` |
| `mplapack::mpblas_<b>_opt` | Gone; use `mplapack::mplapack_<b>_opt` |
| `mplapack::mpblas_dd_opt_cuda` | Gone; use `mplapack::mplapack_dd_opt_cuda` |
| `mplapack::mpblas_binary128_opt_opencl` | Gone; use `mplapack::mplapack_binary128_opt_opencl` |

The same one-target rule applies to CMake. Per-flavor pkg-config files are named
`mplapack_<flavor>.pc`; the aggregate `mplapack.pc` and separate `mpblas`
products are no longer installed.

---

This document covers all changes in MPLAPACK 2.1.0 that require action from users or
downstream package maintainers.

---

## 1. C++17 Now Required

**Impact:** All users building MPLAPACK from source.

### What changed

MPLAPACK 2.1.0 requires C++17 or later. The codebase uses C++17 features including
`std::optional`, structured bindings, `if constexpr`, and `std::atomic` with aggregate
initialization. These are not available in C++14 or earlier.

In 2.0.x, C++11/14 was sufficient for most build configurations.

### Who needs to act

Anyone building MPLAPACK with GCC older than 7, or with any compiler that does not
support C++17.

### Migration

Upgrade to GCC 7 or later (GCC 9+ recommended). The configure step will fail with an
explicit error if C++17 is not available.

```sh
# Verify your compiler supports C++17
g++ --version   # must be 7.0 or later
g++ -std=c++17 -x c++ - <<'EOF'
#include <optional>
int main() { std::optional<int> x = 1; return 0; }
EOF
echo "C++17 OK"

# Configure as usual — C++17 is enabled automatically
export CXX=g++ CC=gcc FC=gfortran
./configure --prefix=$HOME/MPLAPACK ...
```

On macOS, ensure the FSF GCC from MacPorts is used, not Apple Clang:

```sh
export CXX=g++-mp-14
export CC=gcc-mp-14
export FC=gfortran-mp-14
```

---

## 2. MPFR Exponent Range Now Auto-Adjusted with Precision

**Impact:** Users who set `MPLAPACK_MPFR_PRECISION` via environment variable, or who rely
on specific MPFR exponent range (emin/emax) values at runtime.

### What changed

When `MPLAPACK_MPFR_PRECISION` is set, MPLAPACK now automatically adjusts `emin`/`emax`
proportionally to the requested precision:

```
emax =  prec * 64   (clamped to mpfr_get_emax_max())
emin = -prec * 64   (clamped to mpfr_get_emin_min())
```

This adjustment is applied only when no explicit exponent range is given via environment
variables. If `MPLAPACK_MPFR_EMIN` / `MPLAPACK_MPFR_EMAX` are set explicitly, they take
precedence and no automatic adjustment is made.

In 2.0.x, the exponent range was not adjusted when precision changed, which could cause
underflow or overflow in high-precision computations where the default MPFR exponent range
was too narrow for the requested precision.

### Who needs to act

- Users who set `MPLAPACK_MPFR_PRECISION` and relied on the old (fixed) emin/emax.
- Users whose code assumes specific emin/emax values after initialization.

### Migration

If you need to preserve the old fixed exponent range regardless of precision, set
`MPLAPACK_MPFR_EMIN` and `MPLAPACK_MPFR_EMAX` explicitly to prevent auto-adjustment:

```sh
# Disable auto-adjustment by pinning the exponent range explicitly
export MPLAPACK_MPFR_PRECISION=256
export MPLAPACK_MPFR_EMIN=-1073741823    # mpfr_get_emin_min() on most platforms
export MPLAPACK_MPFR_EMAX=1073741823     # mpfr_get_emax_max() on most platforms
```

If you want the new auto-adjusted range (recommended), simply set precision and let
MPLAPACK compute emin/emax automatically:

```sh
export MPLAPACK_MPFR_PRECISION=256
# emin/emax will be set to ±16384 (= ±256*64) automatically
```

To inspect the effective exponent range at runtime:

```cpp
#include <mpfr.h>
mpfr_exp_t emin = mpfr_get_emin();
mpfr_exp_t emax = mpfr_get_emax();
```

---

## 3. binary128 and binary80 Library and Type Names Changed

**Impact:** All users who link against the `_Float128` or `_Float64x` precision backends,
use their types directly, or reference them in build scripts, configure flags, or headers.

### What changed

All names have been unified under the `binary128` / `binary80` naming scheme to abstract
away platform-specific type names and improve portability.

#### Configure flags

| 2.0.x | 2.1.0 |
|---|---|
| `--enable-_Float128=yes` | `--enable-binary128=yes` |
| `--enable-_Float64x=yes` | `--enable-binary80=yes` |

#### Library names

| 2.0.x | 2.1.0 |
|---|---|
| `libmplapack__Float128` | `libmplapack_binary128` |
| `libmpblas__Float128` | `libmpblas_binary128` |
| `libmplapack__Float64x` | `libmplapack_binary80` |
| `libmpblas__Float64x` | `libmpblas_binary80` |

#### Header names

| 2.0.x | 2.1.0 |
|---|---|
| `mplapack__Float128.h` | `mplapack_binary128.h` |
| `mpblas__Float128.h` | `mpblas_binary128.h` |
| `mplapack__Float64x.h` | `mplapack_binary80.h` |
| `mpblas__Float64x.h` | `mpblas_binary80.h` |

#### Type names

The platform-specific types (`_Float128`, `__float128`, `_Float64x`, `long double` used
as extended precision) are now abstracted by typedefs. Use these in all new code:

| 2.0.x (platform-specific) | 2.1.0 |
|---|---|
| `_Float128` / `__float128` | `mplapack_binary128_t` |
| `_Float64x` / `long double` (extended) | `mplapack_binary80_t` |

`mplapack_binary128_t` and `mplapack_binary80_t` are defined in `mplapack_config.h` and
resolve to the correct underlying type for each platform and compiler, allowing code to be
written once without `#ifdef` guards for different compilers or platforms.

### Who needs to act

Anyone who:
- Passes `--enable-_Float128` or `--enable-_Float64x` to `./configure`.
- Links against `libmplapack__Float128`, `libmpblas__Float128`, `libmplapack__Float64x`,
  or `libmpblas__Float64x` by name (in Makefiles, CMakeLists, pkg-config, etc.).
- Includes `mplapack__Float128.h`, `mpblas__Float128.h`, `mplapack__Float64x.h`, or
  `mpblas__Float64x.h` directly.
- Uses `_Float128`, `__float128`, or `_Float64x` as the scalar type in MPLAPACK calls.

### Migration

#### Build system / configure

```sh
# Old (2.0.x)
./configure --enable-_Float128=yes --enable-_Float64x=yes ...

# New (2.1.0)
./configure --enable-binary128=yes --enable-binary80=yes ...
```

#### Linking

```sh
# Old (2.0.x)
-lmplapack__Float128 -lmpblas__Float128
-lmplapack__Float64x -lmpblas__Float64x

# New (2.1.0)
-lmplapack_binary128 -lmpblas_binary128
-lmplapack_binary80  -lmpblas_binary80
```

#### Headers

```cpp
// Old (2.0.x)
#include <mplapack__Float128.h>
#include <mpblas__Float128.h>
#include <mplapack__Float64x.h>
#include <mpblas__Float64x.h>

// New (2.1.0)
#include <mplapack_binary128.h>
#include <mpblas_binary128.h>
#include <mplapack_binary80.h>
#include <mpblas_binary80.h>
```

#### Type usage in application code

```cpp
// Old (2.0.x) — platform-specific, required #ifdef for portability
#if defined(__GNUC__)
    _Float128 x = 1.0Q;
#elif defined(__ICC)
    __float128 x = ...;
#endif

// New (2.1.0) — portable, no #ifdef needed
#include <mplapack_config.h>
mplapack_binary128_t x = 1.0Q;
mplapack_binary80_t  y = 1.0L;
```

---

## 4. DD and QD Output Precision Reduced

**Impact:** Users who parse or compare textual output from DD/QD MPLAPACK routines, or who
use reference output files for golden-value testing.

### What changed

| Backend | 2.0.x digits | 2.1.0 digits |
|---|---|---|
| DD | 32 | 30 |
| QD | 64 | 61 |

This was necessary to eliminate spurious last-digit mismatches caused by rounding differences
across MPFR versions during binary-to-decimal conversion.

### Who needs to act

- Anyone with test suites that compare DD/QD output character-for-character or digit-for-digit
  against stored reference values.
- Anyone whose downstream code parses the number of significant digits from DD/QD output.

### Migration

Update your reference/golden-value files by regenerating them with 2.1.0:

```sh
make check
# Collect new reference output from mplapack/test/compare/dd/ and mplapack/test/compare/qd/
```

If you maintain hardcoded precision constants in your own code:

```cpp
// Old (2.0.x)
// DD_PRECISION = 32;
// QD_PRECISION = 64;

// New (2.1.0)
// DD_PRECISION = 30;
// QD_PRECISION = 61;
```

---

## 5. Compiler Restrictions: Clang

**Impact:** Users building with Intel oneAPI (icx/icpx) or Clang/LLVM.

### What changed

| Compiler | binary128 | binary80 |
|---|---|---|
| Clang/LLVM | ❌ Not supported | ✅ Available |
| GCC | ✅ Available | ✅ Available |
| Intel oneAPI (icx/icpx) | ✅ Available | ✅ Available |

`_Float128` / `__float128` and related intrinsics are GCC-specific and are not reliably
supported by other compilers.

### Who needs to act

- Developers using Intel oneAPI who rely on `binary128` or `binary80` precision.
- Developers using Clang who rely on `binary128` precision.

### Migration

Switch to GCC for builds requiring `binary128` or (for oneAPI users) either backend.

```sh
# If using oneAPI: switch to GCC for any binary128/binary80 work
export CC=gcc
export CXX=g++
./configure --enable-binary128=yes --enable-binary80=yes ...

# Clang: binary80 is available; binary128 is not
export CC=clang
export CXX=clang++
./configure --enable-binary80=yes ...
# Do NOT pass --enable-binary128=yes with Clang
```

If your build system detects the compiler and conditionally enables these backends:

```sh
if [ "$(${CXX} --version | head -1 | grep -c 'g++')" = "1" ]; then
    BINARY128_FLAG="--enable-binary128=yes"
fi
```

---

## 6. Fable Conversion Pipeline: Not in Release Tarball

**Impact:** Developers who want to regenerate BLAS/LAPACK C++ sources or test programs
from Fortran, or who build MPLAPACK directly from a Git clone.

### What changed

`fable/` is a top-level standalone component of the repository that provides automated
Fortran-to-C++ conversion. However, **the Fable conversion pipeline is not included in the
release tarball** (`mplapack-2.1.0.tar.xz`). The tarball contains only the pre-generated
C++ sources.

To use Fable (e.g., to regenerate sources or apply your own patches), you must:

1. Clone the repository directly from Git.
2. Expand the bundled LAPACK 3.9.1 source under `external/lapack/` before running the
   conversion scripts.

### Who needs to act

- Developers who want to regenerate library or test C++ sources from Fortran.
- Developers maintaining out-of-tree patches against the Fortran sources.
- Anyone who previously used the Fable Docker build environment.

### Migration

```sh
# Step 1: clone the repository (tarball is not sufficient)
git clone https://github.com/nakatamaho/mplapack
cd mplapack

# Step 2: expand the bundled LAPACK sources under external/
cd external/lapack
tar xvf lapack-3.9.1.tar.gz    # or equivalent expansion step
cd ../..

# Step 3: run the conversion pipeline
bash fable/go.sh          # library routines (BLAS/LAPACK C++ sources + headers + patches)
bash fable/go_testing.sh  # test programs (EIG/LIN/MATGEN C++ sources + headers + patches)
```

The generated C++ sources are placed in their respective directories under `mpblas/` and
`mplapack/`. The 152 patches under `fable/3.9.1/` are applied automatically by the scripts.

---

## 7. Auto-Generated MPBLAS/MPLAPACK Public Headers (Subset)

**Impact:** Developers who include or patch the public headers for GMP, MPFR, or
test-category-specific backends (EIG, LIN, MATGEN).

### What changed

A subset of public headers is now auto-generated by shell scripts rather than being
hand-maintained. The auto-generated headers are:

- `mpblas_gmp.h`, `mplapack_gmp.h`
- Headers for each combination of test category (EIG, LIN, MATGEN) and precision backend

Other headers such as `mplapack_utils_*.h` are **not** auto-generated and remain
hand-maintained source files.

### Who needs to act

- Developers who patch `mpblas_gmp.h`, `mplapack_gmp.h`, or the test-category headers
  directly and re-apply patches on each release.

### Migration

Do not patch auto-generated headers directly — they will be overwritten. Patch the
generator scripts instead:

```sh
# Auto-generated (do NOT patch directly):
# include/mpblas_gmp.h
# include/mplapack_gmp.h
# include/mplapack_{eig,lin,matgen}_*.h  (per-precision variants)

# Patch the generators instead:
vi fable/gen_include_mpblas.sh
vi fable/gen_include_mplapack.sh

# Then regenerate (requires Git clone + external/lapack expanded; see §5):
bash fable/go.sh
```

---

## 8. Benchmark Directory Restructured

**Impact:** Developers or CI scripts that reference benchmark source paths directly.

### What changed

The benchmark subtree was flattened. `benchmark/mpblas/` and `benchmark/mplapack/`
subdirectories have been merged into `benchmark/`.

The `go.*.sh` benchmark runner scripts are now generated from `go.*.sh.in` templates at
build time and installed at install time. They are **not** source-controlled and must not
be added to `EXTRA_DIST`.

### Who needs to act

- CI pipelines that reference `benchmark/mpblas/` or `benchmark/mplapack/` paths.
- Developers who invoke benchmark runner scripts from the source tree.

### Migration

```sh
# Old paths (2.0.x)
# benchmark/mpblas/go.Rgemm_double.sh
# benchmark/mplapack/go.Rgemm_mpfr.sh

# New paths (2.1.0) — generated at build time, installed to prefix
# benchmark/go.Rgemm_double.sh    (build tree, after make)
# $prefix/bin/go.Rgemm_mpfr.sh   (after make install)

make
./benchmark/go.Rgemm_mpfr.sh
```

---

## 9. Test Output Directory Layout Changed

**Impact:** CI pipelines and scripts that parse or archive test output files by path.

### What changed

Test output files are now written into a `CPU-OS-Compiler` subdirectory under each backend's
test directory, determined at runtime by `cpu_os_detection.sh`. The old flat layout is gone.

### Who needs to act

CI scripts or post-processing tools that collect test output files by hardcoded path.

### Migration

```sh
# Old layout (2.0.x)
# mplapack/test/eig/mpfr/eigtst_mpfr.out

# New layout (2.1.0)
# mplapack/test/eig/mpfr/<CPU>-<OS>-<Compiler>/eigtst_mpfr.out
# e.g.: mplapack/test/eig/mpfr/Ryzen3970X-Ubuntu2204-GCC1140/eigtst_mpfr.out

find mplapack/test/eig -name "*.out" | xargs grep FAIL
```

Use `summarize_mplapack_tests.py` for aggregated pass/fail reporting:

```sh
python3 misc/summarize_mplapack_tests.py
```

---

## 10. ILP32/LP64: `mplapackint` Tracks Platform `long`

**Impact:** Users building on or for 32-bit platforms (Debian i386, ILP32) or LP64
platforms where the integer width matters for GMP/MPFR interoperability.

### What changed

`mplapackint` is now unified with the platform `long` type:

- On **ILP32** platforms (e.g., Debian i386): `long` is 32 bits → `mplapackint` is 32 bits.
- On **LP64** platforms (e.g., Linux amd64/arm64): `long` is 64 bits → `mplapackint` is 64 bits.

**Why `long` specifically:** The GMP and MPFR C APIs accept integers only as `long` (e.g.,
`mpfr_set_si`, `mpz_get_si`). Neither `gmpxx.h` nor `mpreal.h` provide overloads for
`int64_t` or `int32_t` independently of `long`. Unifying `mplapackint` with `long` ensures
that MPLAPACK integer arguments can be passed directly to GMP/MPFR without conversion,
and that the correct width is used on both ILP32 and LP64 without custom wrappers at the
call site.

Associated `long`/`unsigned long` I/O overloads have been added to prevent truncation on
ILP32.

### Who needs to act

- Users building on or for ILP32 platforms who pass integer arrays to MPLAPACK routines.
- Developers who assumed `mplapackint` is always 64 bits.

### Migration

Do not assume `mplapackint` is 64 bits. Use `mplapackint` (or `long`) consistently for
all integer arguments passed to MPLAPACK/MPBLAS routines.

```cpp
// Do NOT assume 64-bit integers
// int64_t n = 1000;       // may truncate on ILP32
// int      n = 1000;      // may truncate on LP64 if > INT_MAX

// Correct: use mplapackint (= long on your platform)
#include <mplapack_config.h>
mplapackint n = 1000;

// Verify the width at compile time if needed:
static_assert(sizeof(mplapackint) == sizeof(long),
              "mplapackint must match platform long");
```

---

## 11. Known Issue: lin/dd on AArch64 (ARM)

**Impact:** Users running the `lin/dd` test suite on AArch64 hardware (e.g., Raspberry Pi 4,
AWS Graviton, Apple Silicon with GCC).

### What is observed

On Raspberry Pi 4 (Cortex-A72, AArch64, Ubuntu 24.04, GCC 13.3.0), the `lin/dd` test suite
produces 226 failures (0.03%), compared to ≤ 12 (0.00%) on all tested amd64 platforms.
All other backends and `eig/dd` on the same hardware are within the normal range.

This is **not** a regression in MPLAPACK routines themselves. The root cause is under
investigation and is suspected to be a GCC bug on AArch64 affecting double-double arithmetic.

See: https://github.com/nakatamaho/mplapack/issues/76

### Who is affected

Users who:
- Build and run `make check` with `--enable-dd=yes` on AArch64 hardware.
- Use `lin/dd` test results as a pass/fail gate in CI on AArch64.

### Workaround

There is no fix available in 2.1.0. If you encounter this failure pattern on AArch64,
it is expected and tracked. You may suppress the `lin/dd` failure gate on AArch64 in CI
until the issue is resolved.

```sh
# After make check, confirm with the summarizer:
python3 misc/summarize_mplapack_tests.py
# Expected anomaly on AArch64:
# lin  dd  865775  226  0.03%  2  4
```

---

## Summary of Required Actions by User Type

| User type | Actions required |
|---|---|
| All users building from source | Upgrade to GCC 7+ (C++17 required) (§1) |
| User setting `MPLAPACK_MPFR_PRECISION` | Verify emin/emax behavior; pin explicitly if needed (§2) |
| User of binary128/binary80 backends | Update configure flags, library names, headers, types (§3) |
| User with DD/QD golden-value tests | Regenerate reference files (§4) |
| Intel oneAPI user | Switch to GCC; binary128 and binary80 both broken (§5) |
| Clang user needing binary128 | Switch to GCC (§5) |
| Developer regenerating C++ sources from Fortran | Use Git clone + expand external/lapack (§6) |
| Developer patching auto-generated headers | Patch generators instead of headers (§7) |
| CI pipeline referencing benchmark paths | Update paths under benchmark/ (§8) |
| CI pipeline collecting test output | Update output path globs (§9) |
| ILP32 / LP64 platform developer | Use `mplapackint` (= `long`) for all integer args (§10) |
| AArch64 user running lin/dd tests | Failures are expected and tracked; see §11 |
