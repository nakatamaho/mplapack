# MPLAPACK 2.0.x → 2.1.0 Migration Guide

This document covers all changes in MPLAPACK 2.1.0 that require action from users or
downstream package maintainers.

---

## 1. Extended BLAS Routines Removed

**Impact:** Users who call extended BLAS routines directly, or call MPLAPACK routines that
depended on them.

### What changed

Extended BLAS routines have been removed from `libmpblas`. All MPLAPACK routines that
internally depended on extended BLAS have also been removed.

### Who needs to act

Anyone who:
- Links against `libmpblas` and calls extended BLAS routines directly.
- Links against `libmplapack` and calls LAPACK routines that were removed because they
  depended on extended BLAS.

### Migration

Check your source code for calls to extended BLAS routines. If you used them, you must
either implement equivalent functionality yourself or reorganize the algorithm to use the
standard BLAS subset.

```cpp
// Extended BLAS (REMOVED in 2.1.0) — example pattern
// Rgemm2(...);   // no longer available

// Standard BLAS alternative — still available
Rgemm("N", "N", m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
```

Compile your project against 2.1.0 and fix any undefined-reference linker errors that
appear for extended BLAS symbols.

---

## 2. Compiler Restrictions: oneAPI and Clang

**Impact:** Users building with Intel oneAPI (icx/icpx) or Clang/LLVM.

### What changed

| Compiler | binary128 | binary80 |
|---|---|---|
| Intel oneAPI (icx/icpx) | ❌ Removed | ❌ Removed |
| Clang/LLVM | ❌ Removed | ✅ Available |
| GCC | ✅ Available | ✅ Available |

These backends are not supported due to missing compiler-level support (128-bit floating-point
types, `__float128`, and related intrinsics are GCC-specific).

### Who needs to act

- Developers using Intel oneAPI who relied on `binary128` or `binary80` precision.
- Developers using Clang who relied on `binary128` precision.

### Migration

Switch to GCC for builds requiring `binary128` or (if using oneAPI) `binary80`.

```sh
# If using oneAPI and needing binary128/binary80: switch to GCC
export CC=gcc
export CXX=g++
./configure --enable-binary128=yes --enable-binary80=yes ...

# Clang: binary80 is available; binary128 is not
export CC=clang
export CXX=clang++
./configure --enable-binary80=yes ...
# Do NOT pass --enable-binary128=yes with Clang
```

If your build system detects the compiler and conditionally enables these backends, add a
GCC check before passing `--enable-binary128`:

```sh
if [ "$(${CXX} --version | head -1 | grep -c 'g++')" = "1" ]; then
    BINARY128_FLAG="--enable-binary128=yes"
fi
```

---

## 3. DD and QD Output Precision Reduced

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

## 4. Fable as Top-Level Standalone Component

**Impact:** Developers who build MPLAPACK from source, package maintainers, and anyone
who integrates the MPLAPACK source tree into larger build systems.

### What changed

`fable/` is now a top-level directory and an independent component of the repository.
It is no longer buried inside another subdirectory. The component includes:

- A modified version of Fable (Fortran-to-C++ converter)
- A modified version of FEM (Fortran Emulator, used in test programs only)
- `fable/3.9.1/`: 152 patches for LAPACK 3.9.1 sources that Fable cannot convert automatically
- `go.sh` and `go_testing.sh` driver scripts

### Who needs to act

- Packagers who enumerate source directories explicitly.
- Developers who maintain out-of-tree patches against the MPLAPACK directory structure.
- Build system integrators.

### Migration

Update any path references from the old location to `fable/` at the top level:

```sh
# Old: fable sources were not at top level
# New: top-level component
ls fable/
# go.sh  go_testing.sh  3.9.1/  ...
```

If you drive code generation yourself, use the provided scripts:

```sh
# Regenerate library routines (BLAS/LAPACK C++ sources + headers + patches)
cd ~/mplapack
bash fable/go.sh

# Regenerate test programs (EIG/LIN/MATGEN C++ sources + headers + patches)
cd ~/mplapack
bash fable/go_testing.sh
```

---

## 5. Auto-Generated Public Headers

**Impact:** Developers who manually maintain or patch `mpblas*.h` / `mplapack*.h` public headers.

### What changed

Public headers for MPBLAS and MPLAPACK are now generated by shell scripts
(`gen_include_mpblas.sh`, `gen_include_mplapack.sh`) rather than being hand-maintained.
The old template files `mplapack_utils__Float64x.h.in` and `mplapack_utils__Float128.h.in`
are no longer used.

### Who needs to act

- Anyone who patches the public headers directly and re-applies patches on each release.
- Anyone whose build system treats these headers as stable source files.

### Migration

Do not patch the generated headers directly. Patch the generator scripts or the template
inputs instead, so your changes survive regeneration.

```sh
# Do NOT edit these directly — they are overwritten by go.sh / go_testing.sh:
# include/mpblas_*.h
# include/mplapack_*.h

# Instead, modify the generator:
vi fable/gen_include_mpblas.sh
vi fable/gen_include_mplapack.sh

# Then regenerate:
bash fable/go.sh
```

---

## 6. Benchmark Directory Restructured

**Impact:** Developers or CI scripts that reference benchmark source paths directly.

### What changed

The benchmark subtree was flattened. `benchmark/mpblas/` and `benchmark/mplapack/`
subdirectories have been merged into `benchmark/`.

The `go.*.sh` benchmark runner scripts are now generated from `go.*.sh.in` templates at
build time and installed at install time. They are **not** source-controlled and must not
be added to `EXTRA_DIST`.

### Who needs to act

- CI pipelines that reference `benchmark/mpblas/` or `benchmark/mplapack/` paths.
- Developers who invoke benchmark runner scripts from their source tree rather than the
  installed prefix.

### Migration

```sh
# Old paths (2.0.x)
# benchmark/mpblas/go.Rgemm_double.sh
# benchmark/mplapack/go.Rgemm_mpfr.sh

# New paths (2.1.0) — generated at build time, installed to prefix
# benchmark/go.Rgemm_double.sh    (build tree, after make)
# $prefix/bin/go.Rgemm_mpfr.sh   (after make install)

# To run benchmarks after build:
make
./benchmark/go.Rgemm_mpfr.sh
```

---

## 7. Test Output Directory Layout Changed

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

# Update glob patterns to include the subdirectory level:
find mplapack/test/eig -name "*.out" | xargs grep FAIL
```

Use `summarize_mplapack_tests.py` for aggregated pass/fail reporting:

```sh
python3 misc/summarize_mplapack_tests.py
```

---

## 8. ILP32 Platform: `mplapackint` Changed to `long`

**Impact:** Users building on 32-bit platforms (Debian i386) or cross-compiling for ILP32.

### What changed

`mplapackint` is now unified with the platform `long` type. On ILP32 platforms this means
`mplapackint` is 32 bits. Associated `long`/`unsigned long` I/O overloads have been added to
prevent truncation.

### Who needs to act

Users building on or for 32-bit (ILP32) platforms who pass integer arrays to MPLAPACK routines.

### Migration

Verify that integer arrays passed to MPLAPACK routines match `mplapackint` (i.e., `long` on
your platform). On ILP32 platforms, do not assume 64-bit integers.

```cpp
// Check mplapackint size at compile time
#include <mplapack_config.h>
static_assert(sizeof(mplapackint) == sizeof(long),
              "mplapackint must match platform long");
```

---

## Summary of Required Actions by User Type

| User type | Actions required |
|---|---|
| Application developer using standard BLAS/LAPACK routines | None for most users; check for extended BLAS usage |
| Application developer using extended BLAS | Remove extended BLAS calls (§1) |
| Intel oneAPI user | Switch to GCC for binary128/binary80 (§2) |
| Clang user needing binary128 | Switch to GCC (§2) |
| User with DD/QD golden-value tests | Regenerate reference files (§3) |
| Packager / build system integrator | Update fable/ paths and benchmark paths (§4, §6) |
| Developer patching public headers | Patch generators instead of headers (§5) |
| CI pipeline collecting test output | Update output path globs (§7) |
| i386 / ILP32 platform developer | Verify integer type assumptions (§8) |
