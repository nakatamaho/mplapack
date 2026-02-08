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

## 9. References (implementation context)

Primary commits that introduced the renames and mode selection:

* Rename floating-point types + I/O macro selection:
  [https://github.com/nakatamaho/mplapack/commit/3d70c98d324b8d302c4e4e179724cd1dfc1404a8](https://github.com/nakatamaho/mplapack/commit/3d70c98d324b8d302c4e4e179724cd1dfc1404a8)

* Type aliasing cleanup in public umbrella headers:
  [https://github.com/nakatamaho/mplapack/commit/fa64d200dd1a9f705767b372d07bc44171515a7d](https://github.com/nakatamaho/mplapack/commit/fa64d200dd1a9f705767b372d07bc44171515a7d)

* Mode-based architecture for binary80/binary128 and Rlamch-related refactor:
  [https://github.com/nakatamaho/mplapack/commit/f92d11d085e78055b6918aa1c863c4ad4008f7e3](https://github.com/nakatamaho/mplapack/commit/f92d11d085e78055b6918aa1c863c4ad4008f7e3)

