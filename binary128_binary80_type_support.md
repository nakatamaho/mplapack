# MPLAPACK binary128 / binary80 Type Support Matrix

## Table 1: All Platforms

| Platform | Arch | Compiler | binary128 config | binary80 config |
|----------|------|----------|-----------------|-----------------|
| Alpine 3.19 | x86_64 | GCC 13.2.1 (Alpine 13.2.1_git20231014) | `__float128` + libquadmath | `long double` (binary80) |
| Debian 12 | amd64 | GCC 14.2.0 | `_Float128` + strfromf128 | `_Float64x` (binary80) |
| Debian 12 | arm64 | GCC 12.2.0 | `long double` (binary128) | ― |
| Debian 12 | i386 | GCC 12.2.0 (Debian 12.2.0-14+deb12u1) | `__float128` + libquadmath | `long double` (binary80) |
| Debian 12 | ppc64le | GCC 12.2.0 (Debian 12.2.0-14+deb12u1) | `__float128` + libquadmath | ― |
| Debian 12 | ppc64le | GCC 12.2.0 *(duplicate entry ― identical result)* | `__float128` + libquadmath | ― |
| Debian 12 | s390x | GCC 12.2.0 | `long double` (binary128) | ― |
| Debian 13 | amd64 | GCC 14.2.0 | `_Float128` + strfromf128 | `_Float64x` (binary80) |
| Debian 13 | arm64 | GCC 14.2.0 | `_Float128` + strfromf128 | ― |
| Debian 13 | i686 | GCC 14.2.0 | `_Float128` + strfromf128 | `_Float64x` (binary80) |
| Debian 13 | ppc64le | GCC 14.2.0 | `_Float128` + strfromf128 | ― |
| Debian 13 | s390x | GCC 14.2.0 | `_Float128` + strfromf128 | ― |
| macOS Sonoma | x86_64 | GCC-MP-14 (Macports) | `__float128` + libquadmath | `long double` (binary80) |
| Ubuntu 22.04 | amd64 | GCC 11.4.0 | `__float128` + libquadmath | *(config string not emitted)* |
| Ubuntu 22.04 | amd64 | x86_64-w64-mingw32-gcc 10-win32 | `__float128` + libquadmath | `long double` (binary80) |
| Ubuntu 24.04 | amd64 | GCC 13.3.0 | `_Float128` + strfromf128 | `_Float64x` (binary80) |
| Ubuntu 24.04 | amd64 | x86_64-w64-mingw32-gcc 13-win32 | `__float128` + libquadmath | `long double` (binary80) |
| Ubuntu 24.04 | arm64 | GCC 13.3.0 | `_Float128` + strfromf128 | ― |

> **Note (Ubuntu 22.04 / GCC 11.4.0):** All binary80 type-support probe lines are present in
> the log (`__float80=yes`, `ldbl_binary80=yes`, I/O and math checks), but the final
> `checking binary80 configuration...` summary line is absent. This is a logging gap in the
> configure script for this compiler version, not an absence of binary80 hardware support.
>
> **Note (Debian 12 ppc64le):** Two log entries exist for this platform/arch/compiler combination
> with identical probe results. The first is labeled with the full Debian package string
> (`gcc (Debian 12.2.0-14+deb12u1) 12.2.0`); the second simply as `gcc-12.2.0`. This is likely
> a duplicate run in the test matrix and carries no distinct information.

---

## Table 2: Unique (binary128 config, binary80 config) Combinations

| # | binary128 config | binary80 config | Representative platform(s) |
|---|-----------------|-----------------|---------------------------|
| A | `__float128` + libquadmath | `long double` (binary80) | Ubuntu 22.04 amd64 / MinGW-w64 GCC 10; Alpine 3.19 x86_64 GCC 13; macOS Sonoma x86_64 GCC-MP-14; Debian 12 i386 GCC 12; Ubuntu 24.04 amd64 / MinGW-w64 GCC 13 |
| B | `__float128` + libquadmath | ― | Debian 12 ppc64le GCC 12 |
| C | `_Float128` + strfromf128 | `_Float64x` (binary80) | Ubuntu 24.04 amd64 GCC 13; Debian 12 amd64 GCC 14; Debian 13 amd64 GCC 14; Debian 13 i686 GCC 14 |
| D | `_Float128` + strfromf128 | ― | Ubuntu 24.04 arm64 GCC 13; Debian 13 arm64/ppc64le/s390x GCC 14 |
| E | `long double` (binary128) | ― | Debian 12 arm64 GCC 12; Debian 12 s390x GCC 12 |

---

## Technical Findings

### 1. GCC version is the primary determinant of the binary128 code path

The transition from `__float128 + libquadmath` (combinations **A/B**) to
`_Float128 + strfromf128` (combinations **C/D**) corresponds tightly to the GCC version:

- GCC  12 (and GCC-MP-14 on macOS, which lacks `strfromf128` in its libquadmath) →
  `__float128 + libquadmath`
- GCC  13 on Linux → `_Float128 + strfromf128`

The Debian 12 amd64 entry with backported GCC 14.2.0 already uses `_Float128 + strfromf128`
despite being Debian 12. This confirms the GCC version, not the distro version, drives the
selection.

### 2. binary80 availability is purely architecture-dependent

binary80 (x87 extended precision) is absent on every non-x86 architecture tested: arm64,
ppc64le, and s390x all show `―`. Any code path depending on binary80 must be guarded by an
architecture check, not an OS or compiler-version check.

### 3. AArch64 and s390x: `long double` on GCC 12, `_Float128` on GCC 14

- Debian 12 arm64/s390x GCC 12 → `long double (binary128)` (combination **E**)
- Debian 13 arm64/s390x GCC 14 → `_Float128 + strfromf128` (combination **D**)

Both ABIs map `long double` to IEEE 754 binary128, so 128-bit precision is available in both
cases. However, the I/O and math function sets differ, which affects portability of formatting
and special-function code.


### 5. MinGW cross-compilation always falls back to combination A

Both the GCC 10-win32 (Ubuntu 22.04) and GCC 13-win32 (Ubuntu 24.04) MinGW cross-compiler
targets produce `__float128 + libquadmath` + `long double (binary80)`, regardless of host GCC
generation. The Windows target runtime lacks `strfromf128`, so the configure script correctly
falls back to libquadmath even when the host GCC is 13.
