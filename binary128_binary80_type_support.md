# MPLAPACK binary128 / binary80 Type Support Matrix

---

## Background: What Is Quadruple Precision (binary128)?

IEEE 754-2008 defines **binary128** ― commonly called *quadruple precision* or *quad precision* ―
as a 128-bit floating-point format with:

| Property | Value |
|----------|-------|
| Total bits | 128 (1 sign + 15 exponent + 112 mantissa fraction) |
| Decimal digits (approx.) | ~34 significant digits |
| Exponent range | approx. 10^－4932 to 10^4932 |
| C type (GCC  7, ISO/IEC TS 18661-3) | `_Float128` |
| Legacy GCC extension | `__float128` |
| `long double` on some ABIs | AArch64, s390x, MIPS (128-bit LD ABI), RISC-V |

Quadruple precision is approximately 2× the mantissa width of double precision (binary64,
~1517 digits), and is used in MPLAPACK as the highest fixed-precision backend before
arbitrary-precision MPFR/GMP kicks in. It is essential for numerical algorithms that require
error-free intermediate accumulation ― e.g., iterative refinement and eigenvalue verification.

### Implementation strategies in practice

| Config string | Meaning |
|---|---|
| `_Float128 + strfromf128` | Standard ISO C `_Float128` type with C23/glibc `strfromf128` for I/O. Preferred path on GCC  13 / glibc  2.29. |
| `__float128 + libquadmath` | GCC extension type with GCC's `libquadmath` for I/O and math. Required when `strfromf128` is absent (GCC  12, macOS Macports, MinGW, RHEL/Rocky  9, Intel oneAPI). |
| `long double (binary128)` | `long double` is IEEE binary128 on this ABI (AArch64, s390x, MIPS64le, RISC-V). No separate quad type; uses standard `long double` math and `%Lg` I/O. |

---

## Background: What Is Extended Precision (binary80)?

IEEE 754-2008 defines **binary80** ― commonly called *extended precision* or *80-bit extended* ―
as an 80-bit floating-point format originally introduced by Intel's x87 FPU:

| Property | Value |
|----------|-------|
| Total bits | 80 (1 sign + 15 exponent + 64 mantissa, with **explicit** integer bit) |
| Decimal digits (approx.) | ~1819 significant digits |
| Exponent range | approx. 10^－4932 to 10^4932 |
| C type (GCC  7, ISO/IEC TS 18661-3) | `_Float64x` |
| Legacy GCC extension | `__float80` |
| `long double` on some ABIs | x86 (IA-32) and x86_64 Linux/macOS/Windows (80-bit LD ABI) |
| Architecture availability | **x86 / x86_64 only** |

Unlike binary128, binary80 is a hardware format specific to the x87 FPU and has no software
fallback on other architectures. Its distinguishing feature is the *explicit* integer bit in
the mantissa ― absent in all other IEEE 754 binary formats ― which requires special handling
in low-level code. In MPLAPACK, binary80 fills the gap between double precision (~1517 digits)
and quadruple precision (~34 digits), and is particularly useful for avoiding catastrophic
cancellation in intermediate steps without the full cost of binary128.

### Implementation strategies in practice

| Config string | Meaning |
|---|---|
| `_Float64x (binary80)` | Standard ISO C `_Float64x` type with C23/glibc `strfromf64x` for I/O and `f64x` math. Preferred path on GCC  13 / glibc  2.29 on x86/x86_64. |
| `long double (binary80)` | `long double` is 80-bit extended on this platform (x86/x86_64 with legacy ABI or older GCC). Uses standard `long double` math and `%Lg` I/O via `snprintf`. |
| ― *(absent)* | Platform has no x87 hardware. binary80 is entirely unavailable. Affects all non-x86 architectures: AArch64, ppc64le, s390x, MIPS64le, RISC-V. |

---

## Table 1: All Platforms (OS-sorted, duplicates removed)

| Platform | Arch | Compiler | binary128 config | binary80 config |
|----------|------|----------|-----------------|-----------------|
| Alpine 3.19 | amd64 | GCC 13.2.1 (Alpine 13.2.1_git20231014) | `__float128` + libquadmath | `long double` (binary80) |
| Debian 12 | amd64 | GCC 14.2.0 | `_Float128` + strfromf128 | `_Float64x` (binary80) |
| Debian 12 | arm64 | GCC 12.2.0 | `long double` (binary128) | ― |
| Debian 12 | i386 | GCC 12.2.0 (Debian 12.2.0-14+deb12u1) | `__float128` + libquadmath | `long double` (binary80) |
| Debian 12 | mips64le | GCC 12.2.0 | `long double` (binary128) | ― |
| Debian 12 | ppc64le | GCC 12.2.0 (Debian 12.2.0-14+deb12u1) | `__float128` + libquadmath | ― |
| Debian 12 | s390x | GCC 12.2.0 | `long double` (binary128) | ― |
| Debian 13 | amd64 | GCC 14.2.0 | `_Float128` + strfromf128 | `_Float64x` (binary80) |
| Debian 13 | arm64 | GCC 14.2.0 | `_Float128` + strfromf128 | ― |
| Debian 13 | i686 | GCC 14.2.0 | `_Float128` + strfromf128 | `_Float64x` (binary80) |
| Debian 13 | ppc64le | GCC 14.2.0 | `_Float128` + strfromf128 | ― |
| Debian 13 | riscv64 | GCC 14.2.0 | `_Float128` + strfromf128 | ― |
| Debian 13 | s390x | GCC 14.2.0 | `_Float128` + strfromf128 | ― |
| macOS Sonoma | x86_64 | GCC-MP-14 (Macports) | `__float128` + libquadmath | `long double` (binary80) |
| Rocky 9 | amd64 | GCC 11.5.0 | `__float128` + libquadmath | *(config string not emitted)* |
| Rocky 9 | arm64 | GCC 11.5.0 | `long double` (binary128) | ― |
| Rocky 9 | ppc64le | GCC 11.5.0 | `__float128` + libquadmath | ― |
| Rocky 9 | s390x | GCC 11.5.0 | `long double` (binary128) | ― |
| Ubuntu 22.04 | amd64 | GCC 11.4.0 | `__float128` + libquadmath | `long double` (binary80) |
| Ubuntu 22.04 | amd64 | x86_64-w64-mingw32-gcc 10-win32 | `__float128` + libquadmath | `long double` (binary80) |
| Ubuntu 24.04 | amd64 | GCC 13.3.0 | `_Float128` + strfromf128 | `_Float64x` (binary80) |
| Ubuntu 24.04 | amd64 | Intel oneAPI 2025.3.2 (icpx) | `__float128` + libquadmath | `long double` (binary80) |
| Ubuntu 24.04 | amd64 | x86_64-w64-mingw32-gcc 13-win32 | `__float128` + libquadmath | `long double` (binary80) |
| Ubuntu 24.04 | arm64 | GCC 13.3.0 | `_Float128` + strfromf128 | ― |

> **Note (Rocky 9 amd64 / GCC 11.5.0):** All binary80 type-support probe lines are present
> (`__float80=yes`, `ldbl_binary80=yes`), but the final `checking binary80 configuration...`
> summary line is absent. This is a configure-script logging gap for GCC 11.x, not a hardware
> absence. The effective configuration is `long double (binary80)` by inference.
>
> **Note (Intel oneAPI 2025.3.2):** The `checking for binary128 type support...` probe line is
> absent ― the configure macro appears to be short-circuited for `icpx`. The binary128
> configuration was nevertheless resolved to `__float128 + libquadmath` via GCC interop, but
> the configure script should handle oneAPI explicitly rather than relying on silent fallthrough.

---

## Table 2: Unique (binary128 config, binary80 config) Combinations

| # | binary128 config | binary80 config | Representative platform(s) |
|---|-----------------|-----------------|---------------------------|
| A | `__float128` + libquadmath | `long double` (binary80) | Alpine 3.19 amd64 GCC 13; Debian 12 i386 GCC 12; macOS Sonoma GCC-MP-14; Rocky 9 amd64 GCC 11.5; Ubuntu 22.04 amd64 GCC 11 + MinGW 10; Ubuntu 24.04 amd64 MinGW 13 + Intel oneAPI 2025 |
| B | `__float128` + libquadmath | ― | Debian 12 ppc64le GCC 12; Rocky 9 ppc64le GCC 11.5 |
| C | `_Float128` + strfromf128 | `_Float64x` (binary80) | Debian 12 amd64 GCC 14; Debian 13 amd64/i686 GCC 14; Ubuntu 24.04 amd64 GCC 13 |
| D | `_Float128` + strfromf128 | ― | Debian 13 arm64/ppc64le/riscv64/s390x GCC 14; Ubuntu 24.04 arm64 GCC 13 |
| E | `long double` (binary128) | ― | Debian 12 arm64/mips64le/s390x GCC 12; Rocky 9 arm64/s390x GCC 11.5 |

---

## Technical Findings

### 1. GCC version is the primary determinant of the binary128 code path

- **GCC  12** and all non-Linux toolchains (Macports, MinGW, Intel oneAPI) →
  `__float128 + libquadmath`. These environments either lack `strfromf128` or supply it
  only in a libquadmath not detected by the configure probe.
- **GCC  13 on Linux** → `_Float128 + strfromf128`. The path is cleanly gated by
  `strfromf128` availability in glibc ( 2.29) and first-class `_Float128` support.
- Distro version is irrelevant: Debian 12 with backported GCC 14 lands in path **C**.

### 2. binary80 is strictly x86/x86_64-only

Without exception, every non-x86 architecture (arm64, ppc64le, s390x, mips64le, riscv64)
produces `―` for binary80. This is a hard ABI/hardware constraint. Any conditional
compilation relying on binary80 must be gated by architecture, not compiler version.

### 3. `long double` as binary128 ― five architectures, two GCC generations

Combination **E** covers AArch64, s390x, MIPS64le, RISC-V, and ppc64le-with-128-bit-LD-ABI
on GCC  12. In all these cases the ABI mandates a 128-bit `long double`, so binary128
arithmetic is available via standard C without any extension. On GCC 14 (Debian 13), these
architectures gain `_Float128 + strfromf128` (combination **D**), providing cleaner type
separation from `long double`.

### 4. ppc64le: libquadmath present, binary80 absent ― a portability trap

Combination **B** (Debian 12 and Rocky 9 ppc64le) is the only case where
`__float128 + libquadmath` is selected but binary80 is absent. ppc64le's `long double`
is either IBM double-double or IEEE binary128; neither is binary80. Any build logic that
assumes "presence of libquadmath implies binary80" will silently fail here.

### 5. MinGW always falls back to combination A regardless of host GCC

Both GCC 10-win32 (Ubuntu 22.04) and GCC 13-win32 (Ubuntu 24.04) resolve to
`__float128 + libquadmath` + `long double (binary80)`. The Windows target runtime does
not provide `strfromf128`, so combination **C** is never reachable from a MinGW build.

### 6. Intel oneAPI 2025: binary128 type probe skipped

`icpx` does not emit the `_Float128 type support` configure probe ― the check is
short-circuited. The fallback to `__float128 + libquadmath` succeeds via GCC runtime
interop, but the configure script should handle oneAPI explicitly rather than relying on
silent fallthrough.

### 7. RISC-V (riscv64) confirmed as combination D

Debian 13 riscv64 GCC 14.2.0 behaves identically to arm64/s390x on the same GCC version:
`_Float128 + strfromf128`, no binary80. RISC-V's `long double` is also IEEE binary128,
consistent with the other 128-bit-LD architectures.

### 8. Rocky 9 mirrors the RHEL 9 / GCC 11 ecosystem

All Rocky 9 entries (amd64, arm64, ppc64le, s390x) use GCC 11.5.0 and consistently fall
into old-style paths (combinations **A**, **B**, **E**). Enterprise distributions tracking
GCC 11.x for their full release cycle do not ship `strfromf128` in glibc and will remain
on the libquadmath path indefinitely without toolchain backports.

