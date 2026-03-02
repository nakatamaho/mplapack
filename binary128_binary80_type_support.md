# MPLAPACK binary128 / binary80 Type Support Matrix

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

## Table 1: All Platforms

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

> **Note (Rocky 9 amd64 / GCC 11.5.0):** All binary80 type-support probe lines are present in the
> log (`__float80=yes`, `ldbl_binary80=yes`, I/O and math checks), but the final
> `checking binary80 configuration...` summary line is absent. This is a logging/configure
> script gap for GCC 11.x, not a hardware absence. The effective configuration is
> `long double (binary80)` by inference.
>
> **Note (Intel oneAPI 2025.3.2):** The `checking for binary128 type support...` probe line is
> absent from the log ― the configure macro appears to have been skipped or short-circuited
> for `icpx`. The binary128 configuration was nevertheless resolved to `__float128 + libquadmath`,
> indicating that `_Float128` is not available but `__float128` (via GCC interop) is.

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
  `__float128 + libquadmath`. These environments either lack `strfromf128` entirely or
  supply it only in a libquadmath that is not detected by the configure probe.
- **GCC  13 on Linux** → `_Float128 + strfromf128`. The combination **C/D** path is
  cleanly gated by `strfromf128` availability in glibc ( 2.29) and GCC built-in
  support for `_Float128` as a first-class ISO type.
- Distro version is irrelevant: Debian 12 with a backported GCC 14 lands in path **C**,
  not **A**.

### 2. binary80 is strictly x86/x86_64-only

Without exception, every non-x86 architecture (arm64, ppc64le, s390x, mips64le, riscv64)
produces `―` for binary80. This is a hard ABI constraint, not a GCC version or distro issue.
Any conditional compilation relying on binary80 must be gated by architecture, not compiler.

### 3. `long double` as binary128 ― five architectures, two GCC generations

Combination **E** (`long double (binary128)`, no binary80) covers AArch64, s390x, MIPS64le,
and RISC-V on GCC  12 (Debian 12, Rocky 9), as well as AArch64 on Rocky 9. In all these
cases the ABI mandates a 128-bit `long double`, so binary128 arithmetic is available via
standard C without any extension. On GCC 14 (Debian 13), these same architectures gain
`_Float128 + strfromf128` (combination **D**), providing cleaner type separation.

### 4. ppc64le: libquadmath present, binary80 absent ― a portability trap

Combination **B** (Debian 12 and Rocky 9 ppc64le) is the only situation where
`__float128 + libquadmath` is selected but binary80 is unavailable. ppc64le's `long double`
is either IBM double-double or IEEE binary128 depending on the ABI; neither is binary80.
Any code or build logic that assumes "presence of libquadmath implies binary80" will silently
break on this architecture.

### 5. MinGW always falls back to combination A regardless of host GCC

Both GCC 10-win32 (Ubuntu 22.04) and GCC 13-win32 (Ubuntu 24.04) MinGW cross-compiler
targets resolve to `__float128 + libquadmath` + `long double (binary80)`. The Windows
target runtime does not provide `strfromf128`, so the `_Float128 + strfromf128` path is
never reachable from a MinGW build regardless of how new the host GCC is.

### 6. Intel oneAPI 2025: binary128 type probe skipped

`icpx` does not emit the `_Float128 type support` configure probe at all ― the check is
apparently short-circuited. The fallback to `__float128 + libquadmath` succeeds via GCC
runtime interop, but the absence of the type-support probe is a configure-script gap that
should be explicitly handled for oneAPI rather than relying on silent fallthrough.

### 7. RISC-V (riscv64) confirmed as combination D

Debian 13 riscv64 GCC 14.2.0 behaves identically to arm64/s390x on the same GCC version:
`_Float128 + strfromf128`, no binary80. RISC-V's `long double` is also IEEE binary128,
consistent with the other 128-bit-LD architectures.

### 8. Rocky 9 mirrors RHEL 9 / older GCC ecosystem

All Rocky 9 entries (amd64, arm64, ppc64le, s390x) use GCC 11.5.0 and consistently fall
into the old-style paths (combinations **A**, **B**, **E**). This is the expected behavior
for enterprise-class distributions that track GCC 11.x for their entire release lifecycle
and do not ship `strfromf128` in glibc.

