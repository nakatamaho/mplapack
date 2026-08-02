# MPLAPACK binary128 / binary80 Type Support Matrix

*Updated from the MPLAPACK 2.1.0 / 2.2.0-dev build matrix — 82 successful configurations
across Alpine, Debian, Ubuntu, Rocky, Fedora, and openSUSE on x86_64, i386/i686,
aarch64, ppc64le, s390x, riscv64, and mips64le, plus Intel oneAPI (icpx) and MinGW-w64
cross targets.*

---

## Background: What Is Quadruple Precision (binary128)?

IEEE 754-2008 defines **binary128** ― commonly called *quadruple precision* or *quad precision* ―
as a 128-bit floating-point format with:

| Property | Value |
|----------|-------|
| Total bits | 128 (1 sign + 15 exponent + 112 mantissa fraction) |
| Decimal digits (approx.) | ~34 significant digits |
| Exponent range | approx. 10^−4932 to 10^4932 |
| C type (GCC ≥ 7, ISO/IEC TS 18661-3) | `_Float128` |
| Legacy GCC extension | `__float128` |
| `long double` on some ABIs | AArch64, s390x, MIPS (128-bit LD ABI), RISC-V |

Quadruple precision is approximately 2× the mantissa width of double precision (binary64,
~15–17 digits), and is used in MPLAPACK as the highest fixed-precision backend before
arbitrary-precision MPFR/GMP kicks in. It is essential for numerical algorithms that require
error-free intermediate accumulation ― e.g., iterative refinement and eigenvalue verification.

### Implementation strategies in practice

| Config string | Meaning |
|---|---|
| `_Float128` + `strfromf128` | Standard ISO C `_Float128` type with C23/glibc `strfromf128` for I/O. Preferred path on **GCC ≥ 13** together with **glibc ≥ 2.26** (effectively ≥ 2.29). |
| `__float128` + `libquadmath` | GCC extension type with GCC's `libquadmath` for I/O and math. Required when `strfromf128` is absent: older GCC (≤ 12), **musl libc (Alpine, any GCC version)**, MinGW, older RHEL/Rocky, Intel oneAPI. |
| `long double (binary128)` | `long double` *is* IEEE binary128 on this ABI (AArch64, s390x, MIPS64le, RISC-V with older GCC). No separate quad type; uses standard `long double` math and `%Lg` I/O. Disappears once GCC ≥ 13 on glibc starts providing proper `_Float128`. |

---

## Background: What Is Extended Precision (binary80)?

IEEE 754-2008 defines **binary80** ― commonly called *extended precision* or *80-bit extended* ―
as an 80-bit floating-point format originally introduced by Intel's x87 FPU:

| Property | Value |
|----------|-------|
| Total bits | 80 (1 sign + 15 exponent + 64 mantissa, with **explicit** integer bit) |
| Decimal digits (approx.) | ~18–19 significant digits |
| Exponent range | approx. 10^−4932 to 10^4932 |
| C type (GCC ≥ 7, ISO/IEC TS 18661-3) | `_Float64x` |
| Legacy GCC extension | `__float80` |
| `long double` on some ABIs | x86 (IA-32) and x86_64 Linux/macOS/Windows (80-bit LD ABI) |
| Architecture availability | **x86 / x86_64 only** |

Unlike binary128, binary80 is a hardware format specific to the x87 FPU and has no software
fallback on other architectures. Its distinguishing feature is the *explicit* integer bit in
the mantissa ― absent in all other IEEE 754 binary formats ― which requires special handling
in low-level code. In MPLAPACK, binary80 fills the gap between double precision (~15–17 digits)
and quadruple precision (~34 digits), and is particularly useful for avoiding catastrophic
cancellation in intermediate steps without the full cost of binary128.

### Implementation strategies in practice

| Config string | Meaning |
|---|---|
| `_Float64x (binary80)` | Standard ISO C `_Float64x` type with glibc `strfromf64x` for I/O. Selected on **GCC ≥ 13 + glibc ≥ 2.29** on x86/x86_64. |
| `long double (binary80)` | `long double` is 80-bit extended on this platform. Uses standard `long double` math and `%Lg` I/O via `snprintf`. Selected on x86/x86_64 for: older GCC, **musl libc (Alpine)**, MinGW, Intel oneAPI. |
| ― *(absent / disabled)* | Platform has no x87 hardware. binary80 is entirely unavailable. Affects *all* non-x86 architectures: aarch64, ppc64le, s390x, mips64le, riscv64. |

---

## Table 1: All Successful Build Configurations

Compiler versions are shown compactly (`gcc 14.2`, `gcc 11.5`, `icpx 2025.3`, etc.).
Configurations grouped by OS family; within each family, rows sorted by OS version then architecture.

| OS / Release | Arch | CPU | IM | Compiler | binary128 config | binary80 config |
|---|---|---|---|---|---|---|
| **Alpine 3.19** | amd64 | x86_64 | LP64 | gcc 13.2 (musl) | `__float128` + libquadmath | `long double (binary80)` |
| **Alpine 3.19** | arm64 | aarch64 | LP64 | gcc 13.2 (musl) | `long double (binary128)` | ― |
| **Alpine 3.22** | i386 | x86_64 | ILP32 | gcc 14.2 (musl) | `__float128` + libquadmath | `long double (binary80)` |
| **Alpine 3.22** | amd64 | x86_64 | LP64 | gcc 14.2 (musl) | `__float128` + libquadmath | `long double (binary80)` |
| **Alpine 3.22** | arm64 | aarch64 | LP64 | gcc 14.2 (musl) | `long double (binary128)` | ― |
| **Alpine 3.22** | ppc64le | ppc64le | LP64 | gcc 14.2 (musl) | `__float128` + libquadmath | ― |
| **Alpine 3.22** | riscv64 | riscv64 | LP64 | gcc 14.2 (musl) | `long double (binary128)` | ― |
| **Alpine 3.22** | s390x | s390x | LP64 | gcc 14.2 (musl) | `long double (binary128)` | ― |
| **Alpine 3.23** | i386 | x86_64 | ILP32 | gcc 15.2 (musl) | `__float128` + libquadmath | `long double (binary80)` |
| **Alpine 3.23** | amd64 | x86_64 | LP64 | gcc 15.2 (musl) | `__float128` + libquadmath | `long double (binary80)` |
| **Alpine 3.23** | arm64 | aarch64 | LP64 | gcc 15.2 (musl) | `long double (binary128)` | ― |
| **Alpine 3.23** | ppc64le | ppc64le | LP64 | gcc 15.2 (musl) | `__float128` + libquadmath | ― |
| **Alpine 3.23** | riscv64 | riscv64 | LP64 | gcc 15.2 (musl) | `long double (binary128)` | ― |
| **Alpine 3.23** | s390x | s390x | LP64 | gcc 15.2 (musl) | `long double (binary128)` | ― |
| **Debian 11** | i386 | i686 | ILP32 | gcc 10.2 | `__float128` + libquadmath | `long double (binary80)` |
| **Debian 11** | amd64 | x86_64 | LP64 | gcc 10.2 | `__float128` + libquadmath | `long double (binary80)` |
| **Debian 11** | arm64 | aarch64 | LP64 | gcc 10.2 | `long double (binary128)` | ― |
| **Debian 12** | i386 | i686 | ILP32 | gcc 12.2 | `__float128` + libquadmath | `long double (binary80)` |
| **Debian 12** | amd64 | x86_64 | LP64 | gcc 12.2 | `__float128` + libquadmath | `long double (binary80)` |
| **Debian 12** | arm64 | aarch64 | LP64 | gcc 12.2 | `long double (binary128)` | ― |
| **Debian 12** | ppc64le | ppc64le | LP64 | gcc 12.2 | `__float128` + libquadmath | ― |
| **Debian 12** | s390x | s390x | LP64 | gcc 12.2 | `long double (binary128)` | ― |
| **Debian 12** | mips64le | mips64 | LP64 | gcc 12.2 | `long double (binary128)` | ― |
| **Debian 13** | i386 | i686 | ILP32 | gcc 14.2 | `_Float128` + strfromf128 | `_Float64x (binary80)` |
| **Debian 13** | amd64 | x86_64 | LP64 | gcc 14.2 | `_Float128` + strfromf128 | `_Float64x (binary80)` |
| **Debian 13** | arm64 | aarch64 | LP64 | gcc 14.2 | `_Float128` + strfromf128 | ― |
| **Debian 13** | ppc64le | ppc64le | LP64 | gcc 14.2 | `_Float128` + strfromf128 | ― |
| **Debian 13** | riscv64 | riscv64 | LP64 | gcc 14.2 | `_Float128` + strfromf128 | ― |
| **Debian 13** | s390x | s390x | LP64 | gcc 14.2 | `_Float128` + strfromf128 | ― |
| **Debian 13** + MinGW-w64 | amd64 (Windows tgt) | x86_64 | LLP64 | x86_64-w64-mingw32-g++ 14 | `__float128` + libquadmath | `long double (binary80)` |
| **Ubuntu 18.04** | i386 | i686 | ILP32 | gcc 7.5 | `__float128` + libquadmath | `long double (binary80)` |
| **Ubuntu 18.04** | amd64 | x86_64 | LP64 | gcc 7.5 | `__float128` + libquadmath | `long double (binary80)` |
| **Ubuntu 18.04** | arm64 | aarch64 | LP64 | gcc 7.5 | `long double (binary128)` | ― |
| **Ubuntu 20.04** | amd64 | x86_64 | LP64 | gcc 9.4 | `__float128` + libquadmath | `long double (binary80)` |
| **Ubuntu 20.04** | arm64 | aarch64 | LP64 | gcc 9.4 | `long double (binary128)` | ― |
| **Ubuntu 22.04** | amd64 | x86_64 | LP64 | gcc 11.4 | `__float128` + libquadmath | `long double (binary80)` |
| **Ubuntu 22.04** | arm64 | aarch64 | LP64 | gcc 11.4 | `long double (binary128)` | ― |
| **Ubuntu 22.04** | ppc64le | ppc64le | LP64 | gcc 11.4 | `__float128` + libquadmath | ― |
| **Ubuntu 22.04** | s390x | s390x | LP64 | gcc 11.4 | `long double (binary128)` | ― |
| **Ubuntu 22.04** + Intel oneAPI | amd64 | x86_64 | LP64 | icpx 2025.3.2 | `__float128` + libquadmath | `long double (binary80)` |
| **Ubuntu 22.04** + MinGW-w64 | amd64 (Windows tgt) | x86_64 | LLP64 | x86_64-w64-mingw32-g++ 10 | `__float128` + libquadmath | `long double (binary80)` |
| **Ubuntu 24.04** | amd64 | x86_64 | LP64 | gcc 13.3 | `_Float128` + strfromf128 | `_Float64x (binary80)` |
| **Ubuntu 24.04** | arm64 | aarch64 | LP64 | gcc 13.3 | `_Float128` + strfromf128 | ― |
| **Ubuntu 24.04** | ppc64le | ppc64le | LP64 | gcc 13.3 | `_Float128` + strfromf128 | ― |
| **Ubuntu 24.04** | s390x | s390x | LP64 | gcc 13.3 | `_Float128` + strfromf128 | ― |
| **Ubuntu 24.04** + Intel oneAPI | amd64 | x86_64 | LP64 | icpx 2025.3.2 | `__float128` + libquadmath | `long double (binary80)` |
| **Ubuntu 24.04** + MinGW-w64 | amd64 (Windows tgt) | x86_64 | LLP64 | x86_64-w64-mingw32-g++ 13 | `__float128` + libquadmath | `long double (binary80)` |
| **Rocky 8** | amd64 | x86_64 | LP64 | gcc 8.5 | `__float128` + libquadmath | `long double (binary80)` |
| **Rocky 8** | arm64 | aarch64 | LP64 | gcc 8.5 | `long double (binary128)` | ― |
| **Rocky 9** | amd64 | x86_64 | LP64 | gcc 11.5 | `__float128` + libquadmath | `long double (binary80)` |
| **Rocky 9** | arm64 | aarch64 | LP64 | gcc 11.5 | `long double (binary128)` | ― |
| **Rocky 9** | ppc64le | ppc64le | LP64 | gcc 11.5 | `__float128` + libquadmath | ― |
| **Rocky 9** | s390x | s390x | LP64 | gcc 11.5 | `long double (binary128)` | ― |
| **Rocky 10** | amd64 | x86_64 | LP64 | gcc 14.3 | `_Float128` + strfromf128 | `_Float64x (binary80)` |
| **Rocky 10** | arm64 | aarch64 | LP64 | gcc 14.3 | `_Float128` + strfromf128 | ― |
| **Rocky 10** | ppc64le | ppc64le | LP64 | gcc 14.3 | `_Float128` + strfromf128 | ― |
| **Rocky 10** | riscv64 | riscv64 | LP64 | gcc 14.3 | `_Float128` + strfromf128 | ― |
| **Rocky 10** | s390x | s390x | LP64 | gcc 14.3 | `_Float128` + strfromf128 | ― |
| **Fedora 42** | amd64 | x86_64 | LP64 | gcc 15.2 | `_Float128` + strfromf128 | `_Float64x (binary80)` |
| **Fedora 42** | arm64 | aarch64 | LP64 | gcc 15.2 | `_Float128` + strfromf128 | ― |
| **Fedora 42** | ppc64le | ppc64le | LP64 | gcc 15.2 | `_Float128` + strfromf128 | ― |
| **Fedora 42** | s390x | s390x | LP64 | gcc 15.2 | `_Float128` + strfromf128 | ― |
| **Fedora 43** | amd64 | x86_64 | LP64 | gcc 15.2 | `_Float128` + strfromf128 | `_Float64x (binary80)` |
| **Fedora 43** | arm64 | aarch64 | LP64 | gcc 15.2 | `_Float128` + strfromf128 | ― |
| **Fedora 43** | ppc64le | ppc64le | LP64 | gcc 15.2 | `_Float128` + strfromf128 | ― |
| **Fedora 43** | s390x | s390x | LP64 | gcc 15.2 | `_Float128` + strfromf128 | ― |
| **openSUSE Leap 15** | amd64 | x86_64 | LP64 | gcc 7.5 | `__float128` + libquadmath | `long double (binary80)` |
| **openSUSE Leap 15** | arm64 | aarch64 | LP64 | gcc 7.5 | `long double (binary128)` | ― |
| **openSUSE Leap 15** | s390x | s390x | LP64 | gcc 7.5 | `long double (binary128)` | ― |
| **openSUSE Leap 16** | amd64 | x86_64 | LP64 | gcc 15.2 | `_Float128` + strfromf128 | `_Float64x (binary80)` |
| **openSUSE Leap 16** | arm64 | aarch64 | LP64 | gcc 15.2 | `_Float128` + strfromf128 | ― |
| **openSUSE Leap 16** | ppc64le | ppc64le | LP64 | gcc 15.2 | `_Float128` + strfromf128 | ― |
| **openSUSE Tumbleweed** | i386 | x86_64 | ILP32 | gcc 15.2 | `_Float128` + strfromf128 | `_Float64x (binary80)` |
| **openSUSE Tumbleweed** | amd64 | x86_64 | LP64 | gcc 15.2 | `_Float128` + strfromf128 | `_Float64x (binary80)` |
| **openSUSE Tumbleweed** | arm64 | aarch64 | LP64 | gcc 15.2 | `_Float128` + strfromf128 | ― |
| **openSUSE Tumbleweed** | ppc64le | ppc64le | LP64 | gcc 15.2 | `_Float128` + strfromf128 | ― |
| **openSUSE Tumbleweed** | riscv64 | riscv64 | LP64 | gcc 15.2 | `_Float128` + strfromf128 | ― |
| **openSUSE Tumbleweed** | s390x | s390x | LP64 | gcc 15.2 | `_Float128` + strfromf128 | ― |

**Total: 77 distinct (OS, arch, toolchain) configurations** (82 raw summaries deduplicated;
the four filename-typo `ensuse-tw_*` entries are duplicates of the `opensuse-tw_*` builds on
arm64/ppc64le/riscv64/s390x and are folded into the openSUSE Tumbleweed rows above).

> **Filename hygiene note.** The release directory contains four mis-named files
> `ensuse-tw_linux-{arm64,ppc64le,riscv64,s390x}_branch.ok` ― the `op` prefix is missing.
> Their configure summaries are byte-identical to the corresponding `opensuse-tw_*` runs.
> This should be fixed in the release harness.

---

## Table 2: Unique `(binary128 config, binary80 config)` Combinations

| # | binary128 config | binary80 config | Characterization |
|---|---|---|---|
| **A** | `__float128` + libquadmath | `long double (binary80)` | x86/x86_64 + *any* of {GCC ≤ 12 on glibc, musl libc (any GCC), MinGW, Intel oneAPI}. |
| **B** | `__float128` + libquadmath | ― | non-x86 with libquadmath path: ppc64le on older-GCC glibc (Debian 12, Ubuntu 22.04, Rocky 9) and **ppc64le on musl (Alpine 3.22/3.23)**. |
| **C** | `_Float128` + strfromf128 | `_Float64x (binary80)` | x86/x86_64 on GCC ≥ 13 + glibc ≥ 2.29 (Debian 13, Ubuntu 24.04, Rocky 10, Fedora 42/43, openSUSE Leap 16, openSUSE Tumbleweed). |
| **D** | `_Float128` + strfromf128 | ― | non-x86 on GCC ≥ 13 + glibc (aarch64, ppc64le, s390x, riscv64 across Debian 13, Ubuntu 24.04, Rocky 10, Fedora 42/43, openSUSE Leap 16, openSUSE Tumbleweed). |
| **E** | `long double (binary128)` | ― | non-x86 with older GCC (≤ 12) *or* with musl, where the ABI's 128-bit `long double` is the only quad path. Covers aarch64, s390x, riscv64, mips64le across Alpine 3.19–3.23, Debian 11/12, Ubuntu 18–22, Rocky 8/9, openSUSE Leap 15. |

No sixth combination occurs in the matrix.

---

## Table 3: C++ Standard Math Exposure Summary

This table summarizes the current configure probes and header behavior for `std::abs`,
`std::sin` / `std::cos` / related scalar math overloads, and complex math. The scalar
`std::math` checks are group probes: if any listed overload is missing, MPLAPACK does not use
the `std::` overload set for that backend and instead keeps its C/libquadmath wrapper path.

| Area | Selected type / mode | Configure macro(s) | Probed functions | Header behavior |
|---|---|---|---|---|
| binary128 scalar abs | `_Float128` or `__float128` | `MPLAPACK_HAVE_STD_ABS_FLOAT128` | `std::abs(_Float128)` or `std::abs(__float128)`, matching the selected binary128 mode | If present, expose `using std::abs`; otherwise define fallback `fabsf128` / `fabsq` wrapper. |
| binary128 scalar math | `_Float128` only | `MPLAPACK_HAVE_STD_MATH_FLOAT128` | `std::sin`, `std::sinh`, `std::cos`, `std::cosh`, `std::atan2`, `std::exp`, `std::log`, `std::log10`, `std::log2`, `std::pow`, `std::sqrt`, `std::ceil`, `std::nextafter`, `std::ldexp` | If present, expose the `std::` overload set; otherwise define `*f128` wrappers. The `__float128` mode stays on libquadmath wrappers, and `long double (binary128)` uses the normal long-double overload set. |
| binary128 complex math | `_Float128` complex | `MPLAPACK_HAVE_STD_COMPLEX_FLOAT128`, `MPLAPACK_HAVE_C_COMPLEX_FLOAT128` | `std::complex<_Float128>`: `abs`, `sqrt`, `sin`, `cos`, `exp`, `log`; C complex availability: `cabsf128`, `csqrtf128`, `csinf128`, `ccosf128`, `cexpf128`, `clogf128` | Prefer `std::complex<_Float128>` math when available; otherwise use manual formulas backed by `*f128` scalar functions. The C complex probe is reported for platform visibility, but installed headers avoid depending on C `_Float128 _Complex` declarations because downstream include order and compiler mode can hide them. |
| binary80 scalar abs | `_Float64x` | `MPLAPACK_HAVE_STD_ABS_FLOAT64X` | `std::abs(_Float64x)` | If present, expose `using std::abs`; otherwise define fallback `fabsf64x` wrapper. `long double (binary80)` uses the normal long-double overload set. |
| binary80 scalar math | `_Float64x` | `MPLAPACK_HAVE_STD_MATH_FLOAT64X` | `std::sin`, `std::sinh`, `std::cos`, `std::cosh`, `std::atan2`, `std::exp`, `std::floor`, `std::log`, `std::log10`, `std::log2`, `std::pow`, `std::sqrt`, `std::nextafter`, `std::ldexp` | If present, expose the `std::` overload set; otherwise define `*f64x` wrappers. The `long double (binary80)` mode uses the normal long-double overload set. |
| binary80 complex math | `_Float64x` complex | none today | none today | The header currently provides only simple algebra helpers such as `pow2` / `pow4` for `std::complex<_Float64x>`; no complex transcendental probe is used today. |

The binary128 scalar group currently includes `ceil`, while the binary80 scalar group includes
`floor`, because those are the global helper wrappers currently exposed by the respective
headers. This is not a mathematical distinction between the formats; it is a probe boundary
chosen to match the functions that would otherwise be defined as global wrappers and could
therefore conflict with C++23 `<math.h>` / `<cmath>` declarations.

---

## Technical Findings

### 1. The binary128 code path is determined by **three orthogonal axes**, not just GCC version

The earlier version of this document said "GCC version is the primary determinant." That was
incomplete. The correct rule is:

```
_Float128 + strfromf128  ⇔  (GCC ≥ 13) ∧ (libc provides strfromf128) ∧ (x86/x86_64 OR non-x86)
__float128 + libquadmath ⇔  else, *if* the target is x86/x86_64 or ppc64le
long double (binary128)  ⇔  else, on archs where long double IS binary128 (aarch64/s390x/riscv64/mips64le)
```

The "libc provides `strfromf128`" predicate is where musl diverges sharply from glibc:
**Alpine 3.19/3.22/3.23 with GCC 13/14/15 all stay on the `__float128` + libquadmath path**,
because musl does not ship `strfromf128` (nor `strfromf64x`). This contradicts any rule of the
form "GCC ≥ 13 → new path." Consequently the full x86_64 matrix is:

| Toolchain | libc | Path |
|---|---|---|
| GCC 7–12 | glibc | A |
| GCC 13+ | glibc | **C** |
| GCC 13+ | musl (Alpine) | **A** |
| icpx 2025.x | glibc | A |
| MinGW-w64 gcc 10/13/14 | mingw-w64 CRT | A |

Any conditional compilation assuming "new GCC → `_Float128` path" will silently mis-detect on
Alpine, MinGW, and Intel oneAPI, and must also check the runtime.

### 2. binary80 is strictly x86/x86_64-only — confirmed without exception

Every non-x86 architecture across all 77 configurations produces `―` for binary80:
aarch64, ppc64le, s390x, mips64le, riscv64. This is a hard ABI/hardware constraint. Any
conditional compilation relying on binary80 must be gated by **architecture**, not compiler
version and not libc. The previous document's finding on this point stands unchanged and is
now corroborated across 15 additional platform entries.

### 3. The `long double (binary128)` path (combination E) is being **retired** by distro upgrades

On non-x86 architectures, the choice between E (`long double` as the quad path) and D
(`_Float128 + strfromf128` as the quad path) is entirely a function of GCC version and libc:

| Arch | Old-GCC path (≤ 12) | New-GCC path (≥ 13, glibc) |
|---|---|---|
| aarch64 | E | D |
| s390x | E | D |
| riscv64 | E | D |
| mips64le | E | (not yet observed on new GCC in this matrix) |
| ppc64le | **B** (note: ppc64le long double is not binary128 — it is IBM double-double or IEEE-128 depending on ABI flags, so it falls through to libquadmath) | D |

Combination E disappears entirely from the Debian 13 / Ubuntu 24.04 / Rocky 10 / Fedora 42+ /
openSUSE Leap 16 / Tumbleweed rows. MPLAPACK code that special-cases "long double IS the quad
type" will become dead on modern non-x86 Linux.

### 4. ppc64le is the **portability trap**: libquadmath present, binary80 absent, and `long double` is neither binary128 nor binary80

Combination **B** — `__float128 + libquadmath` with no binary80 — is exclusive to ppc64le on
older-toolchain builds (Debian 12, Ubuntu 22.04, Rocky 9, **Alpine 3.22/3.23**). Three traps
coexist on this architecture:

- Presence of libquadmath does *not* imply presence of binary80 (as the x86 correlation would
  suggest).
- ppc64le `long double` is either IBM double-double (`-mabi=ibmlongdouble`, still the default
  on many distros) or IEEE binary128 (`-mabi=ieeelongdouble`). Neither is safe to assume.
  Hence ppc64le never takes combination E.
- On GCC ≥ 13 + glibc (Debian 13, Ubuntu 24.04, Rocky 10, Fedora 42+), ppc64le moves cleanly to
  combination D, but Alpine's musl keeps it pinned at B regardless of GCC version.

### 5. MinGW always falls back to combination A regardless of host GCC

Host MinGW-w64 gcc 10 (Ubuntu 22.04), gcc 13 (Ubuntu 24.04), and gcc 14 (Debian 13) all resolve
to `__float128 + libquadmath` + `long double (binary80)`. The Windows target CRT does not
provide `strfromf128`/`strfromf64x`, so combination **C** is unreachable from any MinGW build
on the current matrix. This is the same pattern as musl.

### 6. Intel oneAPI 2025.3.2: now emits a well-formed summary, lands in combination A

The previous document reported that `icpx` short-circuited the `_Float128` probe. With the
current configure logic and oneAPI 2025.3.2, both the Ubuntu 22.04 and Ubuntu 24.04 `icpx`
builds emit a complete summary and resolve to `__float128 + libquadmath` + `long double
(binary80)` — i.e., combination A, the same as old-GCC glibc. The earlier note can be
retired; the current configure handling is adequate.

### 7. Rocky 9 binary80 summary: now emitted normally

The previous document carried a caveat that Rocky 9 / GCC 11.5.0 did not emit the final
`checking binary80 configuration...` summary line. The current build matrix does emit it:
Rocky 9 amd64 shows `Binary80 support: enabled` with `Type: long double (binary80)` in the
same form as every other GCC 11 x86_64 build. The caveat can be removed.

### 8. RISC-V: on musl it regresses to combination E; on glibc (GCC ≥ 13) it is combination D

- **Alpine 3.22/3.23 riscv64** (musl + GCC 14/15): combination **E** — `long double (binary128)`.
  musl does not provide `strfromf128`, so the configure falls through past `_Float128` to the
  `long double` path.
- **Debian 13 riscv64**, **Rocky 10 riscv64**, **openSUSE Tumbleweed riscv64** (glibc + GCC ≥ 14):
  combination **D** — `_Float128 + strfromf128`.

RISC-V therefore demonstrates all three non-x86 paths (B is the one it does not hit, because
its `long double` *is* binary128). It is a good regression target because the path depends on
libc in a way that x86_64 and aarch64 tests will not catch.

### 9. mips64le: configure reports CPU as `mips64` (not `mips64el`)

Debian 12 mips64le (`mips64el-linux-gnuabi64`) resolves to combination **E**. The configure
summary reports `CPU: mips64` rather than `mips64le` or `mips64el`; this is a cosmetic issue
in the CPU detection code but worth noting if downstream code keys off `$target_cpu`.

### 10. Enterprise distros will stay on the old path for their full lifecycle

- **Rocky 8** (GCC 8.5): combination A/E, no path to C/D without toolchain backports.
- **Rocky 9** (GCC 11.5): combination A/B/E.
- **openSUSE Leap 15** (GCC 7.5): combination A/E.
- **Ubuntu 18.04 / 20.04 / 22.04**: combination A/B/E.

Only Rocky 10 (GCC 14.3), Fedora 42/43 (GCC 15.2), openSUSE Leap 16 (GCC 15.2), openSUSE
Tumbleweed (GCC 15.2), Debian 13 (GCC 14.2), and Ubuntu 24.04 (GCC 13.3) reach the `_Float128`
path. MPLAPACK's libquadmath support cannot be dropped in the foreseeable future.

---

## Summary for Conditional Compilation

A correct guard expression for MPLAPACK backend selection on binary128/binary80 is:

```c
/* binary128 backend selection */
#if defined(__GLIBC__) && defined(__STDC_IEC_60559_BFP__) && (__GNUC__ >= 13) \
    && defined(__STRFROMF128_DEFINED__)  /* proxy for glibc >= 2.26/2.29 */
  /* Combination C/D: _Float128 + strfromf128 */
#elif defined(__SIZEOF_FLOAT128__)
  /* Combination A/B: __float128 + libquadmath (covers musl, MinGW, icpx, old GCC,
     and non-x86 ppc64le where long double is neither binary128 nor binary80) */
#elif defined(__aarch64__) || defined(__s390x__) || defined(__riscv) \
   || (defined(__mips64) && !defined(_MIPS_SIM_ABI32))
  /* Combination E: long double IS binary128 */
#else
  #error "no binary128 backend available"
#endif

/* binary80 backend selection */
#if !(defined(__i386__) || defined(__x86_64__))
  /* no binary80 on this architecture — hard constraint */
#elif defined(__GLIBC__) && (__GNUC__ >= 13) && defined(__STRFROMF64X_DEFINED__)
  /* Combination C: _Float64x + strfromf64x */
#else
  /* Combination A: long double (binary80) */
#endif
```

The important invariants are:

1. **Architecture first** for binary80 — there is no software fallback.
2. **libc second** for binary128 — musl blocks the `_Float128` path even on bleeding-edge GCC.
3. **GCC version third** — the textbook rule, but only applicable after the first two.
