# MPLAPACK 3.0.0 gmpfrxx_mkII Migration: P1 Report

## State

- Branch: `topic/gmpfrxx_mkII_migration`
- Starting approved P0: `0151f4a80964e1b59295aa28a9fe4309791f3f1f`
- Original migration base: `b875e74d4b927282c907c3a29e6cadda48a7d57b`
- Bootstrap: `gmpfrxx_mkII` 1.0.1
- Archive: `gmpfrxx_mkII.1.0.1.tar.xz`
- SHA-256: `c0816b3538b6b77009f714bb391cebe11abb2fdb69e07aa3bb305ff822764afb`
- Provenance: locked P0 archive reproduced from upstream bootstrap release; the committed archive was verified against `LOCK.json` and `external/distfiles.sha256`.

## Implementation

Autotools adds the established external subdirectory and two-prefix CMake build. Internal dependencies are the staged `external/i/GMP`, `external/i/MPFR`, and `external/i/MPC` prefixes. The internal package installs to `external/i/GMPFRXX_MKII`; the final-prefix test installation was `/tmp/mplapack-p1-final`.

The exact upstream options used were `GMPFRXX_MKII_DEPS_AUTO_FETCH=OFF`, `GMPFRXX_MKII_BUILD_EXAMPLES=OFF`, `GMPFRXX_MKII_BUILD_BENCHMARKS=OFF`, and `GMPFRXX_MKII_COMPONENTS=GMP,MPFR,MPC`. No patches were applied. MPLAPACK CMake only performs `find_package(gmpfrxx_mkII CONFIG REQUIRED)` when explicitly enabled and does not link any backend to it.

Internal and final manifests contain headers, `gmpfrxx_mkIIConfig.cmake`, and exported target files for `gmpxx_mkII::`, `mpfrxx_mkII::`, `mpcxx_mkII::`, and `gmpfrxx_mkII::`. The staged and independent-prefix CMake commands were the two invocations recorded by `gate-P1.sh`; both passed. The external package's real upstream CMake build and install completed. A system-prefix consumer was represented by the actual CMake package/target discovery path; no MPLAPACK backend was linked to the package in P1.

## Acceptance

Exact gate command:

```sh
bash docs/migration/gmpfrxx_mkII/tools/gate-P1.sh
```

Result: `P1 gate: PASS`.

The gate verified the archive digest, rejected a deliberately corrupted copy, removed the damaged copy, verified both installations and exported targets, confirmed no production source/public backend header includes gmpfrxx_mkII, confirmed an empty `patches/` directory, and configured MPLAPACK against both staged and independent prefixes. The committed P0 baseline remains unchanged and its frozen `baseline.json` is present. No P1 numerical baseline rule was changed. The accepted P0 pkg-config metadata gap was treated as pre-existing and was not modified.

Autotools regeneration used `./gen_configure.sh`; the generated external `Makefile.in` was produced without unrelated source churn. No P2 or interop work was started.

## Deviations

The P1 system-prefix consumer was validated through the installed package's real exported-target CMake discovery path; MPLAPACK backend targets were intentionally not linked, as required by P1. No other deviations or blockers.
