# Announcement draft: MPLAPACK 2.4.0

Subject: MPLAPACK 2.4.0 maintenance release

MPLAPACK 2.4.0, a maintenance release of the MPLAPACK 2.x series, is now
available. MPLAPACK provides BLAS and LAPACK in multiple-precision arithmetic
using GMP, MPFR/MPC, DD, QD, binary80, binary128, and double backends.

Based on LAPACK 3.12.1, this release builds on MPLAPACK 2.3.0 and retains its
self-contained, per-flavor library layout and SONAME 2. It is intended for
users maintaining applications on the 2.x backend interface, rather than
moving to the 3.x series.

This is a small maintenance release. The only user-visible behavior change:

- MPLAPACK no longer sets the MPFR exponent range. Unless specified with
  `MPLAPACK_MPFR_EMIN` / `MPLAPACK_MPFR_EMAX`, the MPFR default range is used.
  Applications relying on the previous automatic precision-dependent limits
  should set their intended limits explicitly.

All other changes are QA and release-engineering matters: dependency linkage
and pkg-config/CMake metadata, installation order, LIN/EIG test-driver
ownership and backend initialization in static test links, and release
validation.

The release QA run completed all 12 Tier 1 jobs and all 4 Tier 2 jobs with
OK status. Coverage includes macOS, Ubuntu amd64/arm64, Windows via
MinGW-w64/Wine, Intel oneAPI, Debian i386, and C++17/20/23/26 buildability.

Numerical QA: LIN reports no failures in 126,515,712 recognized tests.
EIG reports 8 failures in 116,310,045 recognized tests (6 records): Debian i386
double Rse2/Rsvd and binary80 Rsep on macOS 15 and MinGW. These are the same
exceptions reported for 2.3.0. Full records are in the release notes.

Downloads:
https://github.com/nakatamaho/mplapack/releases/tag/v2.4.0

Changes since the latest 2.3.x release (2.3.0):
See the attached `CHANGES.2.4.0.md`.

Release source: `00e0b659d01d073a36a58fa5851fae067fdd3a10`

SHA-256 (`mplapack-2.4.0.tar.xz`):

```text
bfdd8c64498a1b9472b98d3f40bcf94a51fb3e54f2921996cf1eb2c85c6093df
```

MPLAPACK is free software under a 2-clause BSD-style license, supplemental to
the original LAPACK license.

Project: https://github.com/nakatamaho/mplapack
