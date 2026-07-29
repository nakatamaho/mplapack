# MPLAPACK 2.3.0 - Changes since 2.2.1

Release date: 2026-07-29 (base commit: `afb9637445f980b6ab483f826c9244bac67afb49`)

MPLAPACK 2.3.0 is the unified-library release.  The public library layout now
builds one self-contained MPLAPACK library per flavor, containing both MPBLAS
and MPLAPACK symbols.  This is an intentional ABI and link-interface break from
the 2.2.x layout.

This summary is based on 23 commits in `v2.2.1..afb9637445f9`.

## Unified library layout

- Every backend/flavor is now self-contained:
  `libmplapack_<backend>`, `libmplapack_<backend>_opt`,
  `libmplapack_dd_opt_cuda`, and
  `libmplapack_binary128_opt_opencl`.
- Consumers link exactly one unified flavor library.  The separate
  `libmpblas_*` products and the aggregate `mplapack.pc` file were removed.
- Per-flavor pkg-config files and CMake export targets describe the exact
  library and precision dependencies for each flavor.
- Duplicate public symbols are eliminated by a deterministic basename
  shadowing rule: accelerator sources override optimized sources, and
  optimized sources override reference sources.
- Archive duplicate-symbol checks are part of the build/test path.  The checker
  now handles ELF, Mach-O weak external COMDAT/header-inline symbols, and
  MinGW/PE-COFF archives using a relocatable whole-archive link check.
- The optimized LAPACK layer is represented explicitly under
  `mplapack/optimized/<backend>`, even where it is currently empty, so the
  reference, optimized, and accelerator layers stay visible in the source tree.

## Version and ABI

- Package version: 2.3.0.
- The incompatible library-layout change advances libtool version-info from
  `1:0:0` to `2:0:0`.
- CMake shared libraries use version 2.0.0 and SONAME 2.
- Existing applications that linked `libmpblas_*` and `libmplapack_*`
  separately must update link lines to use one unified `libmplapack_*` flavor.

## Autotools and CMake

- Autotools and CMake now both build the unified reference, optimized, and
  accelerator flavor libraries.
- CMake provides build-tree and installed targets as
  `mplapack::mplapack_<flavor>`.
- `MPLAPACK_ENABLE_OPT` controls optimized CPU flavors.
- `MPLAPACK_ENABLE_CUDA` enables the self-contained DD CUDA flavor.
- `MPLAPACK_ENABLE_OPENCL` enables the self-contained binary128 OpenCL flavor.
- `MPLAPACK_CUDA_ARCHITECTURES` selects CUDA device architectures.
- Compare-test generation was updated for unified libraries, and stale
  `libmpblas` references are rejected.
- The old compare Makefile generator was removed.
- The source-manifest drift check compares generated source manifests against
  the reference tree and is registered in both autotools and CTest paths.
- DD compilation continues to preserve `-ffp-contract=off` on relevant
  translation units.

## Accelerator flavors

- `libmplapack_dd_opt_cuda` is now self-contained and links without separate
  MPBLAS/MPLAPACK reference libraries.
- `libmplapack_binary128_opt_opencl` is now self-contained and follows the same
  accelerator-over-optimized-over-reference source shadowing rule.
- The OpenCL `Rgemm` source was renamed to match the public basename used by
  the duplicate-symbol filter.
- CUDA builds include a workaround for the glibc C23 `rsqrt` namespace issue
  seen with newer CUDA toolchains.

## Test harness and numerical notes

- The lin/eig test harnesses are generated for optimized backends in both build
  systems.
- Standard lin/eig test registration uses the non-optimized input set, so
  release tarballs do not require `se2_opt.in` for the standard eig harness.
- The reference-double `se2` run remains at threshold 50 and passes.
- The optimized-double `se2` case was reproduced at test 36, matrix type 15,
  `N=20`, seed `60812,42731,63787,906`: ratio 57.2705 with
  `OMP_NUM_THREADS=1` and 15.2737 with the default thread count on GCC 15.2.0.
  This is treated as bounded summation-order divergence in the optimized BLAS
  path, so optimized `se2` input uses a documented threshold of 100; kernels are
  unchanged.

## Release and QA tooling

- The MPLAPACK test-result summarizer now accepts multiple result roots, prints
  triplet information, and can emit grouped text, CSV, and JSON views.
- Tier2 and tarball release logs use explicit namespaces matching the tier1
  style.
- Debian Dockerfiles use noninteractive apt/tzdata setup.
- openSUSE Leap build dependencies were refreshed, with ppc64le omitted from
  the Leap 15 matrix because binary128 support is unavailable there.
- macOS and MinGW duplicate-symbol checks were hardened for their native object
  and archive formats.
- README, README.cmake, migration notes, examples, benchmarks, packaging files,
  and repo-local agent documentation were updated for the unified layout.

## 2.3.0 test result summary

The 2.3.0 result trees were summarized from:

- `mplapack/test/eig/results/2.3.0/`
- `mplapack/test/lin/results/2.3.0/`

The result set covers 11 platform/compiler triplets:

- `Apple_M4_macos26_gcc-15_2_0`
- `Apple_M4_ubuntu24_04_gcc-13_3_0`
- `Apple_M4_ubuntu26_04_gcc-15_2_0`
- `Core_i7-6920HQ_macos15_gcc-15_2_0`
- `Ryzen_Threadripper_3970X_debian12_gcc-12_2_0`
- `Ryzen_Threadripper_3970X_debian13_gcc-14_2_0`
- `Ryzen_Threadripper_3970X_ubuntu24_04_gcc-13_3_0`
- `Ryzen_Threadripper_3970X_ubuntu24_04_icx-2026_1_0`
- `Ryzen_Threadripper_3970X_ubuntu26_04_gcc-15_2_0`
- `Ryzen_Threadripper_3970X_ubuntu26_04_icx-2026_0_0`
- `Ryzen_Threadripper_3970X_windows_gcc-13_x_x`

Aggregate result counts:

| Category | Triplets | Groups | `.out` files | Recognized tests | Failed tests | Failed records |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| eig | 11 | 96 | 4,032 | 116,310,045 | 8 | 6 |
| lin | 11 | 96 | 384 | 126,515,712 | 0 | 0 |

`lin` passed across all summarized triplets.  `eig` had 8 failed tests out of
116,310,045 recognized tests, an overall failure rate of about 0.0000069%.

`eig` failure records:

| Triplet | Precision | File | Suite | Failed/Total | Fail % |
| --- | --- | --- | --- | ---: | ---: |
| `Ryzen_Threadripper_3970X_debian12_gcc-12_2_0` | double | `double/Rse2.out` | DST | 1/4,440 | 0.0225% |
| `Ryzen_Threadripper_3970X_debian13_gcc-14_2_0` | double | `double/Rse2.out` | DST | 1/4,440 | 0.0225% |
| `Core_i7-6920HQ_macos15_gcc-15_2_0` | binary80 | `binary80/Rsep.out` | DST drivers | 2/13,464 | 0.0149% |
| `Ryzen_Threadripper_3970X_windows_gcc-13_x_x` | binary80 | `binary80/Rsep.out` | DST drivers | 2/13,464 | 0.0149% |
| `Ryzen_Threadripper_3970X_debian12_gcc-12_2_0` | double | `double/Rsvd.out` | DBD | 1/10,260 | 0.0097% |
| `Ryzen_Threadripper_3970X_debian13_gcc-14_2_0` | double | `double/Rsvd.out` | DBD | 1/10,260 | 0.0097% |

Audit note: most eig groups had 32 recognized result files.  The following
Debian groups had 31 recognized result files in the summarizer output:
Debian 12 double, Debian 12 gmp, Debian 13 dd, Debian 13 double, and
Debian 13 gmp.  All lin groups had recognized result files and no failures.

The failure summary can be reproduced from `mplapack/test` with:

```sh
python3 ../../misc/summarize_mplapack_tests.py eig eig/results/2.3.0/* --only-fail
python3 ../../misc/summarize_mplapack_tests.py lin lin/results/2.3.0/* --only-fail
```
