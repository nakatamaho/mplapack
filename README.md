# MPLAPACK — Multi-Precision Linear Algebra Package

MPLAPACK is a multi-precision linear algebra package based on BLAS and LAPACK,
implemented in C++ (C++17 required). It supports a range of high-precision arithmetic
libraries including GMP, MPFR, and QD, allowing users to select the backend that best
suits their requirements. MPLAPACK is free software distributed under a 2-clause
BSD-style license, supplemental to the original LAPACK license.

# News

* 2026-08-29  MPLAPACK 3.0.0 has been released. Available from
  <https://github.com/nakatamaho/mplapack/releases/tag/v3.0.0>. This release
  completes the GMP/MPFR C++ wrapper migration to `gmpfrxx_mkII`, refreshes the
  QD backend through `libQD3`, and hardens CMake, Autotools, MinGW/Wine, macOS,
  Intel oneAPI, and release-tarball workflows. Release QA covered 11
  platform/compiler triplets: all 126,515,712 recognized `lin` tests passed;
  `eig` recorded 30 known non-harmful threshold-edge exceedances in 28 output
  records out of 116,811,306 tests. See
  [CHANGES.3.0.0.md](CHANGES.3.0.0.md) and
  [doc/Release3.0.0.md](doc/Release3.0.0.md).

* 2026-07-29  MPLAPACK 2.3.0 has been released. Available from
  <https://github.com/nakatamaho/mplapack/releases/tag/v2.3.0>. This release
  adopts one self-contained library per backend/flavor, removes the separate
  BLAS-only products, eliminates duplicate public symbols with layered basename
  shadowing, and completes the CMake build/test/install interface for the
  unified layout. See
  [CHANGES.2.3.0.md](CHANGES.2.3.0.md) and [MIGRATION.md](MIGRATION.md).

* 2026-07-21  MPLAPACK 2.2.1 has been released. Available from
  <https://github.com/nakatamaho/mplapack/releases/tag/v2.2.1>. This patch
  release keeps the LAPACK 3.12.1 base from 2.2.0 and focuses on GMP
  transcendental correctness, binary128/binary80 C++ standard compatibility,
  include-order hardening, experimental CMake support, OpenCL binary128
  benchmarks, and release engineering for the expanded Tier 1/2/3 matrix.
  See [CHANGES.2.2.1.md](CHANGES.2.2.1.md).
* 2026-05-12  MPLAPACK 2.2.0 has been released. Available from
  <https://github.com/nakatamaho/mplapack/releases>. This final release is
  based on 2.2.0-rc1: rebased from LAPACK 3.9.1 to 3.12.1, added new
  BLAS/LAPACK routines including `Cgemmtr/Rgemmtr`, `Cgedmd/Rgedmd`,
  `Cgedmdq/Rgedmdq`, `Cgelst/Rgelst`, `Clatrs3/Rlatrs3`,
  `Cgeqp3rk/Rgeqp3rk`, `Cunhr_col/Rorhr_col`,
  `Cungtsqr_row/Rorgtsqr_row`, and `Claqz0-4/Rlaqz0-4`, and improved
  robustness across binary80, binary128, DD, QD, GMP, and MPFR. See
  [CHANGES.2.2.0.md](CHANGES.2.2.0.md). Differences from
  binary64 behavior and other multiprecision-specific failure modes are
  summarized in
  [MULTIPRECISION_FAILURE_MODES_2.2.0.md](MULTIPRECISION_FAILURE_MODES_2.2.0.md).
* 2026-04-09  MPLAPACK 2.1.1 released. Patch release: GCC 15 support
  (`external/gmp` C23 fix, `{C,R}drgvx` IPA-modref workaround), arm64
  promoted to **Tier 1** (Ubuntu arm64 and macOS Apple Silicon),
  musl/Alpine ARM64 build fix, MinGW `laed3` build fix, binary128
  `std::abs(__float128)` ambiguity fix on aarch64-apple-darwin, and a
  latent DD miscompilation fix (propagating `-ffp-contract=off` to all
  DD build targets). OpenBLAS bumped to 0.3.33. ABI/source compatible
  with 2.1.0.
* 2026-03-05  MPLAPACK 2.1.0 released. C++17 now required.
  binary128/binary80 naming unified, MPFR emin/emax auto-adjustment,
  extended build matrix. LAPACK 3.9.1 Fortran sources are now mechanically
  converted to idiomatic C++ via Fable and bundled in the release tarball.
  See [CHANGES.2.1.0.md](CHANGES.2.1.0.md) and [MIGRATION.md](MIGRATION.md).
* 2022-09-12  MPLAPACK 2.0.1 released, featuring CUDA versions of Rgemm (dd) and Rsyrk (dd)
  for Volta and Ampere architectures (~450 GFlops on V100), and Windows DLLs for MinGW-w64.
* 2022-07-26  MPLAPACK 2.0.0 released. All routines (except mixed-precision) functional and
  tested across all supported precisions.

# Capabilities

* **MPBLAS:** All BLAS routines in multiple-precision arithmetic.
* **MPLAPACK:** All LAPACK routines in multiple-precision arithmetic (except mixed-precision routines):
  * Linear Equations
  * Linear Least Squares (LLS) Problems
  * Generalized Linear Least Squares (LSE and GLM) Problems
  * Standard Eigenvalue and Singular Value Problems
  * Symmetric Eigenproblems (SEP)
  * Nonsymmetric Eigenproblems (NEP)
  * Singular Value Decomposition (SVD)
  * Generalized Eigenvalue and Singular Value Problems
  * Generalized Symmetric Definite Eigenproblems (GSEP)
  * Generalized Nonsymmetric Eigenproblems (GNEP)
  * Generalized Singular Value Decomposition (GSVD)

# Supported Precision Backends

* **MPFR + MPC** https://www.mpfr.org/ / http://www.multiprecision.org/mpc/
  (arbitrary precision with IEEE-like rounding; primary backend)
* **GMP** https://gmplib.org/ (arbitrary precision)
* **double** (binary64)
* **DD, QD** https://www.davidhbailey.com/dhbsoftware/
  (DD ≈ binary128, QD ≈ binary256)
* **binary128** (IEEE 754-2018; compiler and platform support is complex —
  see [binary128_binary80_type_support.md](binary128_binary80_type_support.md))
* **binary80** (80-bit extended precision; Intel/AMD x86 only)

# Compiler Support

| Compiler | binary128 | binary80 |
|---|---|---|
| GCC | ✅ Supported | ✅ (x86/x86\_64 only) |
| Intel oneAPI (icx/icpx) | ✅ Supported | ✅ (x86/x86\_64 only) |
| Clang/LLVM | ❌ Not supported | ✅ Supported |

> **Clang users:** `binary128` is not supported. Use GCC for `binary128`.
> **GCC 15:** supported as of 2.1.1. On musl-based distros (e.g. Alpine),
> the binary128 backend stays on `__float128 + libquadmath` because musl does
> not ship `strfromf128`/`strfromf64x`.

# Supported Platforms

| Tier | Guarantee | Platforms |
|---|---|---|
| **Tier 1** | Dedicated release `buildtest` target; `make distcheck` + full test suite | macOS 15/26 (amd64/arm64), Ubuntu 24.04/26.04 (amd64/arm64), Windows / MinGW-w64 via Ubuntu 24.04/26.04 (amd64), Ubuntu 24.04/26.04 Intel oneAPI (amd64), Debian 12/13 i386 |
| **Tier 2** | Docker branch build matrix coverage via `make tier2`; no dedicated release `buildtest` target | Other Debian/Ubuntu architectures and versions, Alpine Linux 3.19–3.23, Rocky Linux 8/9/10, Fedora 42/43, openSUSE Leap 15.6/16.0, openSUSE Tumbleweed |
| **Tier 3** | Build-only Docker coverage for selected non-release-gating checks via `make tier3`; otherwise patches accepted | Ubuntu 26.04 amd64 C++ standard variants; other platforms |

Release test targets are run from `release/`. `make tier1` runs dedicated remote buildtests; `make tier2` runs the Docker branch matrix entries from `release/build-matrix.conf`; `make tier3` runs build-only non-release-gating checks.

| Tier | Make target | OS | CPU | Host | Backend |
|---|---|---|---|---|---|
| Tier 1 | `tier1-macos-arm64` | macOS 26 | arm64 | `172.27.109.40` | SSH |
| Tier 1 | `tier1-macos-amd64` | macOS 15 | amd64 | `172.27.109.97` | SSH |
| Tier 1 | `tier1-ubuntu2404-arm64` | Ubuntu 24.04 | arm64 | `172.27.109.40` | Docker/Colima |
| Tier 1 | `tier1-ubuntu2604-arm64` | Ubuntu 26.04 | arm64 | `172.27.109.40` | Docker/Colima |
| Tier 1 | `tier1-ubuntu2404-amd64` | Ubuntu 24.04 | amd64 | `172.27.109.80` | Docker |
| Tier 1 | `tier1-ubuntu2604-amd64` | Ubuntu 26.04 | amd64 | `172.27.109.80` | Docker |
| Tier 1 | `tier1-ubuntu2404-mingw64-amd64` | Windows via Ubuntu 24.04 | amd64 | `172.27.109.80` | Docker + MinGW64/Wine |
| Tier 1 | `tier1-ubuntu2604-mingw64-amd64` | Windows via Ubuntu 26.04 | amd64 | `172.27.109.80` | Docker + MinGW64/Wine |
| Tier 1 | `tier1-ubuntu2404-inteloneapi-amd64` | Ubuntu 24.04 | amd64 | `172.27.109.80` | Docker + Intel oneAPI |
| Tier 1 | `tier1-ubuntu2604-inteloneapi-amd64` | Ubuntu 26.04 | amd64 | `172.27.109.80` | Docker + Intel oneAPI |
| Tier 1 | `tier1-debian12-i386` | Debian 12 | i386 | `172.27.109.80` | Docker |
| Tier 1 | `tier1-debian13-i386` | Debian 13 | i386 | `172.27.109.80` | Docker |
| Tier 3 | `tier3-ubuntu2604-cxxstd-arm64` | Ubuntu 26.04 | arm64 | `172.27.109.40` | Docker/Colima build-only |
| Tier 3 | `tier3-ubuntu2604-cxxstd-amd64` | Ubuntu 26.04 | amd64 | `172.27.109.80` | Docker build-only |

Dedicated release buildtest scripts are in `release/`:
```
release/buildtest_tier1_macos_amd64.sh
release/buildtest_tier1_macos_arm64.sh
release/buildtest_tier1_mingw64_amd64.sh
release/buildtest_tier1_ubuntu_amd64.sh
release/buildtest_tier1_ubuntu_arm64.sh
release/buildtest_tier1_ubuntu_inteloneapi_amd64.sh
release/buildtest_tier1_debian_i386.sh
release/buildtest_tier3_ubuntu2604_cxxstd_arm64.sh
release/buildtest_tier3_ubuntu2604_cxxstd_amd64.sh
```

# How to Build and Install

## Prerequisites

* GCC / G++ / GFortran (C++17 is required; **GCC 15 supported as of 2.1.1**)
* `wget` or `curl` for release tarball builds, or `git` for Git checkout builds
* Standard autotools (`autoconf`, `automake`, `libtool`) are required only when building from a Git checkout
* `ccache` is recommended. The `misc/reconfig.*.sh` helper scripts use it by default.

All third-party libraries (GMP, MPFR, MPC, QD, OpenBLAS, dlfcn-win32) are bundled and
built automatically. No separate installation of these libraries is required. For benchmark
builds, `--with-openblas=system` uses the distribution OpenBLAS (`-lopenblas`) instead of
building the bundled OpenBLAS; pass `--with-openblas="-L/path -lopenblas"` for a custom
OpenBLAS location.

Release tarballs already include generated autotools files, so `autoreconf` is not needed.
Git checkout builds must run one of the `misc/reconfig.*.sh` scripts before `make`; these
scripts generate local build files, run `autoreconf`, and execute `./configure` with the
usual MPLAPACK options.

## Linux (amd64 / arm64)

### Release tarball

```sh
mkdir -p $HOME/tmp && cd $HOME/tmp
wget https://github.com/nakatamaho/mplapack/releases/download/v2.3.0/mplapack-2.3.0.tar.xz
tar xvf mplapack-2.3.0.tar.xz
cd mplapack-2.3.0
export CXX=g++ CC=gcc FC=gfortran
./configure \
    --prefix=$HOME/MPLAPACK \
    --enable-gmp=yes \
    --enable-mpfr=yes \
    --enable-binary128=yes \
    --enable-qd=yes \
    --enable-dd=yes \
    --enable-double=yes \
    --enable-test=yes \
    --enable-benchmark=yes
make -j$(nproc)
make install
```

### Git checkout

```sh
mkdir -p $HOME/tmp && cd $HOME/tmp
git clone https://github.com/nakatamaho/mplapack.git
cd mplapack
bash misc/reconfig.ubuntu24.04.sh
make -j$(nproc)
make install
```

For Ubuntu 22.04 or 26.04, use `misc/reconfig.ubuntu22.04.sh` or
`misc/reconfig.ubuntu26.04.sh` instead. The Linux reconfig scripts enable
`binary80` automatically on Intel/AMD x86 and omit it on non-x86 CPUs.

To also enable `binary80` manually in a tarball build (Intel/AMD x86 only), add
`--enable-binary80=yes`:

```sh
./configure \
    --prefix=$HOME/MPLAPACK \
    --enable-binary80=yes \
    --enable-binary128=yes \
    --enable-gmp=yes \
    --enable-mpfr=yes \
    --enable-qd=yes \
    --enable-dd=yes \
    --enable-double=yes \
    --enable-test=yes \
    --enable-benchmark=yes
```

## macOS (Intel or Apple Silicon; using MacPorts)

FSF GCC is required. The default Apple Clang does not support `binary128`.
On Apple Silicon (arm64), omit `--enable-binary80=yes` (binary80 is x86-only).

### Release tarball

```sh
sudo port install gcc15 coreutils git gsed
mkdir -p $HOME/tmp && cd $HOME/tmp
wget https://github.com/nakatamaho/mplapack/releases/download/v2.3.0/mplapack-2.3.0.tar.xz
tar xvf mplapack-2.3.0.tar.xz
cd mplapack-2.3.0
export CXX=g++-mp-15 CC=gcc-mp-15 FC=gfortran-mp-15
./configure \
    --prefix=$HOME/MPLAPACK \
    --enable-gmp=yes \
    --enable-mpfr=yes \
    --enable-binary128=yes \
    --enable-qd=yes \
    --enable-dd=yes \
    --enable-double=yes \
    --enable-test=yes \
    --enable-benchmark=yes
# On Intel macs, additionally pass --enable-binary80=yes
make -j$(sysctl -n hw.logicalcpu)
make install
```

### Git checkout

```sh
sudo port install gcc15 coreutils git gsed ccache autoconf automake libtool
mkdir -p $HOME/tmp && cd $HOME/tmp
git clone https://github.com/nakatamaho/mplapack.git
cd mplapack
bash misc/reconfig.macOS.sh
make -j$(sysctl -n hw.logicalcpu)
make install
```

## Windows (MinGW-w64 cross-compile on Ubuntu)

### Release tarball

```sh
sudo apt-get install gcc-mingw-w64-x86-64 g++-mingw-w64-x86-64 gfortran-mingw-w64-x86-64
mkdir -p $HOME/tmp && cd $HOME/tmp
wget https://github.com/nakatamaho/mplapack/releases/download/v2.3.0/mplapack-2.3.0.tar.xz
tar xvf mplapack-2.3.0.tar.xz
cd mplapack-2.3.0
export CXX=x86_64-w64-mingw32-g++
export CC=x86_64-w64-mingw32-gcc
export FC=x86_64-w64-mingw32-gfortran
./configure \
    --host=x86_64-w64-mingw32 \
    --prefix=$HOME/MPLAPACK \
    --enable-gmp=yes \
    --enable-mpfr=yes \
    --enable-binary128=yes \
    --enable-binary80=yes \
    --enable-qd=yes \
    --enable-dd=yes \
    --enable-double=yes \
    --enable-test=yes
make -j$(nproc)
make install
```

### Git checkout

```sh
sudo apt-get install git autoconf automake libtool ccache \
    gcc-mingw-w64-x86-64 g++-mingw-w64-x86-64 gfortran-mingw-w64-x86-64
mkdir -p $HOME/tmp && cd $HOME/tmp
git clone https://github.com/nakatamaho/mplapack.git
cd mplapack
bash misc/reconfig.ubuntu24.04.mingw64.sh
make -j$(nproc)
make install
```

## Verifying the Build

```sh
make check
```

Test results are summarized automatically by `misc/summarize_mplapack_tests.py`.

For an MPFR installation, the shared libraries record their MPC, MPFR, and
GMP runtime dependencies.  The public MPLAPACK MPFR headers include the
gmpfrxx_mkII MPC/MPFR/GMP headers, so these are expressed as `Requires:` in
the MPFR pkg-config files.  The optimized library additionally exposes its
OpenMP link flag through `Libs.private:`.  Both ordinary and static consumers
can therefore use the package metadata directly:

```sh
pkg-config --cflags --libs mplapack_mpfr
pkg-config --cflags --libs --static mplapack_mpfr
misc/check_mpfr_pkgconfig.sh shared
misc/check_mpfr_pkgconfig.sh static
```

On ELF systems, `misc/check_mpfr_elf_dependencies.sh` checks the installed
MPFR library's `DT_NEEDED` entries and runtime relocations.  These checks must
be run with the temporary installation's library directories in the loader
path; no consumer-side preloading or manually added MPFR dependency flags are
required.

On macOS, QD/DD shared-load QA uses a small native `dlopen()` probe. During
`make check`, only that probe carries the bundled QD build directory as an
RPATH; it neither links nor preloads QD. After `make install`, run
`make installcheck` to load the installed reference and optimized libraries
with a separate probe that has no build-directory RPATH or DYLD library
overrides. `make distcheck` runs both stages automatically. CMake provides
corresponding build-tree and installed-library tests. These checks do not
add build-directory RPATHs to installed MPLAPACK libraries.

## Experimental CMake Build

MPLAPACK also has an experimental CMake build. The autotools build above remains
the primary release path. Unlike autotools, the CMake build does not download or
build bundled third-party libraries; GMP, MPFR, MPC, QD, BLAS, and LAPACK are
found from the system when the corresponding backends or tests are enabled. See
[README.cmake.md](README.cmake.md) for the complete option list.

For a dependency-light build that verifies the CMake path itself, build only the
double and binary128 backends:

```sh
cmake -S . -B build-cmake \
    -DCMAKE_BUILD_TYPE=Release \
    -DMPLAPACK_CXX_STANDARD=23 \
    -DMPLAPACK_CXX_EXTENSIONS=OFF \
    -DMPLAPACK_ENABLE_GMP=OFF \
    -DMPLAPACK_ENABLE_MPFR=OFF \
    -DMPLAPACK_ENABLE_QD=OFF \
    -DMPLAPACK_ENABLE_DD=OFF \
    -DMPLAPACK_ENABLE_DOUBLE=ON \
    -DMPLAPACK_ENABLE_BINARY128=ON \
    -DMPLAPACK_ENABLE_BINARY80=OFF
cmake --build build-cmake -j
```

On x86/x86_64 systems, add `-DMPLAPACK_ENABLE_BINARY80=ON` to include the
binary80 backend. To run the CMake test targets, enable tests and use CTest:

```sh
cmake -S . -B build-cmake \
    -DCMAKE_BUILD_TYPE=Release \
    -DMPLAPACK_BUILD_TESTS=ON
cmake --build build-cmake -j
ctest --test-dir build-cmake --output-on-failure
```

Install from a CMake build tree with:

```sh
cmake --install build-cmake --prefix $HOME/MPLAPACK
```

# Fable — Fortran-to-C++ Conversion Pipeline

> **Note:** The Fable conversion pipeline is **not included in the release tarball**.
> It is available only via Git clone. Regenerating sources also requires expanding the
> bundled LAPACK 3.9.1 source under `external/lapack/` before running the scripts.

`fable/` is a top-level standalone component providing automated Fortran-to-C++ conversion
of LAPACK 3.9.1 sources via [Fable](https://cci.lbl.gov/fable/) and FEM (Fortran Emulator).

```sh
# Step 1: clone the repository
git clone https://github.com/nakatamaho/mplapack
cd mplapack

# Step 2: expand bundled LAPACK sources
cd external/lapack
tar xvf lapack-3.9.1.tar.gz
cd ../..

# Step 3: run the conversion pipeline
bash fable/go.sh          # library routines (BLAS/LAPACK C++ sources + headers + patches)
bash fable/go_testing.sh  # test programs (EIG/LIN/MATGEN)
```

# MPLAPACK Test Results

* https://github.com/nakatamaho/mplapack/tree/master/mplapack/test/lin/results
* https://github.com/nakatamaho/mplapack/tree/master/mplapack/test/eig/results

# MPLAPACK Benchmark Results

* https://github.com/nakatamaho/mplapack/tree/master/benchmark/results/2022

# Manual

* https://arxiv.org/abs/2109.13406v2
* https://raw.githubusercontent.com/nakatamaho/mplapack/master/doc/manual/manual.pdf

```bibtex
@misc{2109.13406v2,
  Author = {Maho Nakata},
  Title  = {MPLAPACK version 2.0.1 user manual},
  Year   = {2022},
  Eprint = {arXiv:2109.13406v2},
}
```

# Movies

* https://www.youtube.com/watch?v=M76wHwckNPU (created by Ge Baolai)

# Slides

* https://github.com/nakatamaho/mplapack/blob/master/doc/presentation/2023-06-01%20CMSI%E6%95%99%E8%82%B2%E8%A8%88%E7%AE%97%E7%A7%91%E5%AD%A6%E6%8A%80%E8%A1%93%E7%89%B9%E8%AB%96A%20%E7%AC%AC7%E5%9B%9E%20%E4%B8%AD%E7%94%B0%E7%9C%9F%E7%A7%80.pdf (in Japanese)
* https://github.com/nakatamaho/mplapack/blob/v2.0/doc/presentation/20211128_%E7%B2%BE%E5%BA%A6%E4%BF%9D%E8%A8%98meeting.pdf (in Japanese)

# MPLAPACK Release Process

## MPLAPACK 2.3.0 Release Process

2.3.0 is an ABI and library-layout release on top of 2.2.1. Highlights:

- **Unified self-contained libraries.** Each backend/flavor now installs one
  self-contained `libmplapack_*` library containing both MPBLAS and MPLAPACK
  symbols. The separate `libmpblas_*` products were removed.
- **Layered symbol ownership.** Accelerator sources override optimized sources,
  which override reference sources, and duplicate public symbols are rejected
  by archive checks on ELF, Mach-O, and MinGW/PE-COFF targets.
- **Autotools/CMake parity.** Both build systems generate the unified
  reference, optimized, and accelerator flavor libraries, with per-flavor
  pkg-config and CMake targets.
- **Accelerator flavors.** DD CUDA and binary128 OpenCL optimized libraries are
  self-contained and follow the same symbol-shadowing rule as CPU flavors.
- **Release validation.** The 2.3.0 lin/eig result set covers 11
  platform/compiler triplets. `lin` passed all 126,515,712 recognized tests;
  `eig` recorded 8 failures out of 116,310,045 recognized tests.

See [CHANGES.2.3.0.md](CHANGES.2.3.0.md) and [MIGRATION.md](MIGRATION.md) for
the full change summary and migration notes.

## MPLAPACK 2.2.1 Release Process

2.2.1 is a patch release on top of 2.2.0. No public BLAS/LAPACK routine set
is intentionally changed. Highlights:

- **GMP transcendental hardening.** Native GMP-backed `log2`, `log1p`,
  `log`, `exp`, `expm1`, `sin`/`cos`, `atan`, `atan2`, and `pow` helpers
  were added and tested against MPFR oracles.
- **C++ standard selection.** `--with-cxx-standard` now supports C++17,
  C++20, C++23, and C++26, including GNU and strict modes.
- **binary128/binary80 math probes.** `std::abs`, scalar math, complex math,
  and C-complex availability are now configure/CMake probe results rather
  than language-mode assumptions.
- **Header include-order hardening.** Installed binary128/binary80 BLAS and
  LAPACK headers include `mplapack_config.h` before standard headers where
  feature-test macros matter, and include `<complex>` where declarations use
  `std::complex`.
- **Experimental CMake build.** CMake can build/install MPLAPACK, package
  sources/binaries with CPack, and run selected MPBLAS/MPLAPACK tests and
  benchmarks. Autotools remains the primary release path.
- **OpenCL binary128 benchmarks.** OpenCL binary128 RGEMM support and
  benchmark targets were added and gated on a usable OpenCL setup.
- **Release engineering.** Remote Tier 1/2/3 execution, dist tarball caching,
  per-host job serialization, ccache reporting, source snapshot handling, and
  cleanup were tightened for the expanded release matrix.

See [CHANGES.2.2.1.md](CHANGES.2.2.1.md) for the full change summary.

## MPLAPACK 2.1.1 Release Process

2.1.1 is a patch release. No API/ABI changes vs. 2.1.0. Highlights:

- **GCC 15 support.** `external/gmp` patched for the C23 default
  (`void g(){}` → properly typed prototype); patch is gated on
  `GCC >= 15` at configure time. `{C,R}drgvx` test drivers carry an
  `__attribute__((optimize("O1")))` workaround for a g++-15 IPA-modref
  miscompilation.
- **arm64 promoted to Tier 1**: Ubuntu arm64 and macOS Apple Silicon.
- **musl/Alpine ARM64** build fix (`M_PIl` fallback).
- **MinGW** `lapack/laed3` build fix (`LAMC3` prototype, `max` macro guard).
- **binary128** `sign()` ambiguity on aarch64-apple-darwin / MacPorts GCC 15
  fixed (unqualified `abs()`).
- **DD miscompilation hardening**: `-ffp-contract=off` (and `-fp-model strict`
  for `icpx`) is now propagated to *all* DD build targets — reference and
  optimized BLAS/LAPACK, tests, and benchmarks. CUDA DD targets are
  intentionally excluded (`nvcc` rejects the flag).
- **`configure` fix**: a stale `CXXFLAGS="$SAVE_CXXFLAGS"` line was
  silently clearing user-supplied `CXXFLAGS` (e.g. `-fsanitize=address`).
- **`std::abs(__float128)` probe** tightened to require an exact-overload
  match instead of trusting `__SIZEOF_FLOAT128__`.
- **dd/qd test comparisons** switched from `diff` to `misc/num_diff.py`
  (relative tolerance `1e-30` for dd, `1e-60` for qd) to absorb 1-bit
  FMA-induced rounding noise.
- **`make distcheck`** fixes: `external/openblas/Makefile` moved into the
  benchmark-conditional `AC_CONFIG_FILES` block; dd/qd test scripts use
  `${srcdir}` for VPATH builds; missing `misc/` scripts added to `EXTRA_DIST`.
- **OpenBLAS** updated to 0.3.33 (also fixes the macOS arm64 build).
- **`configure` build summary**: now reports compiler commands and versions,
  target CPU and integer model, OpenMP runtime, enabled backends,
  binary128/binary80 type / I/O / math / literal-suffix / interop, and
  `std::abs(__float128)` availability.
- **Release build matrix** broadened across 77 (OS, arch, toolchain)
  configurations: Alpine 3.19–3.23, Debian 11–13, Ubuntu 18.04–24.04,
  Rocky 8–10, Fedora 42/43, openSUSE Leap 15/16, openSUSE Tumbleweed,
  on x86_64, i386/i686, aarch64, ppc64le, s390x, riscv64, mips64le,
  plus Intel oneAPI (`icpx` 2025.3.2) and MinGW-w64 cross targets.
- **Enterprise lifecycle note**: Rocky 8/9, openSUSE Leap 15, and
  Ubuntu 18/20/22 stay on the `libquadmath` path for their full support
  window — `libquadmath` support cannot be dropped yet.

## MPLAPACK 2.1.0 Release Process

### Tier-S Representative Gate Matrix (Release Blockers)

Tier 1 platforms run the full pipeline including `make distcheck`. Tier 2 platforms run the Docker branch matrix from `release/build-matrix.conf`.

| # | Tier | OS | Arch | Compiler | binary80 | binary128 | Required tasks | Date |
|---:|:---:|:---|:---|:---|:---:|:---:|:---|:---|
| 1 | 1 | macOS Intel Sonoma | amd64 | GCC (MacPorts) | ✅ | ✅ | `make distcheck` + examples | - |
| 2 | 1 | macOS Apple Silicon | arm64 | GCC (MacPorts) | N/A | ✅ | `make distcheck` + examples | 2.1.1 |
| 3 | 1 | Windows | amd64 | GCC (MinGW-w64) | ✅ | ✅ | `make distcheck` + examples | - |
| 4 | 1 | Ubuntu 22.04 | amd64 | GCC | ✅ | ✅ | `make distcheck` + examples | - |
| 5 | 1 | Ubuntu 24.04 | amd64 | GCC | ✅ | ✅ | `make distcheck` + examples | - |
| 6 | 1 | Ubuntu 24.04 | arm64 | GCC | N/A | ✅ | `make distcheck` + examples | 2.1.1 |
| 7 | 2 | Debian 12 | i386 | GCC | ✅ | N/A | build only | - |
| 8 | 2 | Rocky 8/9/10 | amd64 | GCC | ✅ | ✅ | build only | - |
| 9 | 2 | Alpine 3.19–3.23 | amd64, arm64, riscv64 | GCC | ✅ (x86) | ✅ | build only | - |
| 10 | 2 | Debian 12/13 | ppc64le, s390x, riscv64, mips64le | GCC | N/A | ✅ | build only | - |
| 11 | 2 | Debian 13 | amd64 | GCC | ✅ | ✅ | build only | - |

### Tier Policy

> **Tier 1 (release blockers):** `make distcheck` must pass on all Tier 1 platforms.
> **Tier 2 (matrix coverage):** Docker branch build matrix coverage; not release-blocking.
> **Tier 3 (patches accepted):** no CI coverage.

#### CPU Architecture Tiers

| Tier | Architectures | Expectation |
|:---:|:---|:---|
| 1 | amd64 (macOS, Windows, Ubuntu); arm64 (macOS, Ubuntu) | `make distcheck` + examples |
| 2 | i386, ppc64le, s390x, riscv64, mips64le | build only |
| 3 | others | build-only best-effort |

#### Compiler Tiers

| Tier | Compilers | Expectation |
|:---:|:---|:---|
| 1 | GCC 11–15 (native), GCC (MinGW-w64) | Must be green |
| 2 | Clang | Build only; binary128 N/A |
| — | Intel oneAPI | binary128 and binary80 broken (2024+); oneAPI 2023 worked but no longer readily available; https://github.com/nakatamaho/mplapack/issues/77 |

#### Feature Tiers (Precision)

| Tier | Feature | Primary coverage targets | Notes |
|:---:|:---|:---|:---|
| 1 | binary80 | amd64, i386, Windows (MinGW-w64) | N/A on arm64; Clang supported; oneAPI broken |
| 1 | binary128 | amd64, arm64, macOS (GCC), Windows (MinGW-w64) | N/A on Clang; oneAPI broken |

## MPLAPACK 3.0.0 Release Process

3.0.0 completes the GMP/MPFR C++ wrapper migration and the QD dependency
refresh, with release validation across Linux, macOS, MinGW/Wine, Intel oneAPI,
CMake, and source-tarball workflows.

The release uses a clean source snapshot generated from the release baseline and
the accompanying checksum files.  See [CHANGES.3.0.0.md](CHANGES.3.0.0.md)
for the release notes and [doc/Release3.0.0.md](doc/Release3.0.0.md) for the
QA and publication checklist.

## Old Release Schedules

* version 2.0.1: https://github.com/nakatamaho/mplapack/blob/master/doc/Release2.0.1.md
* version 2.0.0: https://github.com/nakatamaho/mplapack/blob/master/doc/Release2.0.0.md
* version 1.0.0: https://github.com/nakatamaho/mplapack/blob/master/doc/Release1.0.0.md

# History

* 2026-07-29  MPLAPACK 2.3.0 released. Unified self-contained libraries,
  separate MPBLAS products removed, ABI advanced to 2, Autotools/CMake layout
  parity completed, and release validation broadened across 11 triplets.
* 2026-07-21  MPLAPACK 2.2.1 released. Patch release: GMP transcendental
  hardening, C++17/20/23/26 selection, binary128/binary80 math probes,
  include-order fixes, experimental CMake support, OpenCL binary128
  benchmarks, and release engineering updates.
* 2026-05-12  MPLAPACK 2.2.0 released. LAPACK reference rebased from 3.9.1
  to 3.12.1; new BLAS/LAPACK routines and arithmetic-params infrastructure
  added.
* 2026-04-09  MPLAPACK 2.1.1 released. GCC 15 support, arm64 promoted to Tier 1
  (Ubuntu arm64, macOS Apple Silicon), DD `-ffp-contract=off` propagation fix,
  binary128 / MinGW / musl build fixes, OpenBLAS 0.3.33. ABI compatible with 2.1.0.
* 2026-03-05  MPLAPACK 2.1.0 released. binary128/binary80 naming unified, MPFR emin/emax
  auto-adjustment, extended build matrix (Alpine, Rocky, Debian i386, CUDA 13.1.1).
  LAPACK 3.9.1 Fortran sources mechanically converted to idiomatic C++ via Fable and
  bundled in the release tarball.
* 2022-09-12  MPLAPACK 2.0.1 released.
* 2022-07-26  MPLAPACK 2.0.0 released.
* 2022-06-14  MPLAPACK 2.0.0 alpha released.
* 2021-11-01  1.0.1 release. Fixed DD and QD arithmetic with Intel oneAPI.
* 2021-10-01  1.0.0 release. All real LAPACK routines available; SVD and non-symmetric
  eigenproblem solvers added.
* 2021-04-11  0.9.3 release. CentOS 7 AArch64 support.
* 2021-04-01  0.9.0 release. Renamed to mplapack.
* 2012-12-25  MPACK 0.8.0. NVIDIA C2050 support for Rgemm (double-double).
* 2010-08-20  MPACK 0.6.7. Condition number estimators added (Rgecon, Rpocon).
  License changed to 2-clause BSD.
* 2009-11-24  MPACK 0.6.0.
* 2008-07-15  mpack-0.0.1. Now configurable and installable.
* 2008-06-24  Project page created.

# Old Page

http://mplapack.sourceforge.net/

# Acknowledgment

This work has been supported by:
The Special Postdoctoral Researchers' Program of RIKEN (2008, 2009),
Grant-in-Aid for Scientific Research (B) 21300017 from the Japan Society for the Promotion
of Science (2009–2011), Microsoft Research CORE6 (2010), JSPS KAKENHI Grant no. 18H03206,
and TIS inc.

M.N. would like to thank Dr. Imamura Toshiyuki, Dr. Nakasato Naohito, Dr. Fujisawa Katsuki,
Dr. Kouya Tomonori, Dr. Takahashi Daisuke, Dr. Goto Kazushige, Dr. Himeno Ryutaro,
Dr. Hishimuna Toshiaki, Dr. Katagiri Takahiro, Dr. Ogita Takeshi, Dr. Kashiwagi Masaaki,
Dr. Yuasa Fukuko, Dr. Ishikawa Tadashi, Dr. Geshi Masaaki, and Mr. Minato Yuichiro
for warm encouragement.

# Citation

```bibtex
@misc{2109.13406v2,
  Author = {Maho Nakata},
  Title  = {MPLAPACK version 2.0.1 user manual},
  Year   = {2022},
  Eprint = {arXiv:2109.13406v2},
}
```

# Contact

NAKATA Maho <maho.nakata@gmail.com> <maho@riken.jp>
