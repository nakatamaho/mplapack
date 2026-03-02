# MPLAPACK — Multi-Precision Linear Algebra Package

MPLAPACK is a multi-precision linear algebra package based on BLAS and LAPACK,
implemented in C++. It supports a range of high-precision arithmetic libraries
including GMP, MPFR, and QD, allowing users to select the backend that best suits
their requirements. MPLAPACK is free software distributed under a 2-clause BSD-style
license, supplemental to the original LAPACK license.

# News

* 2026-03-02  MPLAPACK 2.1.0 released.
* 2022-09-12  MPLAPACK 2.0.1 released, featuring CUDA versions of Rgemm (dd) and Rsyrk (dd)
  for Volta and Ampere architectures (~450 GFlops on V100), and Windows DLLs for MinGW-w64.
* 2022-07-26  MPLAPACK 2.0.0 released. All routines (except mixed-precision) functional and
  tested across all supported precisions.

# Release Checksums

## MPLAPACK 2.1.0

| File | Algorithm | Checksum |
|:---|:---|:---|
| `mplapack-2.1.0.tar.gz` | MD5 | `d436f2fc61f6f010ab7be156c2b949e2` |
| `mplapack-2.1.0.tar.gz` | SHA-256 | `45b8a147b0d5c49f609cfe54f3f1f1a5d4ea0b77aba2b37d5072221cf50be496` |

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
* **binary128** (via glibc or libquadmath; auto-detected)
* **binary80** (xtended precision; Intel/AMD x86 only)

# Supported Platforms

| Tier | Guarantee | Platforms |
|---|---|---|
| **Tier 1** | Functionality guaranteed (build + full test suite) | macOS Intel Sonoma (amd64), Windows / MinGW-w64 (amd64), Ubuntu 22.04 (amd64), Ubuntu 24.04 (amd64) |
| **Tier 2** | Build guaranteed | Alpine Linux 3.19 (amd64, arm64), Rocky Linux 8/9 (amd64), Debian 12 (amd64, arm64, i386, mips64le), Debian 13 (amd64), Ubuntu 22.04 (arm64), Ubuntu 24.04 (arm64) |
| **Tier 3** | Patches accepted; no CI coverage | Other platforms |

# How to Build and Install

## Prerequisites

* GCC / G++ / GFortran
* Standard autotools: `autoconf`, `automake`, `libtool`
* `wget` or `curl` (to fetch the tarball)

All third-party libraries (GMP, MPFR, MPC, QD, OpenBLAS, dlfcn-win32) are bundled and
built automatically. No separate installation of these libraries is required.

## Linux (amd64 / arm64)

```sh
cd $HOME/tmp
wget https://github.com/nakatamaho/mplapack/releases/download/v2.1.0/mplapack-2.1.0.tar.gz
tar xvf mplapack-2.1.0.tar.gz
cd mplapack-2.1.0
export CXX=g++
export CC=gcc
export FC=gfortran
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

To enable `binary80` (Intel/AMD x86 only), add `--enable-binary80=yes`:

```sh
./configure \
    --prefix=$HOME/MPLAPACK \
    --enable-binary80=yes \
    --enable-gmp=yes \
    --enable-mpfr=yes \
    --enable-binary128=yes \
    --enable-qd=yes \
    --enable-dd=yes \
    --enable-double=yes \
    --enable-test=yes \
    --enable-benchmark=yes
```

## macOS (Intel; using MacPorts)

FSF GCC is required. Clang (the default macOS compiler) does neither support `_Float128` nor `__float128` 

```sh
sudo port install gcc14 coreutils git gsed
cd $HOME/tmp
wget https://github.com/nakatamaho/mplapack/releases/download/v2.1.0/mplapack-2.1.0.tar.gz
tar xvf mplapack-2.1.0.tar.gz
cd mplapack-2.1.0
export CXX=g++-mp-14
export CC=gcc-mp-14
export FC=gfortran-mp-14
./configure \
    --prefix=$HOME/MPLAPACK \
    --enable-gmp=yes \
    --enable-mpfr=yes \
    --enable-binary128=yes \
    --enable-binary80=yes \
    --enable-qd=yes \
    --enable-dd=yes \
    --enable-double=yes \
    --enable-test=yes \
    --enable-benchmark=yes
make -j$(sysctl -n hw.logicalcpu)
make install
```

## Windows (MinGW-w64 cross-compile on Ubuntu)

```sh
sudo apt-get install gcc-mingw-w64-x86-64 g++-mingw-w64-x86-64 gfortran-mingw-w64-x86-64
cd $HOME/tmp
wget https://github.com/nakatamaho/mplapack/releases/download/v2.1.0/mplapack-2.1.0.tar.gz
tar xvf mplapack-2.1.0.tar.gz
cd mplapack-2.1.0
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

## Verifying the Build

```sh
make check
```

Test results are summarized automatically by `misc/summarize_mplapack_tests.py`.

# Fable — Fortran-to-C++ Conversion Pipeline

MPLAPACK 2.1.0 includes [Fable](https://cci.lbl.gov/fable/) as a top-level standalone
component (`fable/`). BLAS/LAPACK library routines and test programs (EIG/LIN/MATGEN)
are generated from LAPACK 3.9.1 Fortran sources via automated conversion.

To regenerate library routines or test programs from source:

```sh
bash fable/go.sh          # library routines (BLAS/LAPACK C++ sources + headers + patches)
bash fable/go_testing.sh  # test programs (EIG/LIN/MATGEN C++ sources + headers + patches)
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

## MPLAPACK 2.1.0 Release Process

### Tier-S Representative Gate Matrix (Release Blockers)

| # | Tier | OS | Arch | Compiler | binary80 | binary128 | Required tasks | Date |
|---:|:---:|:---|:---|:---|:---:|:---:|:---|:---|
| 1 | S | Ubuntu 24.04 | amd64 | GCC | ✅ | ✅ | build / tests / `make check` / examples | - |
| 2 | S | Ubuntu 22.04 | amd64 | GCC | ✅ | ✅ | build / tests / `make check` / examples | - |
| 3 | S | Ubuntu 24.04 | arm64 | GCC | N/A | ✅ | build / tests / `make check` / examples | - |
| 4 | S | Debian 12 | i386 | GCC | ✅ | N/A | build / tests / `make check` / examples | - |
| 5 | S | Rocky 8/9 | amd64 | GCC | ✅ | ✅ | build / tests / `make check` / examples | - |
| 6 | S | macOS Intel Sonoma | amd64 | GCC (MacPorts) | ✅ | ✅ | build / tests / `make check` / examples | - |
| 7 | S | Windows | amd64 | GCC (MinGW-w64) | ✅ | ✅ | build / tests / `make check` / examples | - |

### Tier Policy

> **Tier 1 (release blockers):** must be green on the representative gate matrix above.
> **Tier 2 (build guarantee):** tracked but not release-blocking.
> **Tier 3 (patches accepted):** no CI coverage.

#### CPU Architecture Tiers

| Tier | Architectures | Expectation |
|:---:|:---|:---|
| 1 | amd64, arm64, i386 | build + tests + `make check` + examples |
| 2 | mips64le, riscv64 | build-only by default |

#### Compiler Tiers

| Tier | Compilers | Expectation |
|:---:|:---|:---|
| 1 | GCC (native), GCC (MinGW-w64) | Must be green |
| 2 | Clang, Intel oneAPI | Build only; binary128 N/A for Clang; binary128+binary80 N/A for oneAPI |

#### Feature Tiers (Precision)

| Tier | Feature | Primary coverage targets | Notes |
|:---:|:---|:---|:---|
| 1 | binary80 | amd64, i386, Windows (MinGW-w64) | N/A on arm64, Clang is supported, oneAPI is not |
| 1 | binary128 | amd64, arm64, macOS (GCC), Windows (MinGW-w64) | N/A on Clang and oneAPI |

## MPLAPACK 3.0.0 Release Process

| Action | Date | Status | Description |
|---|---|---|---|
| Optimized implementations as default | | | |
| Add template version | | | mockup: https://github.com/nakatamaho/mplapack-template |
| Add gmpfrxx | | | https://math.berkeley.edu/~wilken/code/gmpfrxx/ |
| Add OpenBLAS for double benchmark | | | |
| Update to LAPACK 3.12.1 | | | Patches already bundled in 2.1.0 |
| FMA for QD, DD | | | |
| Add more benchmarks (Rsyev, Rgesvd, etc.) | | | |
| Add QA program for BLAS | | | |
| Take benchmark on A100 (Rgemm, Rsyrk dd) | | | |
| Python integration | | | |
| Octave integration | | | |
| Mixed-precision routines | | | |
| lp64/ilp64/llp64/ilp32 cleanup | | | |
| Eliminate compiler warnings | | | |

## Old Release Schedules

* version 2.0.1: https://github.com/nakatamaho/mplapack/blob/master/doc/Release2.0.1.md
* version 2.0.0: https://github.com/nakatamaho/mplapack/blob/master/doc/Release2.0.0.md
* version 1.0.0: https://github.com/nakatamaho/mplapack/blob/master/doc/Release1.0.0.md

# History

* 2026/03/02  MPLAPACK 2.1.0 released. Fable standalone, automated Fortran-to-C++ pipeline,
  extended build matrix (Alpine, Rocky, Debian i386, CUDA 13.1.1), numerous bug fixes.
* 2022/09/12  MPLAPACK 2.0.1 released.
* 2022/07/26  MPLAPACK 2.0.0 released.
* 2022/06/14  MPLAPACK 2.0.0 alpha released.
* 2021/11/01  1.0.1 release. Fixed DD and QD arithmetic with Intel oneAPI.
* 2021/10/01  1.0.0 release. All real LAPACK routines available; SVD and non-symmetric
  eigenproblem solvers added.
* 2021/04/11  0.9.3 release. CentOS 7 AArch64 support.
* 2021/04/01  0.9.0 release. Renamed to mplapack.
* 2012/12/25  MPACK 0.8.0. NVIDIA C2050 support for Rgemm (double-double).
* 2010/08/20  MPACK 0.6.7. Condition number estimators added (Rgecon, Rpocon).
  License changed to 2-clause BSD.
* 2009/11/24  MPACK 0.6.0.
* 2008/07/15  mpack-0.0.1. Now configurable and installable.
* 2008/06/24  Project page created.

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
