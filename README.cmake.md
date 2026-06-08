# Building MPLAPACK with CMake

This is a modern CMake build for MPLAPACK that lives alongside the existing
autotools build. It produces, for each enabled precision backend, a pair of
libraries — `mpblas_<backend>` (BLAS) and `mplapack_<backend>` (LAPACK) — with
namespaced targets and an installed CMake package, so downstream CMake projects
can consume MPLAPACK via `find_package(mplapack)` or `FetchContent`.

Unlike the autotools build, the CMake build **does not download or build the
precision dependencies**. GMP, MPFR, MPC and QD are located with
`find_package`; install their development packages first (e.g. on Debian/Ubuntu:
`libgmp-dev libmpfr-dev libmpc-dev libqd-dev`).

## Quick start

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
cmake --install build --prefix /path/to/prefix
```

## Backend options

| Option | Default | Notes |
|---|---|---|
| `MPLAPACK_ENABLE_DOUBLE`    | ON  | binary64 reference |
| `MPLAPACK_ENABLE_MPFR`      | ON  | needs MPFR, MPC, GMP (C++) |
| `MPLAPACK_ENABLE_GMP`       | ON  | needs GMP (C++) |
| `MPLAPACK_ENABLE_QD`        | ON  | needs QD |
| `MPLAPACK_ENABLE_DD`        | ON  | needs QD; built with `-ffp-contract=off` |
| `MPLAPACK_ENABLE_BINARY128` | ON  | GCC/Intel; not supported by Clang |
| `MPLAPACK_ENABLE_BINARY80`  | OFF | x86/x86_64 only |

The binary128/binary80 type, I/O and math modes are auto-detected at configure
time (mirroring `configure.ac`) and written into the generated
`mplapack_config.h`.

Other options: `MPLAPACK_BUILD_EXAMPLES`, `MPLAPACK_BUILD_TESTS`,
`MPLAPACK_BUILD_BENCHMARKS` (all OFF by default), and the standard
`BUILD_SHARED_LIBS`.

## Consuming via find_package

```cmake
find_package(mplapack REQUIRED)
target_link_libraries(myapp PRIVATE mplapack::mplapack_mpfr)
```

Linking `mplapack::mplapack_<backend>` transitively provides the MPBLAS layer,
the include directories, the backend-selecting `-D___MPLAPACK_BUILD_WITH_<...>___`
definition, and the precision dependencies — so your sources can simply:

```cpp
#include <mpblas_mpfr.h>
#include <mplapack_mpfr.h>
```

## Consuming via FetchContent

```cmake
include(FetchContent)
set(MPLAPACK_ENABLE_MPFR ON CACHE BOOL "" FORCE)   # pick the backends you need
FetchContent_Declare(mplapack
  GIT_REPOSITORY https://github.com/nakatamaho/mplapack
  GIT_TAG        master)
FetchContent_MakeAvailable(mplapack)

target_link_libraries(myapp PRIVATE mplapack::mplapack_mpfr)
```

Disable the backends you do not need (each backend compiles ~1000 sources).

## Tests

`-DMPLAPACK_BUILD_TESTS=ON` builds two suites and registers them with CTest:

- **`mplapack/test`** — the LAPACK-derived `lin`/`eig` drivers per backend, each
  fed its `.in` dataset on stdin.
- **`mpblas/test`** — one program per BLAS routine per backend, comparing against
  a reference Fortran BLAS. This suite needs the **MPFR backend** (its random/
  comparison support always uses MPFR) and a reference BLAS located via
  `find_package(BLAS)` (e.g. reference BLAS, OpenBLAS or FlexiBLAS); it is
  skipped if either is missing.
- **`mplapack/test/compare`** — one program per MPLAPACK routine per backend,
  comparing against the reference Fortran LAPACK (`find_package(LAPACK)`); same
  MPFR/BLAS requirements as above. A few cases are numerically borderline against
  the double-precision reference and may report a mismatch (a property of the
  tests themselves, not the build).

```sh
cmake -S . -B build -DMPLAPACK_BUILD_TESTS=ON
cmake --build build -j
ctest --test-dir build
```

## Scope vs. the autotools build

Not ported: the bundled third-party builds (GMP/MPFR/MPC/QD/OpenBLAS), the
CUDA-optimized kernels (`*/optimized/*/cuda`), the `_opt` optimized-library
variants, the Fable Fortran-to-C++ regeneration pipeline, and the
reference-BLAS/OpenBLAS comparison benchmarks.
