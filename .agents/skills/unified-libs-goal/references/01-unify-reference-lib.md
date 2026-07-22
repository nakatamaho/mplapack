# GOAL 01: Merge libmpblas_<b> + libmplapack_<b> into a single libmplapack_<b>

Prerequisite: read .agents/skills/unified-libs-goal/references/00-overview.md. Work on branch topic/unified-libs.

## Scope

Reference layer only. After this goal, for every enabled backend <b> exactly
one library `libmplapack_<b>` is built and installed, containing all sources
from `mpblas/reference/` and `mplapack/reference/`. `libmpblas_<b>` ceases to
exist in both build systems.

## Tasks — autotools

1. In `mplapack/reference/Makefile.am`, extend `MPLAPACK_SOURCES` (or the
   per-backend `_SOURCES`) to also include every `../../mpblas/reference/*.cpp`.
   Prefer generating/maintaining a single shared source list; if the file uses
   an explicit list, add the mpblas files explicitly (78 files) rather than
   wildcards — automake does not expand wildcards portably.
2. Transplant everything `mpblas/reference/Makefile.am` attaches per backend
   (CPPFLAGS backend macro, DD_CXXFLAGS / -ffp-contract handling, GMP/MPFR/QD
   include dirs, LDFLAGS incl. macOS/MinGW conditionals, -version-info) onto
   the merged libmplapack_<b> in `mplapack/reference/Makefile.am`, keeping the
   stricter variant when both layers differ.
3. Turn `mpblas/reference/Makefile.am` into a dist-only Makefile.am:
   `lib_LTLIBRARIES` removed, sources kept via EXTRA_DIST (they are compiled
   from mplapack/reference now). Keep the directory in SUBDIRS so `make dist`
   still ships the sources.
4. Grep the whole tree for `-lmpblas` and `libmpblas` and update every link
   line that refers to reference libraries: `mpblas/test/`, `mplapack/test/`,
   `examples/`, `benchmark/` where applicable. mpblas tests now link
   `-lmplapack_<b>`.
5. configure.ac: remove/adjust any AC_CONFIG_FILES or conditionals that only
   served libmpblas_<b> creation; keep the Makefile outputs for the dist-only
   directories.

## Tasks — CMake

1. In `cmake/MplapackBackends.cmake`, collapse `mplapack_add_backend()` to
   create ONE target `mplapack_<backend>` from
   `${MPBLAS_SOURCES};${MPLAPACK_SOURCES}` with the union of the current
   properties of the two targets. Delete the mpblas_<backend> target and its
   alias. Keep alias `mplapack::mplapack_<backend>`.
2. Top-level CMakeLists.txt: drop mpblas_* from MPLAPACK_INSTALL_TARGETS,
   per-backend target_link_libraries lines now apply once (e.g. GMP::GMPXX on
   mplapack_gmp only), dd's -ffp-contract=off applies to the merged target.
3. Update `cmake/MplapackInstall.cmake` / `MplapackPackage.cmake` /
   `mplapackConfig.cmake.in` / `cmake/mplapack.pc.cmake.in` so exported
   targets, pc Libs lines, and COMPONENT lists no longer mention mpblas_*.
4. Update CMake test/example targets that link mpblas_* to link mplapack_*.

## New verification asset

Create `misc/check_unique_symbols.sh` exactly as specified in 00-overview.md,
chmod +x, and wire it:
- autotools: a `check-local` hook (or a script target) in the toplevel that
  runs it over `$(libdir)`-staged or .libs/ archives of every built backend.
- CMake: `add_test(NAME unique_symbols_<b> COMMAND ...)` per backend when
  MPLAPACK_BUILD_TESTS=ON.

## Acceptance criteria

- `./gen_configure.sh && ./configure <minimal flags for dd,double> && make -j`
  produces only libmplapack_dd.*, libmplapack_double.* (no libmpblas_*).
- `cmake -B build -DMPLAPACK_BUILD_TESTS=ON` (dd+double on, others off)
  configures, builds, and `ctest` passes including unique_symbols tests.
- At least one existing mpblas test binary and one mplapack test binary link
  and run against the merged library in BOTH build systems.
- `grep -R "mpblas_dd\b" --include=Makefile.am --include=*.cmake .` returns
  no library-producing hits (dist/test-source references are fine).
- One commit: "Merge mpblas into libmplapack (goal 01)".
