# GOAL 02: Single libmplapack_<b>_opt per backend + mplapack/optimized skeleton

Prerequisite: goal 01 merged on topic/unified-libs. Read 00-overview.md.

## Scope

Replace every `libmpblas_<b>_opt` with a self-contained `libmplapack_<b>_opt`
containing: optimized BLAS dispatchers + their helpers, the non-shadowed
reference BLAS files, ALL reference LAPACK files, and (empty for now) the new
`mplapack/optimized/<b>/` layer. CPU only — CUDA/OpenCL are goal 03 and must
be EXCLUDED from these libraries.

Backends with an optimized layer today: dd, qd, mpfr, gmp, double, binary80,
binary128 (see mpblas/optimized/*/). Create libmplapack_<b>_opt for each of
those; backends without an optimized dir get none.

## Shadowing (normative)

Per 00-overview.md: for each dispatcher `mpblas/optimized/<b>/X.cpp`
(currently X ∈ {Raxpy, Rcopy, Rdot, Rgemm}; derive the set by listing the
directory, do NOT hardcode), drop `mpblas/reference/X.cpp` from the _opt
source list. `*_ref.cpp` and `openmp/*_omp.cpp` are helpers: always included,
never shadowing. Everything else comes from the reference lists unchanged.

## Tasks — skeleton for optimized LAPACK

1. Create `mplapack/optimized/README.md` documenting: purpose (backend-specific
   optimized LAPACK drivers), the shadowing rule (a file `X.cpp` here replaces
   `mplapack/reference/X.cpp` in the _opt library of that backend), the
   `X_ref` naming convention if a fallback is kept, and the requirement that
   every added file keep the public signature identical to the reference one.
2. Create `mplapack/optimized/<b>/` for each optimized backend containing a
   `.gitkeep` (or the README only at top level plus per-backend .gitkeep) and
   a `Makefile.am` that today builds nothing but is wired into configure.ac
   AC_CONFIG_FILES, so adding a source file later is a one-line change.
3. In both build systems, the _opt source computation MUST already consult
   `mplapack/optimized/<b>/*.cpp` and apply shadowing against
   `mplapack/reference/`, so the skeleton is live, not decorative.
   - CMake: `file(GLOB ... CONFIGURE_DEPENDS mplapack/optimized/<b>/*.cpp)`
     then a helper `mplapack_shadow_sources(<out> <base_list> <override_list>)`
     in MplapackBackends.cmake that removes stem-matching entries. Reuse the
     same helper for the mpblas layer.
   - autotools: per-backend explicit lists in
     `mpblas/optimized/<b>/Makefile.am` (renamed responsibility: it now builds
     libmplapack_<b>_opt). Add a comment block naming the shadowing rule and
     listing which reference files are intentionally omitted.

## Tasks — build systems

1. autotools: in each `mpblas/optimized/<b>/Makefile.am`, rename the product
   to `libmplapack_<b>_opt.la`, append all `../../../mplapack/reference/*.cpp`
   sources (same list goal 01 consolidated — factor it into a common included
   fragment, e.g. `mplapack/reference/sources.am` created in goal 01 or here,
   `include $(top_srcdir)/...` style, to avoid a ~970-entry copy per backend),
   plus `$(top_srcdir)/mplapack/optimized/<b>/*.cpp` entries as they appear
   (initially none). Remove the `SUBDIRS += cuda` hook here (moved in goal 03).
2. CMake: add `mplapack_add_opt_backend(<b>)` in MplapackBackends.cmake:
   target `mplapack_<b>_opt`, alias `mplapack::mplapack_<b>_opt`, OpenMP CXX
   linked PRIVATE (the dispatchers use OpenMP), same PUBLIC macro/include
   properties as the base backend, gated by a new
   `option(MPLAPACK_ENABLE_OPT "Build optimized libraries" ON)` combined with
   the backend option. Install/export it like the base targets.
3. Tests: mpblas/test contains _opt test variants (grep for `_opt`); repoint
   their link lines to `-lmplapack_<b>_opt` / `mplapack::mplapack_<b>_opt`.
   An _opt library must pass the SAME reference test suite binaries — wire at
   least the dd and double _opt test executables into `make check` and ctest.

## Acceptance criteria

- Both build systems produce libmplapack_dd_opt and libmplapack_double_opt
  (others per enabled options); no libmpblas_*_opt anywhere.
- `misc/check_unique_symbols.sh` passes on every _opt archive (this is the
  regression the old layout could not guarantee).
- Linking a test program against ONLY `-lmplapack_dd_opt` (no other mplapack
  libs) resolves: proves self-containment.
- Dropping a file `mplapack/optimized/dd/Rpotrf.cpp` (create a temporary copy
  of the reference one, build, confirm it shadows — check the object list —
  then delete it before committing) demonstrates the skeleton is live.
- One commit: "Unify optimized layer into libmplapack_<b>_opt (goal 02)".
