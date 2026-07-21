# GOAL 05: Versioning, migration docs, and periphery cleanup

Prerequisite: goals 01-04 merged on topic/unified-libs. Read 00-overview.md.

## Scope

Finish the break cleanly: version numbers, migration documentation, and every
peripheral consumer inside the repo.

## Tasks

1. Versioning:
   - configure.ac + CMakeLists.txt project version -> 2.3.0 (keep them in
     sync; note the "keep in sync" comment already present in CMakeLists.txt).
   - libtool -version-info: interfaces removed -> current++, age=0
     (1:0:0 -> 2:0:0); mirror in MPLAPACK_LIB_VERSION/SOVERSION in CMake.
2. MIGRATION.md: add a "2.3.0 unified libraries" section with a link-line
   translation table:
   | old                                   | new                              |
   |---------------------------------------|----------------------------------|
   | -lmplapack_dd -lmpblas_dd             | -lmplapack_dd                    |
   | -lmpblas_dd_opt                       | -lmplapack_dd_opt                |
   | -lmpblas_dd_opt -lmpblas_dd_cuda      | -lmplapack_dd_opt_cuda           |
   | -lmpblas_binary128_opt_opencl (+opt)  | -lmplapack_binary128_opt_opencl  |
   plus the rule "link exactly one mplapack library per backend/flavor" and
   the CMake target rename table (mplapack::mpblas_<b> -> gone).
3. CHANGES.2.3.0.md: new file, factual list of the layout change, the
   duplicate-symbol fix, the shadowing rule, the mplapack/optimized skeleton,
   and the new CMake options (MPLAPACK_ENABLE_OPT/CUDA/OPENCL,
   MPLAPACK_CUDA_ARCHITECTURES).
4. Periphery sweep — update link lines and docs in:
   - examples/ (all Makefile/scripts), benchmark/, docker/ (all images),
     packaging/ (spec/deb rules: file lists lose libmpblas*, gain new names),
     release/ scripts, README.md, README.cmake.md, doc/ (manual mentions of
     libmpblas), mplapack.pc.in.
   - `git grep -n "mpblas"` at the end: remaining hits must be source dirs,
     historical CHANGES, or MIGRATION — nothing operational.
5. Docker smoke: if the docker matrix is runnable in the environment, build
   ONE image (the default ubuntu one) to prove packaging file lists are
   consistent; otherwise mark UNVERIFIED in the commit body.

## Acceptance criteria

- Fresh clone of the branch: autotools full build (default backends) and
  CMake build both succeed; `make distcheck` passes if it passed on master
  (check first; if it was already broken, do not fix here, note it).
- MIGRATION.md table complete and consistent with actual installed file names
  (verify against the install trees, not from memory).
- One commit: "Version bump, migration docs, periphery updates (goal 05)".
- Open items for a follow-up PR description: summarize all UNVERIFIED marks
  from goals 03/05 into .agents/skills/unified-libs-goal/references/REPORT.md.
