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
5. CMake support-tier statement: README.cmake.md must state the tier
   explicitly: CMake builds are fully exercised by the test suite as of goal
   04, but autotools remains the canonical system for release tarballs and
   packaging, and the CMake interface (option names, export targets) may
   still change without notice (experimental). Do not silently drop or keep
   the word "experimental" without this sentence.
6. CTest manifest gate: add an add_test invocation running
   misc/check_source_manifests.sh with the source dir, so manifest drift
   fails both build systems, not only autotools check-local.
7. Resolve the standing CTest failure xeigtstR_double_opt_se2
   (ratio 57.2705 > THRESH 50, recorded in goal 04):
   a. Re-run the driver with OMP_NUM_THREADS=1 and with the default thread
      count; record both ratios. Identify the failing test line (driver,
      matrix type, N) from the xeigtst output.
   b. Confirm the reference double run of the same input passes.
   c. If (b) passes and the exceedance is stable and O(100) at most, the
      cause is summation-order divergence in the optimized BLAS path
      (dispatchers + OpenMP kernels). Then raise THRESH to 100 ONLY in the
      input file consumed by the _opt test runs, with an in-file comment
      recording the observed ratio, cause, compiler, and date. Reference
      inputs keep THRESH 50. Do not modify kernels.
   d. If (b) fails, or the ratio is thread-count-unstable beyond O(100),
      STOP: leave the test red, record findings in the commit body, and do
      not change any threshold.
   Record the outcome in CHANGES.2.3.0.md either way.
8. Docker smoke: if the docker matrix is runnable in the environment, build
   ONE image (the default ubuntu one) to prove packaging file lists are
   consistent; otherwise mark UNVERIFIED in the commit body.

## Acceptance criteria

- Zero failing CTests after task 7 (or an explicit task-7d stop recorded).
- Fresh clone of the branch: autotools full build (default backends) and
  CMake build both succeed; `make distcheck` passes if it passed on master
  (check first; if it was already broken, do not fix here, note it).
- MIGRATION.md table complete and consistent with actual installed file names
  (verify against the install trees, not from memory).
- One commit: "Version bump, migration docs, periphery updates (goal 05)".
- Open items for a follow-up PR description: summarize all UNVERIFIED marks
  from goals 03/05 into .agents/skills/unified-libs-goal/references/REPORT.md.
