# GOAL 03: Accelerator libraries: libmplapack_dd_opt_cuda, libmplapack_binary128_opt_opencl

Prerequisite: goals 01-02 merged on topic/unified-libs. Read 00-overview.md.

## Scope

Two self-contained accelerator flavors, each a superset-replacement of the
corresponding _opt library. Old products `libmpblas_dd_cuda` and
`libmpblas_binary128_opt_opencl` are removed. A consumer links exactly one
library; never an accelerator library together with _opt or reference.

## Symbol resolution (this goal owns the duplicate-symbol fix)

Current hazards to eliminate, by construction:

1. `mpblas/optimized/dd/cuda/Rgemm.cu` defines public `Rgemm`, and the _opt
   dispatcher `mpblas/optimized/dd/Rgemm.cpp` also defines `Rgemm`.
   Resolution: in libmplapack_dd_opt_cuda, cuda/Rgemm.cu SHADOWS
   optimized/dd/Rgemm.cpp (accelerator > optimized). Same for Rsyrk if the
   cuda dir provides a public Rsyrk entry — inspect `Rsyrk.cu` and decide by
   the actual defined symbols (`nm` on the object), not by filename guesses.
   Kernel split files (Rgemm_fermi_*, Rgemm_tesla_*, Rsyrk_[NT][LU]_*,
   *_cuda.cu) are helpers: include, never shadow.
2. OpenCL: inspect `mpblas/optimized/binary128/opencl/Rgemm_opencl.cpp`. If it
   defines only `Rgemm_opencl` (helper naming), write a NEW dispatcher
   `mpblas/optimized/binary128/opencl/Rgemm.cpp` that defines public `Rgemm`,
   choosing between `Rgemm_opencl` and `Rgemm_ref`/CPU path with the same
   size-threshold pattern the dd OpenMP dispatcher uses (SINGLEOROMP-style
   constant; pick a documented threshold, make it overridable via an
   environment variable MPLAPACK_OPENCL_MIN_N with a sane default). If it
   already defines `Rgemm`, treat it as the top-layer entry directly.
3. After building each accelerator library, run
   `misc/check_unique_symbols.sh` on it — this is the acceptance gate.

## Tasks — autotools

1. Move CUDA build logic from `mpblas/optimized/dd/cuda/Makefile.am` into a
   Makefile.am producing `libmplapack_dd_opt_cuda.la`: sources = (goal-02 dd
   _opt list with cuda-shadowed stems removed) + cuda/*.cu compiled via the
   existing cudalt.py rule. Keep NVCC flag plumbing (AMPERE_FLAGS,
   HOPPER_FLAGS, NVCC_GLIBC_FEATURE_FLAGS) as-is; do not modernize nvcc
   handling in this goal.
2. Same for OpenCL: `mpblas/optimized/binary128/opencl/Makefile.am` produces
   `libmplapack_binary128_opt_opencl.la` = binary128 _opt list (opencl-shadow
   applied) + opencl sources; keep OPENCL_CFLAGS/LDFLAGS plumbing.
3. configure.ac: keep ENABLE_CUDA / OpenCL detection; ensure the accelerator
   dirs are in SUBDIRS only under their conditionals; remove the old library
   names everywhere (grep `dd_cuda`, `opt_opencl`).

## Tasks — CMake

1. Add `option(MPLAPACK_ENABLE_CUDA ...)` (default OFF) and
   `option(MPLAPACK_ENABLE_OPENCL ...)` (default OFF).
2. CUDA: guard `enable_language(CUDA)` behind the option; target
   `mplapack_dd_opt_cuda` built from the shadow-resolved source list plus the
   .cu files; set CMAKE_CUDA_ARCHITECTURES from a cache variable
   `MPLAPACK_CUDA_ARCHITECTURES` (default "70;80;90"); link CUDA::cudart via
   `find_package(CUDAToolkit REQUIRED)`. dd needs QD::QD and
   -ffp-contract=off on the C++ TUs as in the base dd target;
   `dd_real_cuda.h` stays a private header.
3. OpenCL: `find_package(OpenCL REQUIRED)` under the option; target
   `mplapack_binary128_opt_opencl` linking OpenCL::OpenCL.
4. Reuse `mplapack_shadow_sources()` from goal 02 for both; no second
   implementation of the shadowing logic may exist.
5. Install/export both targets with aliases `mplapack::mplapack_dd_opt_cuda`
   and `mplapack::mplapack_binary128_opt_opencl`; extend the package config
   so `find_package(mplapack COMPONENTS dd_opt_cuda)`-style consumption works
   (follow whatever component convention MplapackPackage.cmake established;
   if none, add a simple one and document it in README.cmake.md).

## Environment note

CUDA/OpenCL toolchains may be unavailable in the sandbox. In that case:
compile-verify everything that can be verified (CMake configure with the
options OFF must be bit-identical in behavior to before; with options ON,
configure must fail gracefully at find_package with a clear message), do a
dry parse of the autotools files (automake runs without nvcc), and mark the
on-GPU build as UNVERIFIED in the commit message body. Do not fake test
results.

## Acceptance criteria

- Old names gone: `git grep -n "mpblas_dd_cuda\|libmpblas_binary128_opt_opencl"`
  hits only CHANGES/MIGRATION docs.
- With both options OFF, build output of goals 01-02 is unchanged.
- Shadow lists are computed, not hardcoded, and the dd_opt_cuda source list
  visibly excludes optimized/dd/Rgemm.cpp (assert in a comment or a CMake
  message(VERBOSE) trace).
- unique-symbol check wired for the accelerator targets (runs when they are
  built).
- One commit: "Self-contained accelerator libraries + symbol dedup (goal 03)".
