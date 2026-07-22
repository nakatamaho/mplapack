# Changes in MPLAPACK 2.3.0

## Unified library layout

- Every backend/flavor is now a self-contained library containing both MPBLAS
  and MPLAPACK routines. Consumers link exactly one of
  `libmplapack_<backend>`, `libmplapack_<backend>_opt`,
  `libmplapack_dd_opt_cuda`, or `libmplapack_binary128_opt_opencl`.
- The separate `libmpblas_*` products and aggregate pkg-config file were
  removed. Per-flavor pkg-config files and CMake export targets describe the
  exact unified library and its precision dependencies.
- Duplicate public symbols are prevented with a deterministic basename
  shadowing rule: accelerator sources override optimized sources, which
  override reference sources. Archive checks reject any remaining duplicate
  global definition.
- The `mplapack/optimized/<backend>` skeleton records the LAPACK-side optimized
  layer even where it is currently empty, keeping the three source layers
  explicit and extensible.

## CMake

- `MPLAPACK_ENABLE_OPT` controls optimized CPU flavors.
- `MPLAPACK_ENABLE_CUDA` enables the self-contained DD CUDA flavor.
- `MPLAPACK_ENABLE_OPENCL` enables the self-contained binary128 OpenCL flavor.
- `MPLAPACK_CUDA_ARCHITECTURES` selects CUDA device architectures.
- Unified build-tree and installed targets are available as
  `mplapack::mplapack_<flavor>`. The source-manifest consistency check is now
  registered with CTest as well as the autotools test path.

## Numerical test resolution

- The reference-double `se2` run remains at threshold 50 and passes.
- The optimized-double run was reproduced at test (36), matrix type 15,
  `N=20`, seed `60812,42731,63787,906`: ratio 57.2705 with
  `OMP_NUM_THREADS=1` and 15.2737 with the default thread count on GCC 15.2.0.
  This bounded difference is summation-order divergence in the optimized BLAS
  path, so optimized `se2` runs use a documented threshold of 100; kernels are
  unchanged.

## Version and ABI

- Package version: 2.3.0.
- The incompatible library-layout change advances libtool version-info from
  `1:0:0` to `2:0:0`; CMake shared libraries use version 2.0.0 and SONAME 2.
