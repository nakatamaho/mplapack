# P0 CMake Test Level

CMake was used at the locked MPLAPACK SHA to configure, build, install, and
exercise downstream package consumers for all seven backends. It was not used
as the P0 numerical baseline.

The exact configure command was:

```sh
PKG_CONFIG_PATH=/tmp/mplapack-p0-final/lib/pkgconfig cmake \
  -S /tmp/mplapack-p0-base \
  -B /tmp/mplapack-p0-cmake-build \
  -DCMAKE_INSTALL_PREFIX=/tmp/mplapack-p0-cmake-final \
  -DCMAKE_PREFIX_PATH=/tmp/mplapack-p0-final \
  -DMPLAPACK_ENABLE_GMP=ON \
  -DMPLAPACK_ENABLE_MPFR=ON \
  -DMPLAPACK_ENABLE_QD=ON \
  -DMPLAPACK_ENABLE_DD=ON \
  -DMPLAPACK_ENABLE_DOUBLE=ON \
  -DMPLAPACK_ENABLE_BINARY80=ON \
  -DMPLAPACK_ENABLE_BINARY128=ON \
  -DMPLAPACK_ENABLE_OPT=ON \
  -DMPLAPACK_ENABLE_CUDA=OFF \
  -DMPLAPACK_ENABLE_OPENCL=ON \
  -DMPLAPACK_BUILD_EXAMPLES=OFF \
  -DMPLAPACK_BUILD_TESTS=OFF \
  -DMPLAPACK_BUILD_BENCHMARKS=OFF
```

The exact build and install commands were:

```sh
cmake --build /tmp/mplapack-p0-cmake-build -j32
cmake --install /tmp/mplapack-p0-cmake-build
```

All three commands passed. `MPLAPACK_BUILD_TESTS=OFF` is the actual CMake test
level. The complete P0 numerical baseline instead comes from the supported
Autotools MPBLAS, LIN, and EIG test sets for GMP and the ordinary 512-bit MPFR
profile. CMake package consumers are recorded independently in
`abi/manifest.json`.
