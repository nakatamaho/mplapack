# Unified-library migration report

This report collects the open verification items requested by goals 03 and 05.

## UNVERIFIED

- Goal 03: none. Its commit records PASS results for the CUDA and OpenCL
  build, single-library consumer, duplicate-symbol, and runtime checks.
- Goal 05 Docker smoke: UNVERIFIED. The execution environment has no `docker`
  command, so the default Ubuntu image could not be built.

## Baseline observations

- A fresh clone of master at `38797b77` lacks the ignored `ltmain.sh`;
  `./gen_configure.sh` reports an automake error but returns status 0.
  Running `libtoolize --copy --force` first provides the normal bootstrap.
- After that bootstrap, master `make distcheck` fails during `make dist`:
  `benchmark/Makefile` has no `distdir` rule. Goal 05 therefore does not
  attempt an unrelated distcheck fix.

## Goal 05 verification highlights

- A post-commit fresh clone builds all default CMake backends with `-j32`;
  all 12 resulting archives pass the duplicate-symbol gate.
- The same fresh clone, bootstrapped with `libtoolize --copy --force`, builds
  all seven CPU reference and optimized autotools backends with `make -j32`;
  all 14 resulting archives pass the duplicate-symbol gate.
- CMake DD/double reference and optimized archives pass the duplicate-symbol
  gate; all 198 registered CTests pass with `OMP_NUM_THREADS=1`.
- Installed CPU-flavor libraries and pkg-config files use only unified
  `mplapack_<flavor>` names and package version 2.3.0.
- Existing goal 03 install trees contain
  `libmplapack_dd_opt_cuda.a` and
  `libmplapack_binary128_opt_opencl.a`, matching the migration table.
- A CMake shared double build produces `libmplapack_double.so.2.0.0` with
  SONAME `libmplapack_double.so.2`.
