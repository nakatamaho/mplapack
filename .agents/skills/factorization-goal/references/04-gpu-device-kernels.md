# GOAL 04: GPU building blocks for device-resident factorizations (dd, td, qd)

Prerequisites: goal 01. Read 00-overview.md ("GPU design") first.

## Scope

Internal device code only (no new public entry point yet), in
mpblas/optimized/common (headers shared by dd/td/qd) and the dd/td/qd
cuda/ directories.

## Tasks

1. Device matrix residency helper (namespace mplapack_lapack_gpu): allocate,
   upload, download a column-major matrix in the layout the existing
   kernels use; device-pointer variants of the tiled Rgemm, Rtrsm, Rsyrk /
   Rherk launches (refactor gemm_cuda / blas_gpu so that the host-copy
   wrapper and the device-pointer path share the kernels).
2. Device kernels with reference semantics: iRamax / iCamax over a column
   (block reduction, first index on ties), row swap and Rlaswp on a device
   matrix, column scaling with the sfmin test (reciprocal multiply or
   accurate_div), sqrt for the Cholesky diagonal (bit-identical to the
   C++ type's sqrt — add a device version next to dd_div if missing and
   test it against libQD3), the Rgetrf2 and Rpotrf2 recursions driven from
   the host on device data.
3. Host emulation build of all of it (MPLAPACK_GEMM_CUDA_HOST_EMULATION).

## Acceptance criteria

- A device test (CMake, SKIP_RETURN_CODE 77) and its host emulation variant
  (autotools make check and CMake) compare each kernel with the CPU
  reference operation by operation: iRamax indices, Rlaswp results, the
  panel factorization of an m x nb block (Rgetrf2 / Rpotrf2) bit for bit.
- Existing Rgemm / Rblas / Cblas device and host tests still pass (dd, td,
  qd).
- check_unique_symbols.sh passes; no new public symbol (`nm -D` diff).
- One commit: "CUDA: device-resident building blocks for getrf/potrf (factorization goal 04)".
