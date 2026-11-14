// Minimal CUDA runtime for deterministic stack guard error-path unit tests.
#ifndef MPLAPACK_STACK_GUARD_CUDA_STUB_H
#define MPLAPACK_STACK_GUARD_CUDA_STUB_H
#include <cstddef>
enum cudaLimit { cudaLimitStackSize };
enum cudaError_t { cudaSuccess=0, cudaErrorInvalidValue=1 };
cudaError_t cudaDeviceGetLimit(size_t *,cudaLimit);
cudaError_t cudaDeviceSetLimit(cudaLimit,size_t);
cudaError_t cudaGetLastError();
#endif
