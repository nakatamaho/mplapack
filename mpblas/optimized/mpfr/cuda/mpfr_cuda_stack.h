// MPFR-only CUDA stack policy; no change to dd/td/qd or mpc_cuda.
#ifndef MPLAPACK_MPFR_CUDA_STACK_H
#define MPLAPACK_MPFR_CUDA_STACK_H
#include <cstddef>
#ifndef MPLAPACK_MPFR_CUDA_HOST_EMULATION
#include <cuda_runtime.h>
#endif
namespace mplapack_mpfr_cuda {
// Cumulative factor-kernel + RN-call frames (sm_80): below 4 KiB at
// 256/512, below 8 KiB at 768/1024, below 16 KiB at 2048.
__attribute__((visibility("hidden"))) inline size_t fixed_stack_bytes(long p) {
    return p<=512?4096:p<=1024?8192:16384;
}
__attribute__((visibility("hidden"))) inline size_t runtime_stack_bytes(long p) {
    return p<=1024?8192:p<=4096?16384:65536;
}
// Complex panels also call MPC division/hypot. Keep their policy separate
// from real panels and BLAS; stack reservations are scoped by the same guard.
__attribute__((visibility("hidden"))) inline size_t complex_fixed_stack_bytes(long p) {
    return p<=1024?8192:16384;
}
__attribute__((visibility("hidden"))) inline size_t complex_runtime_stack_bytes(long p) {
    return p<=1024?8192:p<=2048?16384:65536;
}
#ifndef MPLAPACK_MPFR_CUDA_HOST_EMULATION
struct __attribute__((visibility("hidden"))) mpfr_stack_guard {
    size_t old=0; bool changed=false; int error=0;
    explicit mpfr_stack_guard(size_t want) {
        if(!want) return;
        error=cudaDeviceGetLimit(&old,cudaLimitStackSize);
        if(!error && old<want) {
            error=cudaDeviceSetLimit(cudaLimitStackSize,want); changed=!error;
        }
    }
    int restore() {
        if(!changed) return 0;
        int result=cudaDeviceSetLimit(cudaLimitStackSize,old);
        if(!result) changed=false;
        else (void)cudaGetLastError();
        return result;
    }
    ~mpfr_stack_guard() { (void)restore(); }
};
#endif
}
#endif
