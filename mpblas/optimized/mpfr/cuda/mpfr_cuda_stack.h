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
// mpc_cuda v0.0.2-dev's cu_mpc_div needs more per-thread stack than
// abd9407's; the fixed-precision budget is raised accordingly.
__attribute__((visibility("hidden"))) inline size_t complex_fixed_stack_bytes(long p) {
    return p<=1024?16384:32768;
}
__attribute__((visibility("hidden"))) inline size_t complex_runtime_stack_bytes(long p) {
    return p<=1024?8192:p<=2048?16384:65536;
}
#ifndef MPLAPACK_MPFR_CUDA_HOST_EMULATION
// The device stack limit is raised to the largest size this process has
// needed so far and left there (restore() is a no-op).  cudaDeviceSetLimit
// synchronizes the device and reconfigures the per-thread stack; under
// driver 580.178.04 on the A100 one call costs tens of milliseconds, so the
// previous per-call set/restore pair (~100-200 ms) dominated every
// checked-precision MPFR BLAS call: the small trisolve device test ran
// 402.7 s and mpfr_cuda_Cblas_device_fixed 4691.8 s; with the limit set
// once per process they run 42.6 s and ~50 s.  The stack allowance has no
// effect on the arithmetic (values, rounding, zero signs); the peak
// allowance is the same as with the per-call guard.
struct __attribute__((visibility("hidden"))) mpfr_stack_guard {
    size_t old=0; bool changed=false; int error=0;
    explicit mpfr_stack_guard(size_t want) {
        if(!want) return;
        error=cudaDeviceGetLimit(&old,cudaLimitStackSize);
        if(!error && old<want) {
            error=cudaDeviceSetLimit(cudaLimitStackSize,want); changed=!error;
        }
    }
    int restore() { return 0; }
    ~mpfr_stack_guard() { }
};
#endif
}
#endif
