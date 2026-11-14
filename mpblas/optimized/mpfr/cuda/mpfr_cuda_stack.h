// MPFR-only CUDA stack policy; no change to dd/td/qd or mpc_cuda.
// Default: raise only. MPLAPACK_MPFR_CUDA_STACK_RESTORE=1 restores the
// original limit per call (read on every guard construction; only exact "1").
#ifndef MPLAPACK_MPFR_CUDA_STACK_H
#define MPLAPACK_MPFR_CUDA_STACK_H
#include <cstddef>
#include <cstdlib>
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
// needed so far and left there by default. cudaDeviceSetLimit
// synchronizes the device and reconfigures the per-thread stack; under
// driver 580.178.04 on the A100 one call costs tens of milliseconds, so the
// previous per-call set/restore pair (~100-200 ms) dominated every
// checked-precision MPFR BLAS call: the small trisolve device test ran
// 402.7 s and mpfr_cuda_Cblas_device_fixed 4691.8 s; with the limit set
// once per process they run 42.6 s and ~50 s.  The stack allowance has no
// effect on the arithmetic (values, rounding, zero signs); the peak
// allowance is the same as with the per-call guard.
__attribute__((visibility("hidden"))) inline bool stack_restore_enabled() {
    const char *v=std::getenv("MPLAPACK_MPFR_CUDA_STACK_RESTORE");
    return v && v[0]=='1' && v[1]=='\0';
}
struct __attribute__((visibility("hidden"))) mpfr_stack_guard {
    size_t old=0; bool changed=false; int error=0;
    const bool restore_per_call=stack_restore_enabled();
    mpfr_stack_guard(const mpfr_stack_guard &)=delete;
    mpfr_stack_guard &operator=(const mpfr_stack_guard &)=delete;
    explicit mpfr_stack_guard(size_t want) {
        if(!want) return;
        error=cudaDeviceGetLimit(&old,cudaLimitStackSize);
        if(!error && old<want) {
            error=cudaDeviceSetLimit(cudaLimitStackSize,want); changed=!error;
        }
    }
    int restore() {
        if(restore_per_call && changed) {
            changed=false;
            // Restoration is best effort, never a computation failure. Clear
            // the runtime error so it cannot poison the next kernel launch.
            if(cudaDeviceSetLimit(cudaLimitStackSize,old)!=cudaSuccess)
                (void)cudaGetLastError();
        }
        return 0;
    }
    ~mpfr_stack_guard() { (void)restore(); }
};
#endif
}
#endif
