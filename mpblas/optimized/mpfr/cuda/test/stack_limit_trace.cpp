// Interpose only in the Cblas test executable, including calls originating
// in the autotools CUDA shared library. Production libraries are unchanged.
#include "stack_limit_trace.h"
#include <cuda_runtime_api.h>
#include <dlfcn.h>
namespace {
thread_local bool active=false;
thread_local size_t previous=0;
thread_local mpfr_stack_trace::counts count{};
}
namespace mpfr_stack_trace {
void begin(size_t initial) { count={initial,0,0};previous=initial;active=true; }
counts end() { active=false;return count; }
}
extern "C" cudaError_t cudaDeviceSetLimit(cudaLimit which,size_t value) {
    using setter=cudaError_t (*)(cudaLimit,size_t);
    static setter real=reinterpret_cast<setter>(dlsym(RTLD_NEXT,"cudaDeviceSetLimit"));
    if(!real) return cudaErrorUnknown;
    cudaError_t error=real(which,value);
    if(active && which==cudaLimitStackSize && error==cudaSuccess) {
        size_t actual=0;
        if(cudaDeviceGetLimit(&actual,cudaLimitStackSize)==cudaSuccess) {
            if(actual>previous) ++count.raises;
            if(actual<previous) ++count.lowers;
            if(actual>count.peak) count.peak=actual;
            previous=actual;
        }
    }
    return error;
}
