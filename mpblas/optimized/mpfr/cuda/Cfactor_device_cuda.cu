// Fixed complex policy; scheduler, panels, swaps and streams are shared.
#include "Cfactor_device_body.cuh"
namespace mplapack_mpfr_cuda {
using namespace mplapack_blas_dev;
bool complex_factor_memory_fits(long m,long n,long prec,bool fixed) {
#ifndef MPLAPACK_MPFR_CUDA_HOST_EMULATION
    size_t free=0,total=0; int device=0; cudaDeviceProp prop{};
    if(cudaGetDevice(&device) || cudaGetDeviceProperties(&prop,device) || cudaMemGetInfo(&free,&total)) return false;
    // CUDA reserves the per-thread stack for every potentially resident thread,
    // not merely our bounded launch pool. Include the full requested reserve.
    size_t overhead=size_t(prop.multiProcessorCount)*prop.maxThreadsPerMultiProcessor*(fixed?complex_fixed_stack_bytes(prec):complex_runtime_stack_bytes(prec));
    overhead+=size_t(2*factor_pool)*(16384*((prec+1023)/1024)+sizeof(size_t));
    overhead+=1048576+(size_t(std::min(m,n))+(size_t(m)+127)/128+2)*sizeof(long);
    size_t nl=(prec+63)/64;
    if(!fixed) { size_t tier=2; while(tier<nl) tier*=2; nl=tier; }
    size_t bytes=2*(nl*8+24);
    return free>=overhead && size_t(n)<=(free-overhead)/size_t(m)/bytes;
#else
    (void)m; (void)n; (void)prec; (void)fixed; return true;
#endif
}
#define INSTANTIATE(PB) template int factor_complex<factor_fixed_tag<PB>>(const factor_shape &,long,factor_complex_value<PB> *,const factor_complex_value<PB> &,long *,long &,factor_times &);
MPLAPACK_MPFR_CUDA_FIXED_PRECISIONS(INSTANTIATE)
#undef INSTANTIATE
}
