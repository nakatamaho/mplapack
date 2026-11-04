#include "Rinverse_device_body.cuh"
#include "Rblas1_device_body.cuh"
namespace mplapack_mpfr_inverse {
thread_local bool failed_host=false;
#if defined(__CUDACC__)
__device__ bool failed_device[8192];
#endif
bool memory_fits(size_t count,size_t bytes,long n,long prec,bool fixed) {
#ifndef MPLAPACK_MPFR_CUDA_HOST_EMULATION
    size_t free=0,total=0; int device=0; cudaDeviceProp prop{};
    if(cudaGetDevice(&device) || cudaGetDeviceProperties(&prop,device) || cudaMemGetInfo(&free,&total)) return false;
    size_t stack=inverse_stack_bytes(prec,fixed),inherited=0;
    if(cudaDeviceGetLimit(&inherited,cudaLimitStackSize)) return false;
    size_t overhead=size_t(prop.multiProcessorCount)*prop.maxThreadsPerMultiProcessor*std::max(stack,inherited);
    overhead+=size_t(factor_pool)*(16384*((prec+1023)/1024)+sizeof(size_t));
    overhead+=1048576+size_t(n)*sizeof(long);
    return free>=overhead && count<=(free-overhead)/bytes;
#else
    (void)count;(void)bytes;(void)n;(void)prec;(void)fixed; return true;
#endif
}
#define INSTANTIATE(PB) \
template int run<factor_fixed_tag<PB>,false>(const shape &,long,value<factor_fixed_tag<PB>,false> *,value<factor_fixed_tag<PB>,false> *,const long *,long &,factor_times &); \
template int run<factor_fixed_tag<PB>,true>(const shape &,long,value<factor_fixed_tag<PB>,true> *,value<factor_fixed_tag<PB>,true> *,const long *,long &,factor_times &);
MPLAPACK_MPFR_CUDA_FIXED_PRECISIONS(INSTANTIATE)
#undef INSTANTIATE
#define INSTANTIATE_BLAS1(PB) \
template int run_blas1<factor_fixed_tag<PB>,false>(mplapack_blas1::operation,long,long,value<factor_fixed_tag<PB>,false> *,value<factor_fixed_tag<PB>,false> *,value<factor_fixed_tag<PB>,false>,mplapack_blas1::state<value<factor_fixed_tag<PB>,false>> &); \
template int run_blas1<factor_fixed_tag<PB>,true>(mplapack_blas1::operation,long,long,value<factor_fixed_tag<PB>,true> *,value<factor_fixed_tag<PB>,true> *,value<factor_fixed_tag<PB>,false>,mplapack_blas1::state<value<factor_fixed_tag<PB>,false>> &);
MPLAPACK_MPFR_CUDA_FIXED_PRECISIONS(INSTANTIATE_BLAS1)
#undef INSTANTIATE_BLAS1
}
