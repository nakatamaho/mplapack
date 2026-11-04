// BF CUDA reduction driver: staged publication; no host numerical work.
#ifndef MPLAPACK_BLAS1_CUDA_COMMON_CUH
#define MPLAPACK_BLAS1_CUDA_COMMON_CUH
#include "Rblas_dev_types.h"
#include "Rlapack_sqrt_common.h"
#include "Rblas1_jobs_common.h"
#include <vector>
#include <type_traits>
#if defined(__CUDACC__) && !defined(MPLAPACK_GEMM_CUDA_HOST_EMULATION)
#include <cuda_runtime.h>
#endif
namespace mplapack_blas1 __attribute__((visibility("hidden"))) {
#if defined(__CUDACC__) && !defined(MPLAPACK_GEMM_CUDA_HOST_EMULATION)
template<class J> __global__ void kernel(J job,long count) {
    for(long i=long(blockIdx.x)*blockDim.x+threadIdx.x;i<count;i+=long(blockDim.x)*gridDim.x) job(i);
}
struct launcher {
    int error=0;
    template<class J> void launch(J job,long n) {
        if(error) return;
        long blocks=(n+127)/128;if(blocks>65535) blocks=65535;
        kernel<<<unsigned(blocks),128>>>(job,n);
        error=cudaGetLastError();
    }
};
template<class V> struct buffer {
    V *p=nullptr;
    ~buffer() { cudaFree(p); }
    int allocate(size_t n) { return cudaMalloc(&p,n*sizeof(V)); }
    int up(const V *h,size_t n) { return cudaMemcpy(p,h,n*sizeof(V),cudaMemcpyHostToDevice); }
    int down(V *h,size_t n) { return cudaMemcpy(h,p,n*sizeof(V),cudaMemcpyDeviceToHost); }
};
#else
struct launcher {
    int error=0;
    template<class J> void launch(J job,long n) { for(long i=0;i<n;++i) job(i); }
};
template<class V> struct buffer {
    std::vector<V> v;V *p=nullptr;
    int allocate(size_t n) { v.resize(n);p=v.data();return 0; }
    int up(const V *h,size_t n) { std::copy(h,h+n,p);return 0; }
    int down(V *h,size_t n) { std::copy(p,p+n,h);return 0; }
};
#endif
template<class Ops> struct root {
    using R=mplapack_blas_dev::dreal<Ops>;
    MPLAPACK_DEV_HD static R sqrt(const R &x) { return R(mplapack_lapack_gpu::square_root<Ops>(x.v)); }
};
template<class Ops,bool Complex> int run_bf(operation op,long n,
    const typename std::conditional<Complex,typename mplapack_blas_dev::cgpu_types<Ops>::DC,typename mplapack_blas_dev::cgpu_types<Ops>::DR>::type *x,
    const typename std::conditional<Complex,typename mplapack_blas_dev::cgpu_types<Ops>::DC,typename mplapack_blas_dev::cgpu_types<Ops>::DR>::type *y,
    state<typename mplapack_blas_dev::cgpu_types<Ops>::DR> &result) {
    using R=typename mplapack_blas_dev::cgpu_types<Ops>::DR;
    using X=typename std::conditional<Complex,typename mplapack_blas_dev::cgpu_types<Ops>::DC,R>::type;
    long count=Complex && op==norm?2*n:n;
    const bool product=op==dot || op==dotc || op==dotu;
#if defined(__CUDACC__) && !defined(MPLAPACK_GEMM_CUDA_HOST_EMULATION)
    size_t free=0,total=0;
    if(cudaMemGetInfo(&free,&total) || free/sizeof(state<R>)<size_t(count) ||
       free-size_t(count)*sizeof(state<R>)<size_t(n)*sizeof(X)*(product?2:1)) return 1;
#endif
    buffer<X> dx,dy;buffer<state<R>> a;
    if(dx.allocate(n) || a.allocate(count) || dx.up(x,n) || (product && (dy.allocate(n) || dy.up(y,n)))) return 1;
    launcher ctx;
    schedule<R,X,Complex,root<Ops>>(ctx,op,n,dx.p,dy.p,a.p);
    if(ctx.error) return ctx.error;
    return a.down(&result,1);
}
}
#endif
