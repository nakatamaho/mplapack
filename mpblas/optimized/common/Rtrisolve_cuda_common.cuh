#ifndef MPLAPACK_TRISOLVE_CUDA_COMMON_CUH
#define MPLAPACK_TRISOLVE_CUDA_COMMON_CUH
#include "Rblas1_cuda_common.cuh"
#include "Rtrisolve_jobs.h"
namespace mplapack_trisolve __attribute__((visibility("hidden"))) {
template<class Ops,bool Complex> int run_bf(const shape &s,
    const typename std::conditional<Complex,typename mplapack_blas_dev::cgpu_types<Ops>::DC,typename mplapack_blas_dev::cgpu_types<Ops>::DR>::type *a,
    typename std::conditional<Complex,typename mplapack_blas_dev::cgpu_types<Ops>::DC,typename mplapack_blas_dev::cgpu_types<Ops>::DR>::type *x) {
    using R=typename mplapack_blas_dev::cgpu_types<Ops>::DR;
    using V=typename std::conditional<Complex,typename mplapack_blas_dev::cgpu_types<Ops>::DC,R>::type;
    using A=typename std::conditional<Complex,Carithmetic<V>,Rarithmetic<V>>::type;
    size_t count=s.count();
#if defined(__CUDACC__) && !defined(MPLAPACK_GEMM_CUDA_HOST_EMULATION)
    size_t free=0,total=0;
    if(cudaMemGetInfo(&free,&total) || (count+2*size_t(s.n)+1)>free/sizeof(V)) return 1;
#endif
    mplapack_blas1::buffer<V> da,dx,p;
    if(da.allocate(count) || dx.allocate(s.n) || p.allocate(s.n+1) || da.up(a,count) || dx.up(x,s.n)) return 1;
    mplapack_blas1::launcher ctx;
    schedule<V,A>(ctx,s,da.p,dx.p,p.p);
    if(ctx.error) return ctx.error;
    return dx.down(x,s.n);
}
}
#endif
