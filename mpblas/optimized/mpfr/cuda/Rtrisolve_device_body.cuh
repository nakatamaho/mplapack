#ifndef MPLAPACK_MPFR_TRISOLVE_DEVICE_BODY_CUH
#define MPLAPACK_MPFR_TRISOLVE_DEVICE_BODY_CUH
#include "Rtrisolve_kernel_cuda.h"
#include "Rinverse_device_body.cuh"
namespace mplapack_mpfr_inverse {
template<class Tag,bool Complex> int run_tri(const mplapack_trisolve::shape &s,long prec,value<Tag,Complex> *a,value<Tag,Complex> *x) {
    using V=value<Tag,Complex>;
    using A=typename std::conditional<Complex,mplapack_trisolve::Carithmetic<V>,mplapack_trisolve::Rarithmetic<V>>::type;
#ifndef MPLAPACK_MPFR_CUDA_HOST_EMULATION
    int devices=0;if(cudaGetDeviceCount(&devices) || !devices) return 1;
    constexpr bool fixed=std::is_same<Tag,factor_fixed_tag<256>>::value || std::is_same<Tag,factor_fixed_tag<512>>::value ||
        std::is_same<Tag,factor_fixed_tag<768>>::value || std::is_same<Tag,factor_fixed_tag<1024>>::value || std::is_same<Tag,factor_fixed_tag<2048>>::value;
    mpfr_stack_guard stack(inverse_stack_bytes(prec,fixed));if(stack.error) return stack.error;
    factor_buffer<char> arena;factor_buffer<size_t> tops;
    size_t slab=16384*((prec+1023)/1024);
    if(arena.allocate(factor_pool*slab) || tops.allocate(factor_pool) || cudaMemset(tops.p,0,factor_pool*sizeof(size_t))) return 1;
    complex_arena_guard guard(arena.p,slab,tops.p,prec);if(guard.error) return guard.error;
#else
    long previous=mplapack_mpfr_rt_prec_h;mplapack_mpfr_rt_prec_h=prec;
    struct precision_guard { long p;~precision_guard() { mplapack_mpfr_rt_prec_h=p; } } guard{previous};
#endif
    factor_buffer<V> da,dx,p;
    if(da.allocate(s.count()) || dx.allocate(s.n) || p.allocate(s.n+1) || da.upload(a,s.count()) || dx.upload(x,s.n)) return 1;
    launcher ctx;if(ctx.error) return ctx.error;
#ifndef MPLAPACK_MPFR_CUDA_HOST_EMULATION
    // The factor launcher uses nonblocking streams. Finish the pageable H2D
    // uploads and its constructor's status initialization before consumption.
    if(cudaStreamSynchronize(nullptr)) return 1;
#endif
    mplapack_trisolve::schedule<V,A>(ctx,s,da.p,dx.p,p.p);
    ctx.finish();if(ctx.error) return ctx.error;
    if(dx.download(x,s.n)) return 1;
#ifndef MPLAPACK_MPFR_CUDA_HOST_EMULATION
    if(guard.restore() || stack.restore()) return 1;
#endif
    return 0;
}
}
#endif
