// Reuse inverse-local checked arithmetic, scoped arena, stack and launch pool.
#ifndef MPLAPACK_MPFR_BLAS1_DEVICE_BODY_CUH
#define MPLAPACK_MPFR_BLAS1_DEVICE_BODY_CUH
#include "Rblas1_kernel_cuda.h"
#include "Rinverse_device_body.cuh"
#include "../../common/Rblas1_jobs_common.h"
namespace mplapack_mpfr_inverse {
template<int PB> MPLAPACK_MPFR_DEV_OP factor_real<PB> blas1_sqrt(const factor_real<PB> &a) {
    factor_real<PB> r;
    r.valid=a.valid && !cu_fp::cu_sqrt_rn<PB>(r.v,a.v);
    return r;
}
template<int NL> MPLAPACK_MPFR_DEV_OP rreal<NL> blas1_sqrt(const rreal<NL> &a) {
    rreal<NL> r;
#ifdef MPLAPACK_RT_ACTIVE
    MPLAPACK_RT_DECL(x);MPLAPACK_RT_DECL(z);MPLAPACK_RT_VIEW(x,a);
    MPLAPACK_RT(mpfr_custom_init_set)(z,MPLAPACK_RTU(MPFR_ZERO_KIND),0,rt_prec(),r.m);
    factor_arena_reset();MPLAPACK_RT(mpfr_sqrt)(z,x,MPLAPACK_RTU(MPFR_RNDN));MPLAPACK_RT_TAKE(r,z);
#endif
    return r;
}
template<class Tag> struct blas1_root {
    using R=value<Tag,false>;
    MPLAPACK_DEV_HD static R sqrt(const R &a) { return R(blas1_sqrt(a.v)); }
};
template<class V,class R> struct blas1_element_job {
    bool axpy;const V *x;V *y;R alpha;
    MPLAPACK_DEV_HD void operator()(long i) const {
        if(axpy) y[i]=y[i]+V(alpha)*x[i];else y[i]=x[i];
    }
};
template<class Tag,bool Complex> int run_blas1(mplapack_blas1::operation op,long n,long prec,
    value<Tag,Complex> *x,value<Tag,Complex> *y,value<Tag,false> alpha,mplapack_blas1::state<value<Tag,false>> &result) {
    using namespace mplapack_blas1;
    using R=value<Tag,false>;using V=value<Tag,Complex>;using S=state<R>;
    bool element=op==axpy || op==copy,product=op==dot || op==dotc || op==dotu;
    size_t count=Complex && op==mplapack_blas1::norm?2*size_t(n):size_t(n);
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
    factor_buffer<V> dx,dy;factor_buffer<S> a;
    if(dx.allocate(n) || dx.upload(x,n) || ((element || product) && (dy.allocate(n) || dy.upload(y,n))) || (!element && a.allocate(count))) return 1;
    launcher ctx;if(ctx.error) return ctx.error;
    if(element) ctx.launch(blas1_element_job<V,R>{op==axpy,dx.p,dy.p,alpha},n);
    else schedule<R,V,Complex,blas1_root<Tag>>(ctx,op,n,dx.p,dy.p,a.p);
    ctx.finish();if(ctx.error) return ctx.error;
    if(element?dy.download(y,n):a.download(&result,1)) return 1;
#ifndef MPLAPACK_MPFR_CUDA_HOST_EMULATION
    if(guard.restore() || stack.restore()) return 1;
#endif
    return 0;
}
}
#endif
