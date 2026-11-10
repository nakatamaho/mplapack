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
// Goal 13 one-pass iamax for MPFR. Each lane scans one contiguous chunk of x
// (coalesced reads), then the partials are combined hierarchically with the
// min-index tie rule (order independent, so the reference index results).
template<class Tag,bool Complex,class V,class R> struct blas1_iamax_scan {
    const V *x;long n;long chunk;mplapack_blas1::state<R> *out;
    MPLAPACK_DEV_HD void operator()(long lane) const {
        using namespace mplapack_blas1;
        state<R> best;bool visited=false;
        const long i0=lane*chunk;
        if(i0>=n) { out[lane]=best;return; }
        const long i1=i0+chunk<n?i0+chunk:n;
        for(long i=i0;i<i1;++i) {
            state<R> s;
            if constexpr(Complex) s=complex_leaf<R,V>(iamax,x[i],V(0.0),i);
            else s=leaf<R>(iamax,x[i],R(0.0),i);
            if(!visited || s.re>best.re) { best=s;visited=true; }
        }
        out[lane]=best;
    }
};
template<class R> struct blas1_iamax_combine {
    const mplapack_blas1::state<R> *in;long m;long group;mplapack_blas1::state<R> *out;
    MPLAPACK_DEV_HD void operator()(long lane) const {
        using namespace mplapack_blas1;
        state<R> best;bool visited=false;
        const long i0=lane*group;
        if(i0>=m) { out[lane]=best;return; }
        const long i1=i0+group<m?i0+group:m;
        for(long i=i0;i<i1;++i)
            if(!visited || in[i].re>best.re || (in[i].re==best.re && in[i].index<best.index)) { best=in[i];visited=true; }
        out[lane]=best;
    }
};
template<class R> struct blas1_iamax_finish {
    long m;const mplapack_blas1::state<R> *in;mplapack_blas1::state<R> *out;
    MPLAPACK_DEV_HD void operator()(long) const {
        using namespace mplapack_blas1;
        state<R> best;bool visited=false;
        for(long i=0;i<m;++i)
            if(!visited || in[i].re>best.re || (in[i].re==best.re && in[i].index<best.index)) { best=in[i];visited=true; }
        *out=best;
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
    if(element) {
        if(dx.allocate(n) || dx.upload(x,n) || dy.allocate(n) || dy.upload(y,n)) return 1;
    } else if(op==mplapack_blas1::iamax) {
        // One pass (goal 13): each pool lane scans one contiguous chunk, then
        // the pool partials are combined hierarchically with the min-index tie
        // rule.  Same exact scores as the reference scan, so the index is the
        // reference's.
        if(dx.allocate(n) || dx.upload(x,n) || a.allocate(2*factor_pool)) return 1;
    } else {
        if(dx.allocate(n) || dx.upload(x,n) || ((element || product) && (dy.allocate(n) || dy.upload(y,n))) || a.allocate(count)) return 1;
    }
    launcher ctx;if(ctx.error) return ctx.error;
    if(element) ctx.launch(blas1_element_job<V,R>{op==axpy,dx.p,dy.p,alpha},n);
    else if(op==mplapack_blas1::iamax) {
        const long chunk=(n+factor_pool-1)/factor_pool;
        ctx.launch(blas1_iamax_scan<Tag,Complex,V,R>{dx.p,n,chunk,a.p},factor_pool);
        const long g1=64;
        const long m1=(factor_pool+g1-1)/g1;
        ctx.launch(blas1_iamax_combine<R>{a.p,factor_pool,g1,a.p+factor_pool},m1);
        ctx.launch(blas1_iamax_finish<R>{m1,a.p+factor_pool,a.p},1);
    }
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
