// Complex arithmetic policy for the real milestone's residency/scheduler.
// Kept MPFR-local: neither CPU nor dd/td/qd sees this specialization.
#ifndef MPLAPACK_MPFR_CFACTOR_DEVICE_BODY_CUH
#define MPLAPACK_MPFR_CFACTOR_DEVICE_BODY_CUH
#include "Cfactor_kernel_cuda.h"
#include <type_traits>
#define MPLAPACK_MPFR_FACTOR_COMPLEX
#include "Rfactor_device_body.cuh"
#undef MPLAPACK_MPFR_FACTOR_COMPLEX
#include "../../common/blas_dev/blas_dev.h"
namespace mplapack_mpfr_cuda {
namespace {
using namespace mplapack_blas_dev;
// MPFR's complex absolute value is a correctly rounded hypot, not a square
// root of two separately rounded products. It is used only by the LU sfmin
// test; iCamax uses RCabs1 instead.
#define CFACTOR_HYPOT(P,U,LIMB) \
    LIMB limbs[3][real_t<PB>::N]; P##mpfr_t x,y,z; \
    const real_t<PB> *v[2]={&a.re_.v,&a.im_.v}; P##mpfr_ptr refs[2]={x,y}; \
    for(int q=0;q<2;++q) { \
        for(int t=0;t<real_t<PB>::N;++t) limbs[q][t]=v[q]->m[t]; \
        P##mpfr_custom_init_set(refs[q],v[q]->is_zero()?U##MPFR_ZERO_KIND:(v[q]->sign>0?U##MPFR_REGULAR_KIND:-U##MPFR_REGULAR_KIND),v[q]->exp,PB,limbs[q]); \
    } \
    P##mpfr_custom_init_set(z,U##MPFR_ZERO_KIND,0,PB,limbs[2]); \
    P##mpfr_hypot(z,x,y,U##MPFR_RNDN); \
    int kind=P##mpfr_custom_get_kind(z); \
    if(kind!=U##MPFR_ZERO_KIND && kind!=U##MPFR_REGULAR_KIND) return false; \
    r=DC(); \
    if(kind==U##MPFR_REGULAR_KIND) { \
        r.re_.v.sign=1; r.re_.v.exp=P##mpfr_custom_get_exp(z); \
        const LIMB *p=(const LIMB*)P##mpfr_custom_get_significand(z); \
        for(int t=0;t<real_t<PB>::N;++t) r.re_.v.m[t]=p[t]; \
    } \
    return true;
template<class Tag> struct complex_factor_traits;
template<int PB> struct complex_factor_traits<factor_fixed_tag<PB>> {
    using DR=factor_real<PB>; using DC=factor_complex_value<PB>;
    static DC host_one(int sign,long) { return DC(sign); }
    L3_HD static bool fits(const DR &v) { return v.valid && (v.v.is_zero() || (v.v.exp>=-1073741823L && v.v.exp<=1073741823L)); }
    MPLAPACK_MPFR_DEV_OP static bool magnitude(DC &r,const DC &a) {
#if defined(MPLAPACK_MPFR_CDIV_HOST) && !defined(__CUDA_ARCH__)
        CFACTOR_HYPOT(,,mp_limb_t)
#elif defined(MPLAPACK_MPFR_CDIV_DEVICE) && defined(__CUDA_ARCH__)
        factor_arena_reset();
        CFACTOR_HYPOT(cu_,CU_,cu_mp_limb_t)
#else
        (void)r; (void)a; return false;
#endif
    }
    L3_HD static bool root(DC &r,const DC &a) {
        r=DC(); return cu_fp::cu_sqrt_rn<PB>(r.re_.v,a.re_.v)==0;
    }
};
#undef CFACTOR_HYPOT
template<int NL> struct complex_factor_traits<factor_runtime_tag<NL>> {
    using DR=rreal<NL>; using DC=factor_runtime_complex<NL>;
    static DC host_one(int sign,long prec) {
        DC v; v.re_.kind=sign*3; v.re_.exp=1;
        for(int i=0;i<NL;++i) v.re_.m[i]=0;
        v.re_.m[(prec+63)/64-1]=mplapack_rt_limb(1)<<63;
        return v;
    }
    L3_HD static bool fits(const DR &v) {
        int k=v.kind<0?-v.kind:v.kind;
        return k==2 || (k==3 && v.exp>=-1073741823L && v.exp<=1073741823L);
    }
    MPLAPACK_MPFR_DEV_OP static bool magnitude(DC &r,const DC &a) {
#ifdef MPLAPACK_RT_ACTIVE
        MPLAPACK_RT_DECL(x); MPLAPACK_RT_DECL(y); MPLAPACK_RT_DECL(z);
        MPLAPACK_RT_VIEW(x,a.re_); MPLAPACK_RT_VIEW(y,a.im_);
        DR t;
        MPLAPACK_RT(mpfr_custom_init_set)(z,MPLAPACK_RTU(MPFR_ZERO_KIND),0,rt_prec(),t.m);
        factor_arena_reset();
        MPLAPACK_RT(mpfr_hypot)(z,x,y,MPLAPACK_RTU(MPFR_RNDN));
        MPLAPACK_RT_TAKE(t,z); r=DC(t); return fits(t);
#else
        (void)r; (void)a; return false;
#endif
    }
    MPLAPACK_MPFR_DEV_OP static bool root(DC &r,const DC &a) {
#ifdef MPLAPACK_RT_ACTIVE
        MPLAPACK_RT_DECL(x); MPLAPACK_RT_DECL(z);
        MPLAPACK_RT_VIEW(x,a.re_); DR t;
        MPLAPACK_RT(mpfr_custom_init_set)(z,MPLAPACK_RTU(MPFR_ZERO_KIND),0,rt_prec(),t.m);
        factor_arena_reset();
        MPLAPACK_RT(mpfr_sqrt)(z,x,MPLAPACK_RTU(MPFR_RNDN));
        MPLAPACK_RT_TAKE(t,z); r=DC(t); return fits(t);
#else
        (void)r; (void)a; return false;
#endif
    }
};
template<class Tag> struct complex_factor_engine {
    using Traits=complex_factor_traits<Tag>;
    using DR=typename Traits::DR; using V=typename Traits::DC;
    struct M { V *p; long ld; };
    // The bounded physical-thread pool also provides arenas for the fixed
    // engine's MPC division. The two streams have disjoint arena indices.
    static constexpr bool runtime=true;
    using Update=complex_factor_engine;
    mutable bool ok=true;
    L3_HD Update update() const { return *this; }
    L3_HD void thread(long) { factor_arena_reset(); }
    L3_HD static M shift(M a,long n) { a.p+=n; return a; }
    L3_HD V reg(int) const { return V(); }
    L3_HD V at(const M &a,long i,long j) const { return a.p[i+j*a.ld]; }
    L3_HD void checked(const V &v) const { if(!Traits::fits(v.re_) || !Traits::fits(v.im_)) ok=false; }
    L3_HD void put(M &a,long i,long j,const V &v) const { checked(v); a.p[i+j*a.ld]=v; }
    L3_HD void set(V &r,const V &a) const { r=a; }
    L3_HD void zero(V &r) const { r=V(); }
    L3_HD void add(V &r,const V &a,const V &b) const { r=a+b; checked(r); }
    L3_HD void mul(V &r,const V &a,const V &b) const { r=a*b; checked(r); }
    L3_HD void div(V &r,const V &a,const V &b) const { if(iszero(b)) { ok=false; r=V(); } else { r=a/b; checked(r); } }
    L3_HD void absolute(V &r,const V &a) const { r=V(RCabs1(a)); checked(r); }
    L3_HD void magnitude(V &r,const V &a) const { if(!Traits::magnitude(r,a)) ok=false; checked(r); }
    L3_HD int compare(const V &a,const V &b) const { return a.re_<b.re_?-1:a.re_>b.re_?1:0; }
    L3_HD bool iszero(const V &a) const { return a==V(); }
    L3_HD void root(V &r,const V &a) const { if(!Traits::root(r,a)) ok=false; checked(r); }
    template<class X> L3_HD void trsm(const l3_shape &s,X &x,long lane) const {
        const char *side=s.left?"L":"R", *uplo=s.upper?"U":"L";
        const char *trans=s.trans?"C":"N", *diag=s.nounit?"N":"U";
        if(s.left) {
            V *b=x.c.p+lane*x.c.ld;
            Ctrsm_dev<DR,V>(side,uplo,trans,diag,s.m,1,V(1),x.a.p,x.a.ld,b,x.c.ld);
            for(long i=0;i<s.m;++i) checked(b[i]);
        } else {
            V *b=x.c.p+lane;
            Ctrsm_dev<DR,V>(side,uplo,trans,diag,1,s.n,V(1),x.a.p,x.a.ld,b,x.c.ld);
            for(long j=0;j<s.n;++j) checked(b[j*x.c.ld]);
        }
    }
    template<class X> L3_HD void herk(const l3_shape &s,X &x,long lane) const {
        if(lane%s.n) return;
        long j=lane/s.n;
        Cherk_columns_dev<DR,V>(s.upper?"U":"L",s.trans?"C":"N",s.n,s.k,DR(-1),x.a.p,x.a.ld,DR(1),x.c.p,x.c.ld,j,j+1);
        for(long i=s.upper?0:j;i<(s.upper?j+1:s.n);++i) checked(x.c.p[i+j*x.c.ld]);
    }
    L3_HD void gemm(const gemm_shape &s,M al,M be,M a,M b,M c,long i,long j) const {
        Cgemm_dev<DR,V>(s.nota?"N":"C",s.notb?"N":"C",1,1,s.k,al.p[0],
                       a.p+(s.nota?i:i*a.ld),a.ld,b.p+(s.notb?j*b.ld:j),b.ld,be.p[0],c.p+i+j*c.ld,c.ld);
        checked(c.p[i+j*c.ld]);
    }
};
#ifndef MPLAPACK_MPFR_CUDA_HOST_EMULATION
struct complex_arena_guard {
    char *base=mpc_cuda_arena_base; size_t slab=mpc_cuda_arena_slab; size_t *top=mpc_cuda_arena_top;
    long host_precision=mplapack_mpfr_rt_prec_h,device_precision=0;
    int error=0; bool changed=true,device_known=false;
    complex_arena_guard(char *b,size_t s,size_t *t,long prec) {
        error=cudaMemcpyFromSymbol(&device_precision,mplapack_mpfr_rt_prec_d,sizeof(long));
        device_known=!error;
        if(!error) error=cudaMemcpyToSymbol(mplapack_mpfr_rt_prec_d,&prec,sizeof(long));
        mpc_cuda_arena_base=b; mpc_cuda_arena_slab=s; mpc_cuda_arena_top=t; mplapack_mpfr_rt_prec_h=prec;
    }
    int restore() {
        if(!changed) return 0;
        mpc_cuda_arena_base=base; mpc_cuda_arena_slab=slab; mpc_cuda_arena_top=top;
        mplapack_mpfr_rt_prec_h=host_precision;
        int rc=device_known?cudaMemcpyToSymbol(mplapack_mpfr_rt_prec_d,&device_precision,sizeof(long)):0;
        if(!rc) changed=false;
        return rc;
    }
    ~complex_arena_guard() { (void)restore(); }
};
#endif
}
template<class Tag> int factor_complex(const factor_shape &s,long prec,typename cgpu_types<Tag>::DC *a,
    const typename cgpu_types<Tag>::DC &safe,long *piv,long &info,factor_times &times) {
    using E=complex_factor_engine<Tag>; using V=typename E::V;
    using Clock=std::chrono::steady_clock;
    size_t count=size_t(s.m)*s.n;
#ifndef MPLAPACK_MPFR_CUDA_HOST_EMULATION
    int ndev=0; if(cudaGetDeviceCount(&ndev)!=cudaSuccess || !ndev) return 1;
    // Complex MPC call chains need their own measured stack policy.
    bool fixed=std::is_same<Tag,factor_fixed_tag<256>>::value || std::is_same<Tag,factor_fixed_tag<512>>::value ||
               std::is_same<Tag,factor_fixed_tag<768>>::value || std::is_same<Tag,factor_fixed_tag<1024>>::value ||
               std::is_same<Tag,factor_fixed_tag<2048>>::value;
    factor_stack_guard stack(fixed?complex_fixed_stack_bytes(prec):complex_runtime_stack_bytes(prec)); if(stack.error) return stack.error;
    factor_buffer<char> arena; factor_buffer<size_t> tops;
    size_t slab=16384*((prec+1023)/1024),pool=2*factor_pool;
    if(arena.allocate(pool*slab) || tops.allocate(pool) || cudaMemset(tops.p,0,pool*sizeof(size_t))) return 1;
    complex_arena_guard guard(arena.p,slab,tops.p,prec); if(guard.error) return guard.error;
#else
    long previous_precision=mplapack_mpfr_rt_prec_h; mplapack_mpfr_rt_prec_h=prec;
    struct restore_precision { long p; ~restore_precision() { mplapack_mpfr_rt_prec_h=p; } } guard{previous_precision};
#endif
    factor_buffer<V> matrix,constants;
    if(matrix.allocate(count) || constants.allocate(4)) return 1;
    V c[4]={E::Traits::host_one(1,prec),E::Traits::host_one(-1,prec),safe,V()};
    if(constants.upload(c,4)) return 1;
    factor_context<E> ctx(E{},typename E::M{constants.p,4},s.m,s.n);
    if(ctx.error) return ctx.error;
    auto start=Clock::now();
    if(matrix.upload(a,count)) return 1;
#ifndef MPLAPACK_MPFR_CUDA_HOST_EMULATION
    if(cudaStreamSynchronize(nullptr)!=cudaSuccess) return 1;
#endif
    auto uploaded=Clock::now();
    long result=ctx.factor(s,typename E::M{matrix.p,s.m});
    if(ctx.error) return ctx.error;
    auto factored=Clock::now();
    if(matrix.download(a,count) || (!s.chol && ctx.pivots.download(piv,std::min(s.m,s.n)))) return 1;
    auto downloaded=Clock::now();
#ifndef MPLAPACK_MPFR_CUDA_HOST_EMULATION
    if(guard.restore() || stack.restore()) return 1;
#endif
    info=result;
    times.upload=std::chrono::duration<double>(uploaded-start).count();
    times.factor=std::chrono::duration<double>(factored-uploaded).count();
    times.download=std::chrono::duration<double>(downloaded-factored).count();
    return 0;
}
}
#endif
