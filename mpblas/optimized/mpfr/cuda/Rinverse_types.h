// Inverse-local adapters: reuse factor arithmetic, track every intermediate.
// No changes to the arithmetic or signed-zero conventions of existing BLAS.
#ifndef MPLAPACK_MPFR_INVERSE_TYPES_H
#define MPLAPACK_MPFR_INVERSE_TYPES_H
#include "Cfactor_runtime_types.h"
#include <type_traits>
namespace mplapack_mpfr_inverse __attribute__((visibility("hidden"))) {
using namespace mplapack_blas_dev;
// One sticky flag per physical thread, in the bounded factor launch pool.
// The common MPFR mutex excludes overlapping calls; kernels use one stream.
extern thread_local bool failed_host;
#if defined(__CUDACC__)
extern __device__ bool failed_device[8192];
#endif
MPLAPACK_DEV_HD bool &failure() {
#if defined(__CUDA_ARCH__)
    return failed_device[blockIdx.x*blockDim.x+threadIdx.x];
#else
    return failed_host;
#endif
}
template<int PB> MPLAPACK_DEV_HD bool fits(const factor_real<PB> &a) {
    return a.valid && (a.v.is_zero() || (a.v.exp>=-1073741823L && a.v.exp<=1073741823L));
}
template<int NL> MPLAPACK_DEV_HD bool fits(const rreal<NL> &a) {
    int k=a.kind<0?-a.kind:a.kind;
    return k==2 || (k==3 && a.exp>=-1073741823L && a.exp<=1073741823L);
}
// MPC's exponent/exception state is shared, not thread-local on the device.
// Avoid its underflow/overflow paths rather than racing on exception flags.
// This conservative per-operation guard also catches failures later overwritten
// by reference assignments (e.g. a complex diagonal's imaginary component).
template<int PB> MPLAPACK_DEV_HD bool safe(const factor_real<PB> &a) {
    const long bound=(1073741823L-4L*PB-128)/4;
    return fits(a) && (a.v.is_zero() || (a.v.exp>=-bound && a.v.exp<=bound));
}
template<int NL> MPLAPACK_DEV_HD bool safe(const rreal<NL> &a) {
    int k=a.kind<0?-a.kind:a.kind;
    const long bound=(1073741823L-4L*rt_prec()-128)/4;
    return fits(a) && (k==2 || (a.exp>=-bound && a.exp<=bound));
}
template<int PB> MPLAPACK_MPFR_DEV_OP factor_real<PB> divide(const factor_real<PB> &a,const factor_real<PB> &b) {
    mplapack_mpfr_cuda::l3_fixed_engine<PB,true> e;
    factor_real<PB> r; e.div(r.v,a.v,b.v); r.valid=a.valid && b.valid && e.ok; return r;
}
template<int NL> MPLAPACK_MPFR_DEV_OP rreal<NL> divide(const rreal<NL> &a,const rreal<NL> &b) {
    rreal<NL> r;
#ifdef MPLAPACK_RT_ACTIVE
    MPLAPACK_RT_DECL(x); MPLAPACK_RT_DECL(y); MPLAPACK_RT_DECL(z);
    MPLAPACK_RT_VIEW(x,a); MPLAPACK_RT_VIEW(y,b);
    MPLAPACK_RT(mpfr_custom_init_set)(z,MPLAPACK_RTU(MPFR_ZERO_KIND),0,rt_prec(),r.m);
    factor_arena_reset(); MPLAPACK_RT(mpfr_div)(z,x,y,MPLAPACK_RTU(MPFR_RNDN));
    MPLAPACK_RT_TAKE(r,z);
#endif
    return r;
}
template<int PB> MPLAPACK_DEV_HD void number(factor_real<PB> &r,double d) { r=factor_real<PB>(d); }
template<int NL> MPLAPACK_DEV_HD void number(rreal<NL> &r,double d) {
#ifdef MPLAPACK_RT_ACTIVE
    factor_arena_reset(); rt_set_d(r,d);
#else
    // The scheduler constructs only 0 and +/-1 on the host in CUDA builds.
    r.kind=d==0?2:d<0?-3:3; r.exp=d==0?0:1;
    for(int i=0;i<NL;++i) r.m[i]=0;
    if(d!=0) r.m[(rt_prec()+63)/64-1]=mplapack_rt_limb(1)<<63;
#endif
}
template<class Tag,bool Complex> struct value {
    using R=typename cgpu_types<Tag>::DR;
    using C=typename cgpu_types<Tag>::DC;
    using Base=typename std::conditional<Complex,C,R>::type;
    Base v;
    MPLAPACK_DEV_HD value() {}
    MPLAPACK_DEV_HD value(double d) { R r; number(r,d); if constexpr(Complex) v=C(r); else v=r; }
    MPLAPACK_DEV_HD explicit value(const Base &b):v(b) { check(); }
    template<bool B=Complex,typename std::enable_if<B,int>::type=0>
    MPLAPACK_DEV_HD value(const value<Tag,false> &r):v(r.v) {}
    template<bool B=Complex,typename std::enable_if<B,int>::type=0>
    MPLAPACK_DEV_HD value(const value<Tag,false> &r,const value<Tag,false> &i):v(r.v,i.v) {}
    MPLAPACK_DEV_HD void check() const {
        if constexpr(Complex) { if(!fits(v.re_) || !fits(v.im_)) failure()=true; }
        else if(!fits(v)) failure()=true;
    }
    MPLAPACK_DEV_HD bool safe_operand() const {
        if constexpr(Complex) return safe(v.re_) && safe(v.im_);
        else return safe(v);
    }
    MPLAPACK_DEV_HD value<Tag,false> real() const { return value<Tag,false>(v.re_); }
    MPLAPACK_DEV_HD value<Tag,false> imag() const { return value<Tag,false>(v.im_); }
    MPLAPACK_DEV_HD value operator-() const { return value(-v); }
    MPLAPACK_DEV_HD value operator+() const { return *this; }
    MPLAPACK_MPFR_DEV_OP static value op(const value &a,const value &b,int which) {
        if(!a.safe_operand() || !b.safe_operand()) { failure()=true; return value(); }
        factor_arena_reset();
        if(which==0) return value(a.v+b.v);
        if(which==1) return value(a.v-b.v);
        if(which==2) return value(a.v*b.v);
        if(b==value()) { failure()=true; return value(); }
        if constexpr(Complex) return value(a.v/b.v);
        else return value(divide(a.v,b.v));
    }
#define INVERSE_OP(OP,ASSIGN,N) \
    friend MPLAPACK_DEV_HD value operator OP(const value &a,const value &b) { return op(a,b,N); } \
    MPLAPACK_DEV_HD value &operator ASSIGN(const value &b) { return *this=op(*this,b,N); }
    INVERSE_OP(+,+=,0)
    INVERSE_OP(-,-=,1)
    INVERSE_OP(*,*=,2)
    INVERSE_OP(/,/=,3)
#undef INVERSE_OP
    friend MPLAPACK_DEV_HD bool operator==(const value &a,const value &b) { a.check();b.check();return a.v==b.v; }
    friend MPLAPACK_DEV_HD bool operator!=(const value &a,const value &b) { return !(a==b); }
    friend MPLAPACK_DEV_HD bool operator<(const value &a,const value &b) { return a.v<b.v; }
    friend MPLAPACK_DEV_HD bool operator>(const value &a,const value &b) { return b<a; }
};
template<class Tag> MPLAPACK_DEV_HD value<Tag,false> abs(const value<Tag,false> &a) { return value<Tag,false>(abs(a.v)); }
template<class Tag> MPLAPACK_DEV_HD value<Tag,true> conj(const value<Tag,true> &a) { return value<Tag,true>(conj(a.v)); }
}
#endif
