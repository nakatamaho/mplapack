// Signed-zero, range-checked cu_freal policy for resident complex factors.
// BLAS bodies remain the existing generated reference templates. This does
// not alter the documented +0 convention of the public complex BLAS.
#ifndef MPLAPACK_MPFR_CFACTOR_FIXED_TYPES_H
#define MPLAPACK_MPFR_CFACTOR_FIXED_TYPES_H
#include "Cblas_mpfr_dev_types.h"
#include "blas3_kernel_cuda.h"
namespace mplapack_blas_dev {
// Every factor value owns its limbs (including runtime values). No operand
// survives in the bump arena, so reclaim scratch before each MPFR/MPC call.
MPLAPACK_DEV_HD void factor_arena_reset() {
#if defined(__CUDA_ARCH__)
    mpc_cuda_arena_reset();
#endif
}
template<int PB> struct factor_fixed_tag {};
template<int PB> struct factor_real {
    using F=cu_fp::cu_freal<PB>;
    F v; bool valid=true;
    MPLAPACK_DEV_HD factor_real() { v.set_zero(); v.sign=1; }
    MPLAPACK_DEV_HD factor_real(double d):v(F::from_double(d)) {}
    MPLAPACK_DEV_HD explicit factor_real(const F &x):v(x) {}
    MPLAPACK_DEV_HD factor_real operator-() const { auto r=*this; r.v.sign=-r.v.sign; return r; }
    MPLAPACK_DEV_HD factor_real operator+() const { return *this; }
    MPLAPACK_DEV_HD factor_real &operator+=(const factor_real &r);
    MPLAPACK_DEV_HD factor_real &operator-=(const factor_real &r);
    MPLAPACK_DEV_HD factor_real &operator*=(const factor_real &r);
};
template<int PB> MPLAPACK_MPFR_DEV_OP factor_real<PB> factor_real_op(const factor_real<PB> &a,const factor_real<PB> &b,int op) {
    mplapack_mpfr_cuda::l3_fixed_engine<PB,true> e;
    factor_real<PB> r;
    if(op==0) e.add(r.v,a.v,b.v); else if(op==1) e.sub(r.v,a.v,b.v); else e.mul(r.v,a.v,b.v);
    r.valid=a.valid && b.valid && e.ok; return r;
}
template<int PB> MPLAPACK_DEV_HD factor_real<PB> operator+(const factor_real<PB> &a,const factor_real<PB> &b) { return factor_real_op(a,b,0); }
template<int PB> MPLAPACK_DEV_HD factor_real<PB> operator-(const factor_real<PB> &a,const factor_real<PB> &b) { return factor_real_op(a,b,1); }
template<int PB> MPLAPACK_DEV_HD factor_real<PB> operator*(const factor_real<PB> &a,const factor_real<PB> &b) { return factor_real_op(a,b,2); }
template<int PB> MPLAPACK_DEV_HD factor_real<PB> &factor_real<PB>::operator+=(const factor_real &r) { return *this=*this+r; }
template<int PB> MPLAPACK_DEV_HD factor_real<PB> &factor_real<PB>::operator-=(const factor_real &r) { return *this=*this-r; }
template<int PB> MPLAPACK_DEV_HD factor_real<PB> &factor_real<PB>::operator*=(const factor_real &r) { return *this=*this*r; }
template<int PB> MPLAPACK_DEV_HD bool operator==(const factor_real<PB> &a,const factor_real<PB> &b) { return mp_cmp<PB>(a.v,b.v)==0; }
template<int PB> MPLAPACK_DEV_HD bool operator!=(const factor_real<PB> &a,const factor_real<PB> &b) { return !(a==b); }
template<int PB> MPLAPACK_DEV_HD bool operator<(const factor_real<PB> &a,const factor_real<PB> &b) { return mp_cmp<PB>(a.v,b.v)<0; }
template<int PB> MPLAPACK_DEV_HD bool operator>(const factor_real<PB> &a,const factor_real<PB> &b) { return b<a; }
template<int PB> MPLAPACK_DEV_HD factor_real<PB> abs(const factor_real<PB> &a) { auto r=a; r.v.sign=1; return r; }
template<int PB> struct factor_complex_value {
    using R=factor_real<PB>;
    R re_,im_;
    MPLAPACK_DEV_HD factor_complex_value() {}
    MPLAPACK_DEV_HD factor_complex_value(double d):re_(d) {}
    MPLAPACK_DEV_HD factor_complex_value(double r,double i):re_(r),im_(i) {}
    MPLAPACK_DEV_HD factor_complex_value(const R &r):re_(r) {}
    MPLAPACK_DEV_HD factor_complex_value(const R &r,const R &i):re_(r),im_(i) {}
    MPLAPACK_DEV_HD R &real() { return re_; }
    MPLAPACK_DEV_HD const R &real() const { return re_; }
    MPLAPACK_DEV_HD R &imag() { return im_; }
    MPLAPACK_DEV_HD const R &imag() const { return im_; }
    MPLAPACK_DEV_HD factor_complex_value operator-() const { return {-re_,-im_}; }
    MPLAPACK_DEV_HD factor_complex_value operator+() const { return *this; }
    MPLAPACK_DEV_HD factor_complex_value &operator+=(const factor_complex_value &r) { re_+=r.re_; im_+=r.im_; return *this; }
    MPLAPACK_DEV_HD factor_complex_value &operator-=(const factor_complex_value &r) { re_-=r.re_; im_-=r.im_; return *this; }
    MPLAPACK_DEV_HD factor_complex_value &operator*=(const factor_complex_value &r);
    MPLAPACK_DEV_HD factor_complex_value &operator/=(const factor_complex_value &r);
    MPLAPACK_DEV_HD factor_complex_value &operator+=(const R &r) { return *this+=factor_complex_value(r); }
    MPLAPACK_DEV_HD factor_complex_value &operator-=(const R &r) { return *this-=factor_complex_value(r); }
    MPLAPACK_DEV_HD factor_complex_value &operator*=(const R &r) { return *this*=factor_complex_value(r); }
};
template<int PB> MPLAPACK_DEV_HD bool operator==(const factor_complex_value<PB> &a,const factor_complex_value<PB> &b) { return a.re_==b.re_ && a.im_==b.im_; }
template<int PB> MPLAPACK_DEV_HD bool operator!=(const factor_complex_value<PB> &a,const factor_complex_value<PB> &b) { return !(a==b); }
template<int PB> MPLAPACK_DEV_HD factor_complex_value<PB> operator+(factor_complex_value<PB> a,const factor_complex_value<PB> &b) { return a+=b; }
template<int PB> MPLAPACK_DEV_HD factor_complex_value<PB> operator-(factor_complex_value<PB> a,const factor_complex_value<PB> &b) { return a-=b; }
template<int PB> MPLAPACK_MPFR_DEV_OP factor_complex_value<PB> operator*(const factor_complex_value<PB> &a,const factor_complex_value<PB> &b) {
    using C=factor_complex_value<PB>; using R=factor_real<PB>;
    C r{R(mp_cmul_re<PB>(a.re_.v,a.im_.v,b.re_.v,b.im_.v)),R(mp_cmul_im<PB>(a.re_.v,a.im_.v,b.re_.v,b.im_.v))};
    // FMMS/FMMA exact cancellations are +0 in RNDN. A negative zero only
    // occurs when both exact products are zero with the indicated signs.
    if(r.re_.v.is_zero()) r.re_.v.sign=(a.re_.v.is_zero() || b.re_.v.is_zero()) &&
        (a.im_.v.is_zero() || b.im_.v.is_zero()) && a.re_.v.sign!=b.re_.v.sign && a.im_.v.sign==b.im_.v.sign?-1:1;
    if(r.im_.v.is_zero()) r.im_.v.sign=(a.re_.v.is_zero() || b.im_.v.is_zero()) &&
        (a.im_.v.is_zero() || b.re_.v.is_zero()) && a.re_.v.sign!=b.im_.v.sign && a.im_.v.sign!=b.re_.v.sign?-1:1;
    bool valid=a.re_.valid && a.im_.valid && b.re_.valid && b.im_.valid;
    for(auto *v:{&r.re_,&r.im_}) v->valid=valid && (v->v.is_zero() || (v->v.exp>=-1073741823L && v->v.exp<=1073741823L));
    return r;
}
// Same MPFR custom views and MPC division as Ctrsm, retaining zero kinds
// rather than canonicalizing them; nonfinite quotients invalidate the call.
#define CFACTOR_DIV(P,U,LIMB) \
    LIMB limbs[6][cu_fp::cu_freal<PB>::N]; P##mpc_t x,y,z; \
    const factor_real<PB> *in[4]={&a.re_,&a.im_,&b.re_,&b.im_}; \
    P##mpfr_ptr refs[4]={P##mpc_realref(x),P##mpc_imagref(x),P##mpc_realref(y),P##mpc_imagref(y)}; \
    for(int q=0;q<4;++q) { \
        for(int t=0;t<cu_fp::cu_freal<PB>::N;++t) limbs[q][t]=in[q]->v.m[t]; \
        int kind=in[q]->v.is_zero()?U##MPFR_ZERO_KIND:U##MPFR_REGULAR_KIND; \
        P##mpfr_custom_init_set(refs[q],in[q]->v.sign*kind,in[q]->v.exp,PB,limbs[q]); \
    } \
    P##mpfr_custom_init_set(P##mpc_realref(z),U##MPFR_ZERO_KIND,0,PB,limbs[4]); \
    P##mpfr_custom_init_set(P##mpc_imagref(z),U##MPFR_ZERO_KIND,0,PB,limbs[5]); \
    P##mpc_div(z,x,y,U##MPC_RNDNN); \
    factor_real<PB> *out[2]={&r.re_,&r.im_}; P##mpfr_ptr os[2]={P##mpc_realref(z),P##mpc_imagref(z)}; \
    for(int q=0;q<2;++q) { \
        int kind=P##mpfr_custom_get_kind(os[q]); out[q]->v.set_zero(); out[q]->v.sign=kind<0?-1:1; \
        out[q]->valid=a.re_.valid && a.im_.valid && b.re_.valid && b.im_.valid; \
        if(kind==U##MPFR_REGULAR_KIND || kind==-U##MPFR_REGULAR_KIND) { \
            out[q]->v.exp=P##mpfr_custom_get_exp(os[q]); const LIMB *p=(const LIMB*)P##mpfr_custom_get_significand(os[q]); \
            for(int t=0;t<cu_fp::cu_freal<PB>::N;++t) out[q]->v.m[t]=p[t]; \
            out[q]->valid=out[q]->valid && out[q]->v.exp>=-1073741823L && out[q]->v.exp<=1073741823L; \
        } else if(kind!=U##MPFR_ZERO_KIND && kind!=-U##MPFR_ZERO_KIND) out[q]->valid=false; \
    }
template<int PB> MPLAPACK_MPFR_DEV_OP factor_complex_value<PB> operator/(const factor_complex_value<PB> &a,const factor_complex_value<PB> &b) {
    factor_complex_value<PB> r;
#if defined(MPLAPACK_MPFR_CDIV_HOST) && !defined(__CUDA_ARCH__)
    CFACTOR_DIV(,,mp_limb_t)
#elif defined(MPLAPACK_MPFR_CDIV_DEVICE) && defined(__CUDA_ARCH__)
    factor_arena_reset();
    CFACTOR_DIV(cu_,CU_,cu_mp_limb_t)
#else
    r.re_.valid=r.im_.valid=false;
#endif
    return r;
}
#undef CFACTOR_DIV
template<int PB> MPLAPACK_DEV_HD factor_complex_value<PB> &factor_complex_value<PB>::operator*=(const factor_complex_value &r) { return *this=*this*r; }
template<int PB> MPLAPACK_DEV_HD factor_complex_value<PB> &factor_complex_value<PB>::operator/=(const factor_complex_value &r) { return *this=*this/r; }
template<int PB> MPLAPACK_DEV_HD factor_complex_value<PB> operator+(factor_complex_value<PB> a,const factor_real<PB> &b) { return a+=b; }
template<int PB> MPLAPACK_DEV_HD factor_complex_value<PB> operator-(factor_complex_value<PB> a,const factor_real<PB> &b) { return a-=b; }
template<int PB> MPLAPACK_DEV_HD factor_complex_value<PB> operator*(factor_complex_value<PB> a,const factor_real<PB> &b) { return a*=b; }
template<int PB> MPLAPACK_DEV_HD factor_complex_value<PB> operator*(const factor_real<PB> &a,const factor_complex_value<PB> &b) { return factor_complex_value<PB>(a)*b; }
template<int PB> MPLAPACK_DEV_HD factor_complex_value<PB> conj(const factor_complex_value<PB> &a) { return {a.re_,-a.im_}; }
template<int PB> MPLAPACK_DEV_HD factor_real<PB> RCabs1(const factor_complex_value<PB> &a) { return abs(a.re_)+abs(a.im_); }
template<int PB> struct cgpu_types<factor_fixed_tag<PB>> { using DR=factor_real<PB>; using DC=factor_complex_value<PB>; };
}
#endif
