// Factor-local runtime MPC values. Reuse the BLAS arithmetic, reclaiming
// its per-operation scratch because factor jobs contain many BLAS calls.
#ifndef MPLAPACK_MPFR_CFACTOR_RUNTIME_TYPES_H
#define MPLAPACK_MPFR_CFACTOR_RUNTIME_TYPES_H
#include "Cfactor_fixed_types.h"
namespace mplapack_blas_dev {
template<int NL> struct factor_runtime_tag {};
template<int NL> struct factor_runtime_complex : rcomplex<NL> {
    using R=rreal<NL>; using Base=rcomplex<NL>;
    MPLAPACK_DEV_HD factor_runtime_complex() {}
    MPLAPACK_DEV_HD factor_runtime_complex(double d) { factor_arena_reset(); rt_set_d(this->re_,d); }
    MPLAPACK_DEV_HD factor_runtime_complex(double r,double i) {
        factor_arena_reset(); rt_set_d(this->re_,r);
        factor_arena_reset(); rt_set_d(this->im_,i);
    }
    MPLAPACK_DEV_HD factor_runtime_complex(const R &r):Base(r) {}
    MPLAPACK_DEV_HD factor_runtime_complex(const R &r,const R &i):Base(r,i) {}
    MPLAPACK_DEV_HD factor_runtime_complex operator+() const { return *this; }
    MPLAPACK_DEV_HD factor_runtime_complex operator-() const { return {-this->re_,-this->im_}; }
    MPLAPACK_DEV_HD factor_runtime_complex &operator+=(const factor_runtime_complex &r);
    MPLAPACK_DEV_HD factor_runtime_complex &operator-=(const factor_runtime_complex &r);
    MPLAPACK_DEV_HD factor_runtime_complex &operator*=(const factor_runtime_complex &r);
    MPLAPACK_DEV_HD factor_runtime_complex &operator/=(const factor_runtime_complex &r);
    MPLAPACK_DEV_HD factor_runtime_complex &operator+=(const R &r) { return *this+=factor_runtime_complex(r); }
    MPLAPACK_DEV_HD factor_runtime_complex &operator-=(const R &r) { return *this-=factor_runtime_complex(r); }
    MPLAPACK_DEV_HD factor_runtime_complex &operator*=(const R &r) { return *this*=factor_runtime_complex(r); }
};
template<int NL> MPLAPACK_MPFR_DEV_OP factor_runtime_complex<NL> factor_runtime_op(
    const factor_runtime_complex<NL> &a,const factor_runtime_complex<NL> &b,int op) {
    factor_runtime_complex<NL> r;
    factor_arena_reset(); rt_cop<NL>(r,a,b,op); return r;
}
#define CFACTOR_RT_OP(OP,ASSIGN,CODE) \
template<int NL> MPLAPACK_DEV_HD factor_runtime_complex<NL> operator OP(const factor_runtime_complex<NL> &a,const factor_runtime_complex<NL> &b) { return factor_runtime_op(a,b,CODE); } \
template<int NL> MPLAPACK_DEV_HD factor_runtime_complex<NL> &factor_runtime_complex<NL>::operator ASSIGN(const factor_runtime_complex &r) { return *this=*this OP r; } \
template<int NL> MPLAPACK_DEV_HD factor_runtime_complex<NL> operator OP(const factor_runtime_complex<NL> &a,const rreal<NL> &b) { return a OP factor_runtime_complex<NL>(b); } \
template<int NL> MPLAPACK_DEV_HD factor_runtime_complex<NL> operator OP(const rreal<NL> &a,const factor_runtime_complex<NL> &b) { return factor_runtime_complex<NL>(a) OP b; }
CFACTOR_RT_OP(+,+=,0)
CFACTOR_RT_OP(-,-=,1)
CFACTOR_RT_OP(*,*=,2)
CFACTOR_RT_OP(/,/=,3)
#undef CFACTOR_RT_OP
template<int NL> MPLAPACK_DEV_HD bool operator==(const factor_runtime_complex<NL> &a,const factor_runtime_complex<NL> &b) { return rt_ccmp<NL>(a,b)==0; }
template<int NL> MPLAPACK_DEV_HD bool operator!=(const factor_runtime_complex<NL> &a,const factor_runtime_complex<NL> &b) { return !(a==b); }
template<int NL> MPLAPACK_DEV_HD factor_runtime_complex<NL> conj(const factor_runtime_complex<NL> &a) { return {a.re_,-a.im_}; }
template<int NL> MPLAPACK_DEV_HD rreal<NL> RCabs1(const factor_runtime_complex<NL> &a) {
    factor_arena_reset(); return abs(a.re_)+abs(a.im_);
}
template<int NL> struct cgpu_types<factor_runtime_tag<NL>> { using DR=rreal<NL>; using DC=factor_runtime_complex<NL>; };
}
#endif
