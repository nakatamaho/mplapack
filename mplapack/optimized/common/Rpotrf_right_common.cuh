#ifndef MPLAPACK_BF_POTRF_RIGHT_COMMON_CUH
#define MPLAPACK_BF_POTRF_RIGHT_COMMON_CUH
#include "../../../mpblas/optimized/common/Rpotrf_right_common.h"
#include <type_traits>
namespace mplapack_lapack_gpu {
template<class Ops,bool Complex> struct bf_right_policy {
    using R=dreal<Ops>;
    using V=typename std::conditional<Complex,dcomplex<R>,R>::type;
    struct M { V *p=nullptr; long ld=0; };
    MPLAPACK_DEV_HD V reg(int) const { return V(); }
    MPLAPACK_DEV_HD V at(const M &a,long i,long j) const { return a.p[i+j*a.ld]; }
    MPLAPACK_DEV_HD bool iszero(const V &v) const { return v==V(0); }
    MPLAPACK_DEV_HD void init(V &v,const V &a,bool diag) const {
        if constexpr(Complex) v=diag?V(a.real()):a;
        else v=a;
    }
    MPLAPACK_DEV_HD void scaled(V &v,const V &b,bool rank,const M &) const {
        if constexpr(Complex) v=rank?(-R(1))*conj(b):(-V(1))*conj(b);
        else v=(-R(1))*b;
    }
    MPLAPACK_DEV_HD void product(V &v,const V &a,const V &b) const { v=a*b; }
    MPLAPACK_DEV_HD void accumulate(V &v,const V &term,bool diag) const {
        if constexpr(Complex) { if(diag) v=V(v.real()+term.real()); else v+=term; }
        else v+=term;
    }
    MPLAPACK_DEV_HD void put(M &a,long i,long j,const V &v) const { a.p[i+j*a.ld]=v; }
    MPLAPACK_DEV_HD bool ok() const { return true; }
};
template<class Ops,bool Complex,class Context> static void bf_right_launch(Context &ctx,long n,long nb,long col,long k,long first,long cols,typename bf_right_policy<Ops,Complex>::V *a,long ld) {
    if(ctx.error || cols<=0) return;
    right_chol_job<bf_right_policy<Ops,Complex>> j{};
    j.a={a,ld}; j.n=n; j.nb=nb; j.col=col; j.k=k; j.first=first; j.cols=cols;
#if defined(__CUDACC__)
    ctx.error=right_chol_device(j,ctx.stream());
#else
    ctx.error=right_chol_host(j);
#endif
}
}
#endif
