// Tiny recursive Cholesky leaves, without host feedback between scalar steps.
#ifndef MPLAPACK_POTRF_SMALL_COMMON_CUH
#define MPLAPACK_POTRF_SMALL_COMMON_CUH
#include <type_traits>
namespace mplapack_lapack_gpu {
template<class Ops> struct complex_chol_base;
template<class Ops,class E,int Limit> MPLAPACK_DEV_HD long small_chol_factor(bool upper,long n,E *a,long ld) {
    using R=dreal<Ops>;
    if(!n) return 0;
    if(n==1) {
        long info=0;
        if constexpr(std::is_same<E,R>::value) base_chol<Ops>{a,&info}(0);
        else complex_chol_base<Ops>{a,&info}(0);
        return info;
    }
    if constexpr(Limit>1) {
        const long n1=n/2,n2=n-n1;
        long info=small_chol_factor<Ops,E,Limit/2>(upper,n1,a,ld);
        if(info) return info;
        E *off=upper?a+n1*ld:a+n1,*diag=a+n1+n1*ld;
        if constexpr(std::is_same<E,R>::value) {
            using F=typename Ops::F;
            F scalars[2]={Ops::from_d(1),Ops::from_d(-1)};
            mplapack_blas::gpu_lane_args<Ops> g{};
            g.s.op=mplapack_blas::L3_TRSM; g.s.left=upper; g.s.upper=upper;
            g.s.trans=1; g.s.nounit=1; g.s.alpha_one=1;
            g.s.m=upper?n1:n2; g.s.n=upper?n2:n1;
            g.a=reinterpret_cast<F*>(a); g.c=reinterpret_cast<F*>(off);
            g.al=scalars; g.lda=ld; g.ldc=ld;
            for(long u=0;u<mplapack_blas::blas_lanes(g.s);++u) mplapack_blas::gpu_lane(g,u);
            g.s.op=mplapack_blas::L3_SYRK; g.s.trans=upper; g.s.n=n2; g.s.k=n1;
            g.s.alpha_one=0; g.s.beta_one=1;
            g.a=reinterpret_cast<F*>(off); g.c=reinterpret_cast<F*>(diag);
            g.al=scalars+1; g.be=scalars;
            for(long u=0;u<mplapack_blas::blas_lanes(g.s);++u) mplapack_blas::gpu_lane(g,u);
        } else {
            job_Ctrsm<Ops> solve{};
            solve.side[0]=upper?'L':'R'; solve.uplo[0]=upper?'U':'L'; solve.ta[0]='C'; solve.diag[0]='N';
            solve.m=upper?n1:n2; solve.n=upper?n2:n1; solve.lda=ld; solve.ldb=ld;
            solve.alpha=E(1); solve.A=a; solve.B=off;
            for(long u=0;u<solve.units();++u) solve(u);
            job_Cherk<Ops> rank{};
            rank.uplo[0]=upper?'U':'L'; rank.tr[0]=upper?'C':'N';
            rank.n=n2; rank.k=n1; rank.lda=ld; rank.ldc=ld;
            rank.alpha=R(-1); rank.beta=R(1); rank.A=off; rank.C=diag;
            for(long u=0;u<rank.units();++u) rank(u);
        }
        info=small_chol_factor<Ops,E,Limit/2>(upper,n2,diag,ld);
        return info?info+n1:0;
    }
    return 0; // callers bound n by Limit
}
template<class Ops,class E> struct small_chol_job {
    bool upper; long n; E *a; long ld; long *info;
    MPLAPACK_DEV_HD void operator()(long) const { *info=small_chol_factor<Ops,E,4>(upper,n,a,ld); }
};
}
#endif
