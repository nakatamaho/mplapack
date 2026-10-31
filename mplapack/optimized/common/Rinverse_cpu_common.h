// Copyright (c) 2026 Nakata, Maho. BSD-2-Clause (see Rgetrf.cpp).
// Reference-order blocked inverses with recursive independent-output tasks.
#ifndef MPLAPACK_INVERSE_CPU_COMMON_H
#define MPLAPACK_INVERSE_CPU_COMMON_H
#include "Rlookahead_common.h"
namespace mplapack_inverse {
template<class T> struct ops;
#define MPLAPACK_INVERSE_OPS(T,P,TRANS,RANK) \
template<> struct ops<T> { \
    void tri(bool u,const char *d,INTEGER n,T *a,INTEGER ld,INTEGER &info) const { P##trti2(u?"U":"L",d,n,a,ld,info); } \
    void lauum(bool u,INTEGER n,T *a,INTEGER ld,INTEGER &info) const { P##lauu2(u?"U":"L",n,a,ld,info); } \
    void mm(bool left,bool upper,const char *tr,const char *d,INTEGER m,INTEGER n,T alpha,T *a,T *b,INTEGER ld) const { P##trmm(left?"L":"R",upper?"U":"L",tr,d,m,n,alpha,a,ld,b,ld); } \
    void sm(bool upper,const char *d,INTEGER m,INTEGER n,T *a,T *b,INTEGER ld) const { const T one=T(1); P##trsm("R",upper?"U":"L","N",d,m,n,-one,a,ld,b,ld); } \
    void rank(bool upper,INTEGER n,INTEGER k,T *a,T *c,INTEGER ld) const { RANK(upper?"U":"L",upper?"N":TRANS,n,k,REAL(1),a,ld,REAL(1),c,ld); } \
    void update(bool upper,INTEGER m,INTEGER n,INTEGER k,T *a,T *b,T *c,INTEGER ld) const { P##gemm(upper?"N":TRANS,upper?TRANS:"N",m,n,k,T(1),a,ld,b,ld,T(1),c,ld); } \
    const char *trans() const { return TRANS; } \
};
MPLAPACK_INVERSE_OPS(REAL,R,"T",Rsyrk)
MPLAPACK_INVERSE_OPS(COMPLEX,C,"C",Cherk)
#undef MPLAPACK_INVERSE_OPS

// Split independent output lanes only. In-place triangular multiplication
// may split columns on the left, rows on the right, never its dependent axis.
template<class F> static void ranges(INTEGER lo,INTEGER hi,bool tasks,const F &f) {
    if(!tasks || hi-lo<=32) { if(hi>lo) f(lo,hi); return; }
    INTEGER mid=lo+(hi-lo)/2;
#pragma omp taskgroup
    {
#pragma omp task shared(f) firstprivate(lo,mid,tasks)
        ranges(lo,mid,tasks,f);
        ranges(mid,hi,tasks,f);
    }
}
template<class F> static void team(INTEGER n,INTEGER nb,bool fast,const F &f) {
#ifdef _OPENMP
    if(!fast && mplapack_lapack::la_lookahead(true,n,nb)) {
        const mplapack_blas::blas_thread_context ctx;
        const int threads=int(std::min(long(omp_get_max_threads()),std::max(1L,(long(n)+31)/32)));
#pragma omp parallel num_threads(threads)
        {
            const mplapack_blas::blas_thread_context::scope scope(ctx);
#pragma omp single
            f(true);
        }
        return;
    }
#endif
    f(false);
}
template<class T> static void tri_block(bool upper,const char *diag,INTEGER before,INTEGER jb,
    T *inverse,T *original,T *off,INTEGER ld,bool tasks,const ops<T> &op) {
    ranges(0,jb,tasks,[&](INTEGER lo,INTEGER hi) {
        op.mm(true,upper,"N",diag,before,hi-lo,T(1),inverse,off+lo*ld,ld);
    });
    ranges(0,before,tasks,[&](INTEGER lo,INTEGER hi) {
        op.sm(upper,diag,hi-lo,jb,original,off+lo,ld);
    });
}
template<class T> static void tri_fast(bool upper,const char *diag,INTEGER n,T *a,INTEGER ld,INTEGER &info,const ops<T> &op) {
    if(n<=32) { op.tri(upper,diag,n,a,ld,info); return; }
    const INTEGER n1=n/2,n2=n-n1;
    T *a22=a+n1+n1*ld,*off=a+(upper?n1*ld:n1);
    if(upper) {
        tri_fast(true,diag,n1,a,ld,info,op);
        tri_block(true,diag,n1,n2,a,a22,off,ld,false,op);
        tri_fast(true,diag,n2,a22,ld,info,op);
    } else {
        tri_fast(false,diag,n2,a22,ld,info,op);
        tri_block(false,diag,n2,n1,a22,a,off,ld,false,op);
        tri_fast(false,diag,n1,a,ld,info,op);
    }
}
template<class T> static void tri(bool upper,const char *diag,INTEGER n,T *a,INTEGER ld,INTEGER nb,INTEGER &info) {
    info=0;
    if(Mlsame(diag,"N")) for(INTEGER i=0;i<n;++i) if(a[i+i*ld]==T(0)) { info=i+1; return; }
    const ops<T> op;
    const bool fast=mplapack_lapack::la_fast();
    if(fast) { tri_fast(upper,diag,n,a,ld,info,op); return; }
    if(nb<=1 || nb>=n) { op.tri(upper,diag,n,a,ld,info); return; }
    team(n,nb,false,[&](bool tasks) {
        for(INTEGER step=0;step<(n+nb-1)/nb;++step) {
            const INTEGER j=upper?step*nb:((n-1)/nb-step)*nb;
            const INTEGER jb=std::min(nb,n-j),before=upper?j:n-j-jb;
            T *d=a+j+j*ld;
            if(before) tri_block(upper,diag,before,jb,upper?a:a+j+jb+(j+jb)*ld,d,
                upper?a+j*ld:a+j+jb+j*ld,ld,tasks,op);
            op.tri(upper,diag,jb,d,ld,info);
        }
    });
}
template<class T> static void lauum_fast(bool upper,INTEGER n,T *a,INTEGER ld,INTEGER &info,const ops<T> &op) {
    if(n<=32) { op.lauum(upper,n,a,ld,info); return; }
    const INTEGER n1=n/2,n2=n-n1;
    T *a22=a+n1+n1*ld,*off=a+(upper?n1*ld:n1);
    lauum_fast(upper,n1,a,ld,info,op);
    op.rank(upper,n1,n2,off,a,ld);
    op.mm(!upper,upper,op.trans(),"N",upper?n1:n2,upper?n2:n1,T(1),a22,off,ld);
    lauum_fast(upper,n2,a22,ld,info,op);
}
template<class T> static void lauum(bool upper,INTEGER n,T *a,INTEGER ld,INTEGER nb,INTEGER &info) {
    info=0; const ops<T> op;
    if(mplapack_lapack::la_fast()) { lauum_fast(upper,n,a,ld,info,op); return; }
    if(nb<=1 || nb>=n) { op.lauum(upper,n,a,ld,info); return; }
    team(n,nb,false,[&](bool tasks) {
        for(INTEGER j=0;j<n;j+=nb) {
            const INTEGER jb=std::min(nb,n-j),rest=n-j-jb;
            T *d=a+j+j*ld,*off=upper?a+j*ld:a+j;
            ranges(0,j,tasks,[&](INTEGER lo,INTEGER hi) {
                op.mm(!upper,upper,op.trans(),"N",upper?hi-lo:jb,upper?jb:hi-lo,T(1),d,
                    off+(upper?lo:lo*ld),ld);
            });
            op.lauum(upper,jb,d,ld,info);
            if(rest) {
                // Each output still receives the full reference reduction.
                ranges(0,j,tasks,[&](INTEGER lo,INTEGER hi) {
                    op.update(upper,upper?hi-lo:jb,upper?jb:hi-lo,rest,
                        upper?a+lo+(j+jb)*ld:a+j+jb+j*ld,
                        upper?a+j+(j+jb)*ld:a+j+jb+lo*ld,
                        off+(upper?lo:lo*ld),ld);
                });
                op.rank(upper,jb,rest,upper?a+j+(j+jb)*ld:a+j+jb+j*ld,d,ld);
            }
        }
    });
}
}
#endif
