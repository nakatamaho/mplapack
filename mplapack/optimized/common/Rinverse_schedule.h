// Copyright (c) 2026 Nakata, Maho. BSD-2-Clause (see Rgetrf.cpp).
// Backend-independent reference-order inverse jobs and schedule.
#ifndef MPLAPACK_INVERSE_SCHEDULE_H
#define MPLAPACK_INVERSE_SCHEDULE_H
#include <algorithm>
#include <type_traits>
#include "inverse_dev/inverse_dev.h"
namespace mplapack_inverse_gpu {
using namespace mplapack_blas_dev;
template<class R,class C,bool Complex> struct inverse_arithmetic;
#define INVERSE_DEVICE_ARITHMETIC(FLAG,TYPE,P,TRANS,RANK) \
template<class R,class C> struct inverse_arithmetic<R,C,FLAG> { \
    using E=TYPE; \
    static constexpr bool cx=!std::is_same<R,C>::value; \
    MPLAPACK_DEV_HD static void tri(bool u,bool unit,long n,E *a,long ld,long &info) { P##trti2_dev<R,C>(u?"U":"L",unit?"U":"N",n,a,ld,info); } \
    MPLAPACK_DEV_HD static void lauum(bool u,long n,E *a,long ld,long &info) { P##lauu2_dev<R,C>(u?"U":"L",n,a,ld,info); } \
    MPLAPACK_DEV_HD static void mm(bool left,bool u,char tr,bool unit,long m,long n,E alpha,E *a,E *b,long ld) { char trans[2]={tr,0}; P##trmm_dev<R,C>(left?"L":"R",u?"U":"L",trans,unit?"U":"N",m,n,alpha,a,ld,b,ld); } \
    MPLAPACK_DEV_HD static void sm(bool u,bool unit,long m,long n,E alpha,E *a,E *b,long ld) { P##trsm_dev<R,C>("R",u?"U":"L","N",unit?"U":"N",m,n,alpha,a,ld,b,ld); } \
    MPLAPACK_DEV_HD static void gemm(char ta,char tb,long k,E alpha,E *a,E *b,E *c,long ld) { char x[2]={ta,0},y[2]={tb,0}; P##gemm_dev<R,C>(x,y,1,1,k,alpha,a,ld,b,ld,E(1),c,ld); } \
    MPLAPACK_DEV_HD static void rank(bool u,long n,long k,E *a,E *c,long ld,long j) { RANK<R,C>(u?"U":"L",u?"N":TRANS,n,k,R(1),a,ld,R(1),c,ld,j,j+1); } \
    MPLAPACK_DEV_HD static void gemv(long k,E *a,E *x,E *y,long ld) { const E one(1); P##gemv_dev<R,C>("N",1,k,-one,a,ld,x,1,one,y,1); } \
    static char transpose() { return TRANS[0]; } \
};
INVERSE_DEVICE_ARITHMETIC(false,R,R,"T",Rsyrk_dev)
INVERSE_DEVICE_ARITHMETIC(true,C,C,"C",Cherk_columns_dev)
#undef INVERSE_DEVICE_ARITHMETIC

// Goal 13 (continued) leaves: one block per diagonal block instead of one
// thread for the whole leaf.  Every element keeps the reference sequence of
// terms (the sweeps and the rows run in the reference order; only the
// independent elements inside one step run in parallel), so the results are
// bit-identical.
template<class E> struct trti2_leaf_job {
    // Reference Rtrti2 on one block: columns in the reference order, each
    // column's trmv sweep step propagating in parallel (every element keeps
    // the reference term order: its own sweep step (propagate, then the
    // diagonal scale), then the terms from the later sweep steps in order,
    // then the final scal).
    typedef E e_type;
    E *a; long ld,n; bool upper,unit;
    MPLAPACK_DEV_HD void block_run(int tid,int T,E *scratch=nullptr) const {
        (void)scratch;
        const E one(1.0),zero(0.0);
        if(upper) {
            // reference: for j = 1..n (1-based); 0-based column jc=j-1
            for(long jc=0;jc<n;++jc) {
                if(tid==0 && !unit) a[jc+jc*ld]=one/a[jc+jc*ld];
#if defined(__CUDA_ARCH__)
                __syncthreads();
#endif
                const E ajj=unit?-one:-a[jc+jc*ld];
                E *x=&a[jc*ld]; // column jc, rows 0..jc-1
                // trmv("U","N",n=jc elements): sweep source column js=0..jc-1
                for(long js=0;js<jc;++js) {
#if defined(__CUDA_ARCH__)
                    __syncthreads();
#endif
                    bool skip=true;
                    E temp(zero);
                    if(tid==0 && x[js]!=zero) {
                        skip=false;temp=x[js];
                        if(!unit) x[js]=x[js]*a[js+js*ld];
                    }
#if defined(__CUDA_ARCH__)
                    __shared__ bool s_skip;__shared__ E s_temp;
                    if(tid==0) { s_skip=skip;s_temp=temp; }
                    __syncthreads();
                    if(!s_skip) for(long i=tid;i<js;i+=T) x[i]=x[i]+s_temp*a[i+js*ld];
                    __syncthreads();
#else
                    if(!skip) for(long i=tid;i<js;i+=T) x[i]=x[i]+temp*a[i+js*ld];
#endif
                }
#if defined(__CUDA_ARCH__)
                __syncthreads();
#endif
                for(long i=tid;i<jc;i+=T) x[i]=x[i]*ajj;
#if defined(__CUDA_ARCH__)
                __syncthreads();
#endif
            }
        } else {
            // reference: for j = n..1 (1-based); 0-based column jc=j-1
            for(long jc=n;jc-->0;) {
                if(tid==0 && !unit) a[jc+jc*ld]=one/a[jc+jc*ld];
#if defined(__CUDA_ARCH__)
                __syncthreads();
#endif
                const E ajj=unit?-one:-a[jc+jc*ld];
                if(jc<n-1) {
                    E *x=&a[(jc+1)+jc*ld]; // column jc rows jc+1..n-1
                    const E *as=&a[(jc+1)+(jc+1)*ld];
                    const long m=n-(jc+1);
                    // trmv("L","N",m): sweep source js=m-1..0 (descending)
                    for(long js=m;js-->0;) {
#if defined(__CUDA_ARCH__)
                        __syncthreads();
#endif
                        bool skip=true;
                        E temp(zero);
                        if(tid==0 && x[js]!=zero) {
                            skip=false;temp=x[js];
                            if(!unit) x[js]=x[js]*as[js+js*ld];
                        }
#if defined(__CUDA_ARCH__)
                        __shared__ bool s_skip2;__shared__ E s_temp2;
                        if(tid==0) { s_skip2=skip;s_temp2=temp; }
                        __syncthreads();
                        if(!s_skip2) for(long i=js+1+tid;i<m;i+=T) x[i]=x[i]+s_temp2*as[i+js*ld];
                        __syncthreads();
#else
                        if(!skip) for(long i=js+1+tid;i<m;i+=T) x[i]=x[i]+temp*as[i+js*ld];
#endif
                    }
#if defined(__CUDA_ARCH__)
                    __syncthreads();
#endif
                    for(long i=tid;i<m;i+=T) x[i]=x[i]*ajj;
#if defined(__CUDA_ARCH__)
                    __syncthreads();
#endif
                }
            }
        }
    }
};
template<class E> struct lauu2_leaf_job {
    // Reference Rlauu2 on one block: rows in the reference order, the
    // independent gemv outputs of a row in parallel, each output keeping the
    // reference term order (the beta phase first, then the reference dot).
    // (The complex Clauu2 keeps the serial reference routine for now.)
    typedef E e_type;
    E *a; long ld,n; bool upper;
    MPLAPACK_DEV_HD void block_run(int tid,int T,E *aii) const {
        const E one(1.0),zero(0.0);
        if(upper) {
            // reference 1-based i_ref = i+1: dot over ROW i (columns i..n-1),
            // gemv("N", m=i, n-i): y = column i, rows 0..i-1
            for(long i=0;i<n;++i) {
                if(tid==0) aii[i]=a[i+i*ld];
#if defined(__CUDA_ARCH__)
                __syncthreads();
#endif
                if(i<n-1) {
                    if(tid==0) {
                        E s=E();
                        for(long l=0;l<n-i;++l) {
                            const E v=a[i+(i+l)*ld];
                            s=s+v*v;
                        }
                        a[i+i*ld]=s;
                    }
#if defined(__CUDA_ARCH__)
                    __syncthreads();
#endif
                    // gemv("N"): y(r)=beta*y(r), then y(r)+=alpha*x(l)*A(r,i+1+l)
                    const E beta=aii[i];
                    for(long r=tid;r<i;r+=T) {
                        E temp;
                        if(beta==one) temp=a[r+i*ld];
                        else if(beta==zero) temp=zero;
                        else temp=beta*a[r+i*ld];
                        for(long l=0;l<n-i-1;++l) temp=temp+one*a[r+(i+1+l)*ld]*a[i+(i+1+l)*ld];
                        a[r+i*ld]=temp;
                    }
#if defined(__CUDA_ARCH__)
                    __syncthreads();
#endif
                } else {
                    const E beta=aii[i];
                    if(beta!=one) for(long r=tid;r<=i;r+=T) a[r+i*ld]=beta==zero?zero:beta*a[r+i*ld];
#if defined(__CUDA_ARCH__)
                    __syncthreads();
#endif
                }
            }
        } else {
            // reference: dot over COLUMN i (rows i..n-1), gemv("T", n-i-1, i):
            // y = row i, columns 0..i-1; x = column i, rows i+1..n-1
            for(long i=0;i<n;++i) {
                if(tid==0) aii[i]=a[i+i*ld];
#if defined(__CUDA_ARCH__)
                __syncthreads();
#endif
                if(i<n-1) {
                    if(tid==0) {
                        E s=E();
                        for(long l=0;l<n-i;++l) {
                            const E v=a[(i+l)+i*ld];
                            s=s+v*v;
                        }
                        a[i+i*ld]=s;
                    }
#if defined(__CUDA_ARCH__)
                    __syncthreads();
#endif
                    // gemv("T"): y(l)=beta*y(l), then y(l)+=alpha*dot
                    const E beta=aii[i];
                    for(long l=tid;l<i;l+=T) {
                        E temp;
                        if(beta==one) temp=a[i+l*ld];
                        else if(beta==zero) temp=zero;
                        else temp=beta*a[i+l*ld];
                        E d=E();
                        for(long r=0;r<n-i-1;++r) d=d+a[i+1+r+l*ld]*a[i+1+r+i*ld];
                        a[i+l*ld]=temp+one*d;
                    }
#if defined(__CUDA_ARCH__)
                    __syncthreads();
#endif
                } else {
                    const E beta=aii[i];
                    if(beta!=one) for(long l=tid;l<=i;l+=T) a[i+l*ld]=beta==zero?zero:beta*a[i+l*ld];
#if defined(__CUDA_ARCH__)
                    __syncthreads();
#endif
                }
            }
        }
    }
};
// Goal 13 (continued) tiled left-N TRMM: one block per column of B.  The
// original column is staged in shared memory first (the reference sweep reads
// the original values), then one thread per output element computes it in the
// reference per-element order (own diagonal step first, then the propagated
// terms in the reference direction, with the reference zero skips).  Used by
// the trtri off-diagonal update; the transposed and right-side forms keep
// the lane schedule.
template<class E> struct mm_leftn_col_job {
    typedef E e_type;
    E *a,*b; long m,ld; bool upper,unit; E alpha;
    MPLAPACK_DEV_HD void block_run(int tid,int T,long c,E *xs) const {
        const E zero(0.0);
        for(long k=tid;k<m;k+=T) xs[k]=b[k+c*ld];
#if defined(__CUDA_ARCH__)
        __syncthreads();
#endif
        for(long r=tid;r<m;r+=T) {
            E v=xs[r];
            if(xs[r]!=zero) {
                v=alpha*xs[r];
                if(!unit) v=v*a[r+r*ld];
            }
            if(upper) {
                for(long k=r+1;k<m;++k) if(xs[k]!=zero) v=v+(alpha*xs[k])*a[r+k*ld];
            } else {
                // the lower sweep runs k descending, so element r receives its
                // terms from k = r-1..0 in that order
                for(long k=r-1;k>=0;--k) if(xs[k]!=zero) v=v+(alpha*xs[k])*a[r+k*ld];
            }
            b[r+c*ld]=v;
        }
    }
};
template<class A> struct inverse_job {
    using E=typename A::E;
    E *a=nullptr,*b=nullptr,*c=nullptr;
    long *status=nullptr,*piv=nullptr;
    long m=0,n=0,k=0,ld=0,col=0;
    int op=0; bool upper=false,left=false,unit=false;
    char ta='N',tb='N'; E alpha=E(1);
    MPLAPACK_DEV_HD void operator()(long lane) const {
        if(op==0) {
            *status=0;
            if(!unit) for(long i=0;i<n;++i) if(a[i+i*ld]==E(0)) { *status=i+1; break; }
        } else if(op==1) A::tri(upper,unit,n,a,ld,*status);
        else if(op==2) A::lauum(upper,n,a,ld,*status);
        else if(op==3) A::mm(left,upper,ta,unit,left?m:1,left?1:n,alpha,a,b+(left?lane*ld:lane),ld);
        else if(op==4) A::sm(upper,unit,1,n,alpha,a,b+lane,ld);
        else if(op==5) {
            const long i=lane%m,j=lane/m;
            A::gemm(ta,tb,k,alpha,a+(ta=='N'?i:i*ld),b+(tb=='N'?j*ld:j),c+i+j*ld,ld);
        } else if(op==6) A::rank(upper,n,k,a,c,ld,lane);
        else if(op==7) {
            const long i=lane%ld,j=col+lane/ld;
            if(i>j) { b[i+(j-col)*ld]=a[i+j*ld]; a[i+j*ld]=E(0); }
        } else if(op==8) A::gemv(k,a+lane,b,c+lane,ld);
        else if(op==9) {
            for(long j=n-1;j-->0;) {
                const long p=piv[j]-1;
                if(p!=j) { E t=a[lane+j*ld];a[lane+j*ld]=a[lane+p*ld];a[lane+p*ld]=t; }
            }
        }
    }
};
// Launcher owns status, synchronization and error handling; it must supply
// launch(job,lanes), feedback(status), status_data(), finish(), and error.
template<class A,class Launcher> class inverse_schedule : public Launcher {
  public:
    using J=inverse_job<A>;using E=typename A::E;
    using Launcher::error; using Launcher::launch; using Launcher::feedback;
    J make(int op,E *a,long ld) { J j{};j.op=op;j.a=a;j.ld=ld;j.status=this->status_data();return j; }
    void mm(bool left,bool upper,char tr,bool unit,long m,long n,E alpha,E *a,E *b,long ld) {
        // the left-N forms run as one thread per output element (the
        // reference sweep order per element); the transposed and right-side
        // forms keep the lane schedule
        J j=make(3,a,ld);j.left=left;j.upper=upper;j.ta=tr;j.unit=unit;j.m=m;j.n=n;j.alpha=alpha;j.b=b;
        if(m && n) launch(j,left?n:m);
    }
    void sm(bool upper,bool unit,long m,long n,E alpha,E *a,E *b,long ld) {
        J j=make(4,a,ld);j.upper=upper;j.unit=unit;j.m=m;j.n=n;j.alpha=alpha;j.b=b;
        if(n) launch(j,m);
    }
    void gemm(char ta,char tb,long m,long n,long k,E alpha,E *a,E *b,E *c,long ld) {
        J j=make(5,a,ld);j.ta=ta;j.tb=tb;j.m=m;j.n=n;j.k=k;j.alpha=alpha;j.b=b;j.c=c;
        if(k) launch(j,m*n);
    }
    void rank(bool upper,long n,long k,E *a,E *c,long ld) {
        J j=make(6,a,ld);j.upper=upper;j.n=n;j.k=k;j.c=c;if(k) launch(j,n);
    }
    long tri(bool upper,bool unit,long n,long nb,E *a,long ld) {
        J j=make(0,a,ld);j.n=n;j.unit=unit;launch(j,1);
        long info=feedback(this->status_data());if(error || info) return info;
        if(nb<=1 || nb>=n) {j.op=1;j.upper=upper;launch(j,1);return 0;}
        // (the per-iteration leaf launch below uses the block leaf job)
        const E one(1);
        for(long s=0;s<(n+nb-1)/nb && !error;++s) {
            const long col=upper?s*nb:((n-1)/nb-s)*nb,jb=std::min(nb,n-col);
            E *diag=a+col+col*ld;
            const long before=upper?col:n-col-jb;
            if(before) {
                E *off=upper?a+col*ld:a+col+jb+col*ld;
                E *inv=upper?a:a+col+jb+(col+jb)*ld;
                mm(true,upper,'N',unit,before,jb,one,inv,off,ld);
                sm(upper,unit,before,jb,-one,diag,off,ld);
            }
            if(this->block_leaves) {
                trti2_leaf_job<E> leaf{diag,ld,jb,upper,unit};this->launch_leaf(leaf);
            } else {
                J leaf=make(1,diag,ld);leaf.upper=upper;leaf.unit=unit;leaf.n=jb;launch(leaf,1);
            }
        }
        return 0;
    }
    void lauum(bool upper,long n,long nb,E *a,long ld) {
        if(nb<=1 || nb>=n) {J j=make(2,a,ld);j.upper=upper;j.n=n;launch(j,1);return;}
        for(long col=0;col<n && !error;col+=nb) {
            const long jb=std::min(nb,n-col),rest=n-col-jb;
            E *diag=a+col+col*ld,*off=upper?a+col*ld:a+col;
            mm(!upper,upper,A::transpose(),false,upper?col:jb,upper?jb:col,E(1),diag,off,ld);
            if(A::cx || !this->block_leaves) {
                J leaf=make(2,diag,ld);leaf.upper=upper;leaf.n=jb;launch(leaf,1);
            } else {
                lauu2_leaf_job<E> leaf{diag,ld,jb,upper};this->launch_leaf(leaf);
            }
            if(rest) {
                gemm(upper?'N':A::transpose(),upper?A::transpose():'N',upper?col:jb,upper?jb:col,rest,E(1),
                    upper?a+(col+jb)*ld:a+col+jb+col*ld,upper?a+col+(col+jb)*ld:a+col+jb,off,ld);
                rank(upper,jb,rest,upper?a+col+(col+jb)*ld:a+col+jb+col*ld,diag,ld);
            }
        }
    }
    long getri(long n,long nbtri,long nb,long nbmin,E *a,E *work,long *piv,long ld) {
        long info=tri(true,false,n,nbtri,a,ld);if(error || info) return info;
        const E one(1);
        if(nb<nbmin || nb>=n) {
            for(long col=n;col-->0 && !error;) {
                J copy=make(7,a,ld);copy.b=work;copy.col=col;launch(copy,n);
                if(col+1<n) {J v=make(8,a+(col+1)*ld,ld);v.k=n-col-1;v.b=work+col+1;v.c=a+col*ld;launch(v,n);}
            }
        } else for(long col=((n-1)/nb)*nb;col>=0 && !error;col-=nb) {
            long jb=std::min(nb,n-col);J copy=make(7,a,ld);copy.b=work;copy.col=col;launch(copy,n*jb);
            if(col+jb<n) gemm('N','N',n,jb,n-col-jb,-one,a+(col+jb)*ld,work+col+jb,a+col*ld,ld);
            sm(false,true,n,jb,one,work+col,a+col*ld,ld);
        }
        J swaps=make(9,a,ld);swaps.n=n;swaps.piv=piv;launch(swaps,n);return 0;
    }
};
}
#endif
