// Copyright (c) 2026 Nakata, Maho. BSD-2-Clause (see Rgetrf.cpp).
// Backend-independent reference-order inverse jobs and schedule.
#ifndef MPLAPACK_INVERSE_SCHEDULE_H
#define MPLAPACK_INVERSE_SCHEDULE_H
#include <algorithm>
#include "inverse_dev/inverse_dev.h"
namespace mplapack_inverse_gpu {
using namespace mplapack_blas_dev;
template<class R,class C,bool Complex> struct inverse_arithmetic;
#define INVERSE_DEVICE_ARITHMETIC(FLAG,TYPE,P,TRANS,RANK) \
template<class R,class C> struct inverse_arithmetic<R,C,FLAG> { \
    using E=TYPE; \
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
            J leaf=make(1,diag,ld);leaf.upper=upper;leaf.unit=unit;leaf.n=jb;launch(leaf,1);
        }
        return 0;
    }
    void lauum(bool upper,long n,long nb,E *a,long ld) {
        if(nb<=1 || nb>=n) {J j=make(2,a,ld);j.upper=upper;j.n=n;launch(j,1);return;}
        for(long col=0;col<n && !error;col+=nb) {
            const long jb=std::min(nb,n-col),rest=n-col-jb;
            E *diag=a+col+col*ld,*off=upper?a+col*ld:a+col;
            mm(!upper,upper,A::transpose(),false,upper?col:jb,upper?jb:col,E(1),diag,off,ld);
            J leaf=make(2,diag,ld);leaf.upper=upper;leaf.n=jb;launch(leaf,1);
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
