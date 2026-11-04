// Device-pointer BLAS1 jobs. Host emulation executes these exact jobs.
#ifndef MPLAPACK_BLAS1_JOBS_COMMON_H
#define MPLAPACK_BLAS1_JOBS_COMMON_H
#define MPLAPACK_BLAS1_HD MPLAPACK_DEV_HD
#include "Rblas1_reduction_common.h"
namespace mplapack_blas1 __attribute__((visibility("hidden"))) {
template<class R,class X,bool Complex> struct leaves {
    operation op;long n;const X *x,*y;state<R> *out;
    MPLAPACK_DEV_HD void operator()(long i) const {
        const bool product=op==dot || op==dotc || op==dotu;
        if constexpr(Complex) {
            if(op==norm) out[i]=leaf<R>(norm,i%2?x[i/2].imag():x[i/2].real(),R(0.0),i);
            else out[i]=complex_leaf<R>(op,x[i],product?y[i]:X(0.0),i);
        } else out[i]=leaf<R>(op,x[i],product?y[i]:R(0.0),i);
    }
};
template<class R> struct combine {
    operation op;long n,stride;state<R> *a;
    MPLAPACK_DEV_HD void operator()(long q) const {
        long i=q*stride*2;
        if(i<n-stride) a[i]=merge(op,a[i],a[i+stride]);
    }
};
template<class R,class Root> struct finish {
    operation op;state<R> *a;
    MPLAPACK_DEV_HD void operator()(long) const {
        if(op==norm) a[0].re=a[0].scale*Root::sqrt(a[0].re);
    }
};
template<class R,class X,bool Complex,class Root,class Launcher>
void schedule(Launcher &ctx,operation op,long n,const X *x,const X *y,state<R> *a) {
    long count=Complex && op==norm?2*n:n;
    ctx.launch(leaves<R,X,Complex>{op,n,x,y,a},count);
    for(long stride=1;stride<count;) {
        ctx.launch(combine<R>{op,count,stride,a},(count-1)/stride/2+1);
        if(stride>count/2) break;
        stride*=2;
    }
    ctx.launch(finish<R,Root>{op,a},1);
}
}
#endif
