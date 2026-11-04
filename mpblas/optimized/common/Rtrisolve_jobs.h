// Device-pointer triangular solve schedule. Column order is the BLAS reference
// order; only independent updates/products run in parallel. No host arithmetic.
#ifndef MPLAPACK_TRISOLVE_JOBS_H
#define MPLAPACK_TRISOLVE_JOBS_H
#include "Rblas_dev_base.h"
#include "blas_dev/tri_arithmetic.h"
namespace mplapack_trisolve __attribute__((visibility("hidden"))) {
struct shape {
    long n,k,lda; int storage; // 0 full, 1 band, 2 packed
    bool upper,notrans,conjugate,unit;
    MPLAPACK_DEV_HD long first(long j) const { return upper?(storage==1 && j>k?j-k:0):j+1; }
    MPLAPACK_DEV_HD long last(long j) const { return upper?j:(storage==1 && k<n-1-j?j+k+1:n); }
    MPLAPACK_DEV_HD long index(long i,long j) const {
        if(storage==0) return i+j*lda;
        if(storage==1) return (upper?k+i-j:i-j)+j*lda;
        return upper?j*(j+1)/2+i:j*n-j*(j-1)/2+i-j;
    }
    long count() const { return storage==2?n*(n+1)/2:n*lda; }
};
template<class V,class A> struct pivot_job {
    shape s;long j;const V *a;V *x,*p;
    MPLAPACK_DEV_HD void operator()(long) const {
        if(s.notrans) {
            p[s.n]=x[j]; // retain pre-division zero-skip even after underflow
            if(x[j]!=V(0.0) && !s.unit) x[j]=A::divide(x[j],a[s.index(j,j)]);
        } else {
            V temp=x[j];
            if(s.upper) for(long i=s.first(j);i<s.last(j);++i) temp=temp-p[i];
            else for(long i=s.last(j);i-->s.first(j);) temp=temp-p[i];
            if(!s.unit) {
                V d=a[s.index(j,j)];
                temp=A::divide(temp,s.conjugate?A::conjugate(d):d);
            }
            x[j]=temp;
        }
    }
};
template<class V,class A> struct column_job {
    shape s;long j;const V *a;V *x,*p;
    MPLAPACK_DEV_HD void operator()(long lane) const {
        long i=s.first(j)+lane;
        if(s.notrans) {
            if(p[s.n]!=V(0.0)) x[i]=A::update(x[i],a[s.index(i,j)],x[j]);
        } else {
            V d=a[s.index(i,j)];
            p[i]=A::product(x[i],s.conjugate?A::conjugate(d):d);
        }
    }
};
template<class V,class A,class L> void schedule(L &ctx,const shape &s,const V *a,V *x,V *p) {
    for(long step=0;step<s.n && !ctx.error;++step) {
        long j=s.upper==s.notrans?s.n-1-step:step;
        long count=s.last(j)-s.first(j);
        if(s.notrans) ctx.launch(pivot_job<V,A>{s,j,a,x,p},1);
        if(count) ctx.launch(column_job<V,A>{s,j,a,x,p},count);
        if(!s.notrans) ctx.launch(pivot_job<V,A>{s,j,a,x,p},1);
    }
}
}
#endif
