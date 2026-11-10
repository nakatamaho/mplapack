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
// Goal 13 wavefront solve: one persistent single-block kernel runs all the
// column steps (no launch per column).  Every element keeps the reference
// sequence of terms: NoTrans finalizes x(j) (pre-division zero skip, division)
// and then all pending updates x(i) -= x(j)*A(i,j) run in the thread-strided
// phase; Trans/C forms compute the products in parallel and fold them in the
// reference order on one thread, then divide.  The updates of the independent
// elements of a column run concurrently while the columns themselves proceed
// one after another, exactly as the reference schedule orders them, so the
// result is bit-identical.  Emulation runs block_run(0,1): the same phases
// serially (the sync points are no-ops off the device).
#if defined(__CUDA_ARCH__)
#define MPLAPACK_WAVE_SYNC __syncthreads()
#else
#define MPLAPACK_WAVE_SYNC
#endif
template<class V,class A> struct wave_job {
    typedef V vtype;
    shape s;const V *a;V *x,*p;
    MPLAPACK_DEV_HD const V &at(long i,long j,const V *c,long cap) const { return cap>=s.n?c[i]:a[s.index(i,j)]; }
    MPLAPACK_DEV_HD void prefetch(long j,long tid,long T,V *c) {
        for(long i=s.first(j)+tid;i<s.last(j);i+=T) c[i]=a[s.index(i,j)];
    }
    MPLAPACK_DEV_HD void pivot(long j,const V *dg,long cap) {
        if(s.notrans) {
            p[s.n]=x[j]; // retain pre-division zero-skip even after underflow
            if(x[j]!=V(0.0) && !s.unit) x[j]=A::divide(x[j],at(j,j,dg,cap));
        } else {
            V temp=x[j];
            if(s.upper) for(long i=s.first(j);i<s.last(j);++i) temp=temp-p[i];
            else for(long i=s.last(j);i-->s.first(j);) temp=temp-p[i];
            if(!s.unit) {
                V d=at(j,j,dg,cap);
                temp=A::divide(temp,s.conjugate?A::conjugate(d):d);
            }
            x[j]=temp;
        }
    }
    MPLAPACK_DEV_HD void parallel(long j,long tid,long T,const V *c,long cap) {
        const long i0=s.first(j),i1=s.last(j);
        if(s.notrans) {
            for(long i=i0+tid;i<i1;i+=T)
                if(p[s.n]!=V(0.0)) x[i]=A::update(x[i],at(i,j,c,cap),x[j]);
        } else {
            for(long i=i0+tid;i<i1;i+=T) {
                V d=at(i,j,c,cap);
                p[i]=A::product(x[i],s.conjugate?A::conjugate(d):d);
            }
        }
    }
    // c0/c1 hold whole columns of A (double buffered, indexed by absolute i),
    // dg the diagonal; cap < n disables the staging (global loads instead).
    MPLAPACK_DEV_HD void block_run(int tid,int T,V *c0,V *c1,V *dg,long cap) {
        const bool pf=cap>=s.n && c0 && c1 && dg;
        if(pf) for(long j=tid;j<s.n;j+=T) dg[j]=a[s.index(j,j)];
        const bool uppercol=s.upper==s.notrans;
        long cur=0;
        if(pf) prefetch(uppercol?s.n-1:0,tid,T,c0);
        MPLAPACK_WAVE_SYNC;
        for(long step=0;step<s.n;++step) {
            const long j=uppercol?s.n-1-step:step;
            const long jn=step+1<s.n?(uppercol?j-1:j+1):-1;
            V *ccur=pf?(cur==0?c0:c1):nullptr,*cnext=pf?(cur==0?c1:c0):nullptr;
            if(s.notrans) {
                if(tid==0) pivot(j,dg,cap);
                MPLAPACK_WAVE_SYNC;
                if(s.last(j)>s.first(j)) parallel(j,tid,T,ccur,cap);
                if(pf && jn>=0) prefetch(jn,tid,T,cnext);
                MPLAPACK_WAVE_SYNC;
            } else {
                if(s.last(j)>s.first(j)) parallel(j,tid,T,ccur,cap);
                MPLAPACK_WAVE_SYNC;
                if(tid==0) pivot(j,dg,cap);
                if(pf && jn>=0) prefetch(jn,tid,T,cnext);
                MPLAPACK_WAVE_SYNC;
            }
            cur^=1;
        }
    }
};
template<class V,class A,class L> void schedule(L &ctx,const shape &s,const V *a,V *x,V *p) {
    ctx.launch_wave(wave_job<V,A>{s,a,x,p});
}
}
#endif
