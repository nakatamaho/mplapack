// Device-pointer triangular solve schedule. Column order is the BLAS reference
// order; only independent updates/products run in parallel. No host arithmetic.
#ifndef MPLAPACK_TRISOLVE_JOBS_H
#define MPLAPACK_TRISOLVE_JOBS_H
#include "Rblas_dev_base.h"
#include "blas_dev/tri_arithmetic.h"
#if defined(__CUDACC__) && !defined(MPLAPACK_GEMM_CUDA_HOST_EMULATION)
#include <cuda/atomic>
#endif
#include <cstdio>
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
// Goal 13 (continued) wavefront solve: flag-driven true wavefront in one
// persistent block.  The columns no longer synchronize the whole block; each
// element has one owner thread that applies its terms in the exact reference
// order as soon as the x(j) they use is final, detected by a per-element flag.
// The owner of element i spins on its dependencies (prefetching the A element
// of the next term while it waits), applies each term when the flag appears,
// and at its own turn performs the zero test / division and publishes its
// flag.  The pivot sequence across elements is forced by the dependencies to
// the reference column order, and every element sees exactly the reference
// sequence of terms, so the result is bit-identical.
//
// Flag values: 0 pending; 1 final (terms may use x(j)); 2 final and, for the
// NoTrans forms, the pre-division x(j) was zero — the reference applies no
// terms of such a column at all (the zero test is on the pre-division value,
// retained even if the division itself would underflow to zero).
//
// Dependencies of element i: for the NoTrans forms the columns j whose update
// range covers i (below i for Lower, above for Upper, processed in the
// reference direction); for the Trans/C forms the product range of the own
// column (above i for Lower, below for Upper).  Both cases share: the range
// sits below i and runs ascending iff notrans != upper, above i and
// descending otherwise.  Band storage only narrows the range.
template<class V,class A> struct wave_job {
    typedef V vtype;
    shape s;const V *a;V *x,*p;
    MPLAPACK_DEV_HD bool below() const { return s.notrans!=s.upper; }
    MPLAPACK_DEV_HD long dep_lo(long i) const {
        if(!below()) return i+1;
        return s.storage==1 && i>s.k ? i-s.k : 0;
    }
    MPLAPACK_DEV_HD long dep_hi(long i) const {
        if(below()) return i;
        return s.storage==1 && s.k<s.n-1-i ? i+s.k+1 : s.n;
    }
    MPLAPACK_DEV_HD int flag_get(const int *flags,long i) const {
#if defined(__CUDA_ARCH__)
        // relaxed poll; the acquire load of flag_acquire orders the x(j)
        // reads that follow
        return cuda::atomic_ref<const int,cuda::thread_scope_device>(flags[i]).load(cuda::memory_order_relaxed);
#else
        return flags[i];
#endif
    }
    MPLAPACK_DEV_HD int flag_acquire(const int *flags,long i) const {
#if defined(__CUDA_ARCH__)
        return cuda::atomic_ref<const int,cuda::thread_scope_device>(flags[i]).load(cuda::memory_order_acquire);
#else
        return flags[i];
#endif
    }
    // Publish with release: every write the flag advertises (x(j), division
    // result) is ordered before it.  The flags live in global memory (device
    // scope) because the wavefront spreads the elements over all resident
    // blocks.
    MPLAPACK_DEV_HD void flag_set(int *flags,long i,int v) const {
#if defined(__CUDA_ARCH__)
        cuda::atomic_ref<int,cuda::thread_scope_device>(flags[i]).store(v,cuda::memory_order_release);
#else
        flags[i]=v;
#endif
    }
    // The pivot/parallel pair with the p[] workspace stays for the flat
    // per-column schedule of the MPFR runtime precision (bounded per-lane
    // scratch; the flag wavefront would keep one live value per owner).
    MPLAPACK_DEV_HD void pivot(long j) {
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
    MPLAPACK_DEV_HD void parallel(long j,long lane,long T) {
        const long i0=s.first(j),i1=s.last(j);
        if(s.notrans) {
            for(long i=i0+lane;i<i1;i+=T)
                if(p[s.n]!=V(0.0)) x[i]=A::update(x[i],a[s.index(i,j)],x[j]);
        } else {
            for(long i=i0+lane;i<i1;i+=T) {
                V d=a[s.index(i,j)];
                p[i]=A::product(x[i],s.conjugate?A::conjugate(d):d);
            }
        }
    }
    // Event-driven owner of a contiguous chunk of elements.  The owner keeps
    // every owned element's term stream as current as it can: each pass
    // applies, head first (the element that finalizes first), every term whose
    // x(j) is already final, and pivots every completed element.  Only when a
    // whole pass makes no progress does it spin, on the head's next dependency
    // (prefetching that coefficient while it waits).  This is what keeps the
    // critical path at one term per column even when a thread owns several
    // elements: the terms of a later element are applied as they arrive, not
    // after the earlier element finished, which would put n/2 serial terms per
    // column back on the chain.  The running dot-form accumulation lives in
    // x(e) itself — nobody reads x(e) before its pivot — so no per-element
    // register survives across passes; prog(e) records the applied-term
    // count (owner-private, kept next to the flags).  Every element still
    // receives its terms in the exact reference order, so the result is
    // bit-identical.  With all blocks resident (cooperative launch) the
    // dependency spin is acyclic and cannot deadlock; the caller zeroes the
    // flags and prog before the launch.
    MPLAPACK_DEV_HD void block_run(int tid,int T,int *flags,int *prog) {
        // Strided ownership: element e = tid + k*T.  Consecutive elements land
        // on consecutive threads, so the working threads spread over all
        // resident blocks; a thread owns ceil(n/T) elements staggered T columns
        // apart, which the pass loop keeps current.
        if(tid>=s.n) return;
        const bool asc=below();
        long pending=0;
        {
            long e=tid;
            while(e<s.n) { ++pending;e+=T; }
        }
        const long cnt=pending;
        if(cnt==1) {
            // One element per thread (the common grid case): a tight loop with
            // all state in registers — no pass scan, no prog traffic.  Each
            // iteration prefetches the next coefficient, spins until its x(j)
            // is final (the wait covers the load), applies the term and, at
            // the element's own turn, divides and publishes.
            const long e=tid;
            const long lo=dep_lo(e),hi=dep_hi(e),d=hi-lo;
            const V dg=a[s.index(e,e)];
            if(s.notrans) {
                for(long t=0;t<d;++t) {
                    const long j=asc?lo+t:hi-1-t;
                    const V dd=a[s.index(e,j)];
                    int st;
#if defined(__CUDA_ARCH__)
                    while((st=flag_acquire(flags,j))==0) __nanosleep(32);
#endif
                    if(st==1) x[e]=A::update(x[e],dd,x[j]);
                }
                const V pre=x[e];
                if(pre!=V(0.0) && !s.unit) x[e]=A::divide(pre,dg);
                flag_set(flags,e,pre==V(0.0)?2:1);
            } else {
                for(long t=0;t<d;++t) {
                    const long j=asc?lo+t:hi-1-t;
                    const V dd=a[s.index(j,e)];
                    int st;
#if defined(__CUDA_ARCH__)
                    while((st=flag_acquire(flags,j))==0) __nanosleep(32);
#endif
                    x[e]=x[e]-A::product(x[j],s.conjugate?A::conjugate(dd):dd);
                }
                if(!s.unit) x[e]=A::divide(x[e],s.conjugate?A::conjugate(dg):dg);
                flag_set(flags,e,1);
            }
            return;
        }
        while(pending) {
            bool progress=false;
            for(long k=0;k<cnt;++k) {
                const long e=tid+(asc?k:cnt-1-k)*T;
                if(e>=s.n) continue;
                const long lo=dep_lo(e),hi=dep_hi(e),d=hi-lo;
                long t=prog[e];
                if(s.notrans) {
                    while(t<d) {
                        const long j=asc?lo+t:hi-1-t;
                        const int st=flag_get(flags,j);
                        if(!st) break;
                        if(st==1) {
                            if(flag_acquire(flags,j)==1) x[e]=A::update(x[e],a[s.index(e,j)],x[j]);
                        }
                        ++t;progress=true;
                    }
                    if(t==d && flag_get(flags,e)==0) {
                        const V pre=x[e];
                        if(pre!=V(0.0) && !s.unit) x[e]=A::divide(pre,a[s.index(e,e)]);
                        flag_set(flags,e,pre==V(0.0)?2:1);
                        --pending;progress=true;
                    }
                } else {
                    while(t<d) {
                        const long j=asc?lo+t:hi-1-t;
                        if(!flag_get(flags,j)) break;
                        if(!flag_acquire(flags,j)) break;
                        const V dd=a[s.index(j,e)];
                        x[e]=x[e]-A::product(x[j],s.conjugate?A::conjugate(dd):dd);
                        ++t;progress=true;
                    }
                    if(t==d && flag_get(flags,e)==0) {
                        if(!s.unit) {
                            const V dd=a[s.index(e,e)];
                            x[e]=A::divide(x[e],s.conjugate?A::conjugate(dd):dd);
                        }
                        flag_set(flags,e,1);
                        --pending;progress=true;
                    }
                }
                prog[e]=t;
            }
            if(!progress) {
                // No owned element can advance: the head's next dependency is
                // not ours (an owned dependency of the head would precede the
                // head in finalize order and has been pivoted above), so spin
                // until another thread publishes it.  Prefetch the coefficient
                // of the pending term AND the diagonal of the head, so that on
                // wake the whole column finish (term, division, publish) runs
                // without further global-load latency: this is the critical
                // path of the whole solve.
                const long e=tid+(asc?0:cnt-1)*T;
                const long lo=dep_lo(e),hi=dep_hi(e);
                const long t0=prog[e];
                const long j=asc?lo+t0:hi-1-t0;
                const long d=hi-lo;
                int st=1;
                V dg=V(0.0);
#if defined(__CUDA_ARCH__)
                if(s.notrans) {
                    const V dd=a[s.index(e,j)];
                    dg=a[s.index(e,e)];
                    while((st=flag_acquire(flags,j))==0) __nanosleep(32);
                    if(st==1) x[e]=A::update(x[e],dd,x[j]);
                } else {
                    const V dd=a[s.index(j,e)];
                    dg=a[s.index(e,e)];
                    while((st=flag_acquire(flags,j))==0) __nanosleep(32);
                    x[e]=x[e]-A::product(x[j],s.conjugate?A::conjugate(dd):dd);
                }
#else
                (void)j;
#endif
                prog[e]=t0+1;
                if(t0+1==d) {
                    // the head received its last term on this wake: finish it
                    // with the prefetched diagonal
                    if(s.notrans) {
                        const V pre=x[e];
                        if(pre!=V(0.0) && !s.unit) x[e]=A::divide(pre,dg);
                        flag_set(flags,e,pre==V(0.0)?2:1);
                    } else {
                        if(!s.unit) x[e]=A::divide(x[e],s.conjugate?A::conjugate(dg):dg);
                        flag_set(flags,e,1);
                    }
                    --pending;
                }
            }
        }
    }
};
template<class V,class A,class L> void schedule(L &ctx,const shape &s,const V *a,V *x,V *p) {
    ctx.launch_wave(wave_job<V,A>{s,a,x,p});
}
}
#endif
