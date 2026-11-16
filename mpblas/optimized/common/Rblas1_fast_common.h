// Opt-in CPU BLAS1 reductions. Default paths never enter this code.
#ifndef MPLAPACK_BLAS1_FAST_COMMON_H
#define MPLAPACK_BLAS1_FAST_COMMON_H
#include <vector>
#include <limits>
#include <stdexcept>
#include "Rblas_cpu_common.h"
#include "Rblas1_reduction_common.h"
#include "Rblas1_policy.h"
namespace mplapack_blas1 __attribute__((visibility("hidden"))) {
template<class F> void parallel(long count,F f) {
    const mplapack_blas::blas_thread_context ctx;
#ifdef _OPENMP
#pragma omp parallel if(count>=256 && mplapack_blas::blas_env_long("MPLAPACK_BLAS_OMP",1)!=0)
#endif
    {
        const mplapack_blas::blas_thread_context::scope scope(ctx);
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
        for(long i=0;i<count;++i) f(i);
    }
}
template<class R,class F> state<R> tree(long n,operation op,F make) {
    std::vector<state<R>> a(static_cast<size_t>(n));
    parallel(n,[&](long i) { a[i]=make(i); });
    // In-place sparse tree: only the leftmost element of a subtree is written.
    // Thus a stage never overwrites another worker's input.
    for(long stride=1;stride<n;) {
        const long groups=(n-1)/stride/2+1;
        parallel(groups,[&](long q) {
            const long i=q*stride*2;
            if(i<n-stride) a[i]=merge(op,a[i],a[i+stride]);
        });
        if(stride>n/2) break;
        stride*=2;
    }
    return a[0];
}
// Goal 13 fast-mode order, fixed by the element count alone so that the CPU
// and the GPU return bit-identical results for every thread and launch shape:
// contiguous chunks of 256 elements are split into 32 sub-lanes of 8 elements
// (sub-lane s takes chunk*256 + s + 32*k, k = 0..7, a coalesced pattern on the
// GPU), each sub-lane is summed sequentially in index order, the 32 sub-lane
// states of a chunk are combined by the adjacent-pair tree with strides
// 16,8,4,2,1, and the chunk results are combined by the same adjacent-pair
// tree.  Every element-to-result path has a fixed number of rounding steps,
// independent of the thread count.
// Goal 13: word-level BF arithmetic for the fast in-lane accumulation on
// dd/td/qd (the same libQD3 branch-free word operations the SIMD BLAS uses;
// verified bit-identical to the type operators by the qd_ops test, and to the
// GPU arithmetic by the bit-identity premise).  One double per word.
#if defined(MPLAPACK_BUILD_WITH_DD)
#define MPLAPACK_FAST_WORDS 2
namespace mplapack_blas1 {
inline void fast_mul_w(const REAL &a,const REAL &b,double *r) { mplapack_gemm::dd_bf_mul_w(a.x[0],a.x[1],b.x[0],b.x[1],r[0],r[1]); }
inline void fast_add_w(const double *a,const double *b,double *r) { mplapack_gemm::dd_bf_add_w(a[0],a[1],b[0],b[1],r[0],r[1]); }
// the reference abs keeps (high >= 0) values unchanged, including signed zero
inline void fast_abs_w(const REAL &a,double *r) { if(a.x[0]<0.0) { r[0]=-a.x[0];r[1]=-a.x[1]; } else { r[0]=a.x[0];r[1]=a.x[1]; } }
inline REAL fast_word_real(const double *w) { REAL r;r.x[0]=w[0];r.x[1]=w[1];return r; }
}
#elif defined(MPLAPACK_BUILD_WITH_TD)
#define MPLAPACK_FAST_WORDS 3
namespace mplapack_blas1 {
inline void fast_mul_w(const REAL &a,const REAL &b,double *r) { mplapack_gemm::td_bf_mul_w(a.x[0],a.x[1],a.x[2],b.x[0],b.x[1],b.x[2],r[0],r[1],r[2]); }
inline void fast_add_w(const double *a,const double *b,double *r) { mplapack_gemm::td_bf_add_w(a[0],a[1],a[2],b[0],b[1],b[2],r[0],r[1],r[2]); }
inline void fast_abs_w(const REAL &a,double *r) { if(a.x[0]<0.0) for(long q=0;q<3;++q) r[q]=-a.x[q]; else for(long q=0;q<3;++q) r[q]=a.x[q]; }
inline REAL fast_word_real(const double *w) { REAL r;for(long q=0;q<3;++q) r.x[q]=w[q];return r; }
}
#elif defined(MPLAPACK_BUILD_WITH_QD)
#define MPLAPACK_FAST_WORDS 4
namespace mplapack_blas1 {
inline void fast_mul_w(const REAL &a,const REAL &b,double *r) { mplapack_gemm::qd_bf_mul_w(a.x[0],a.x[1],a.x[2],a.x[3],b.x[0],b.x[1],b.x[2],b.x[3],r[0],r[1],r[2],r[3]); }
inline void fast_add_w(const double *a,const double *b,double *r) { mplapack_gemm::qd_bf_add_w(a[0],a[1],a[2],a[3],b[0],b[1],b[2],b[3],r[0],r[1],r[2],r[3]); }
inline void fast_abs_w(const REAL &a,double *r) { if(a.x[0]<0.0) for(long q=0;q<4;++q) r[q]=-a.x[q]; else for(long q=0;q<4;++q) r[q]=a.x[q]; }
inline REAL fast_word_real(const double *w) { REAL r;for(long q=0;q<4;++q) r.x[q]=w[q];return r; }
}
#endif
template<class R,operation OP,class FX,class FY> state<R> two_level_impl(long count,FX xat,FY yat) {
    using X=decltype(xat(0L));
    constexpr bool C=!std::is_same<X,R>::value;
    constexpr long E=fast_sub_elements;
    const long subs=(count+E-1)/E;
    std::vector<state<R>> a(static_cast<size_t>(subs));
    const mplapack_blas::blas_thread_context ctx;
    const bool use_omp=count>=256 && mplapack_blas::blas_env_long("MPLAPACK_BLAS_OMP",1)!=0;
#ifdef _OPENMP
#pragma omp parallel if(use_omp)
#endif
    {
        const mplapack_blas::blas_thread_context::scope scope(ctx);
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
        for(long j=0;j<subs;++j) {
            const long base=j*E,end=base+E<count?base+E:count;
#ifdef MPLAPACK_BUILD_WITH_MPFR
            // The checked MPFR arithmetic uses the goal 11 state chain (the
            // device jobs of the same backend run exactly this order).
            state<R> s2;bool first2=true;
            for(long i=base;i<end;++i) {
                state<R> t;
                if constexpr(OP==norm) {
                    // the switch below instantiates every op; norm with a
                    // complex element type never runs (the caller passes the
                    // component accessor) but must compile
                    if constexpr(C) t=leaf<R>(norm,xat(i).real(),R(0.0),i);
                    else t=leaf<R>(norm,xat(i),R(0.0),i);
                }
                else if constexpr(OP==asum) {
                    if constexpr(C) t=complex_leaf<R,X>(OP,xat(i),X(0.0),i);
                    else t=leaf<R>(OP,xat(i),R(0.0),i);
                }
                else if constexpr(C) t=complex_leaf<R,X>(OP,xat(i),yat(i),i);
                else t=leaf<R>(OP,xat(i),yat(i),i);
                if(first2) { s2=t;first2=false; }
                else s2=merge(OP,s2,t);
            }
            a[(size_t)j]=s2;
#else
            // The sub-lane result: a fixed pairwise mini-tree over the E slots
            // when the sub-lane is full, the same left-to-right chain over the
            // existing values otherwise.  Identical code runs on the GPU, so
            // CPU and GPU results agree bit for bit for every n.
            if constexpr(OP==dot || OP==dotc || OP==dotu) {
                if constexpr(C) {
                    X v[E];long m=0;
                    for(long i=base;i<end;++i,++m) v[m]=OP==dotc?conj(xat(i))*yat(i):xat(i)*yat(i);
                    if(m==E) for(long st=1;st<E;st*=2) for(long k=0;k<E-st;k+=2*st) v[k]=v[k]+v[k+st];
                    else for(long k=1;k<m;++k) v[0]=v[0]+v[k];
                    a[(size_t)j].re=v[0].real();a[(size_t)j].im=v[0].imag();
                } else {
#ifdef MPLAPACK_FAST_WORDS
                    constexpr long W=MPLAPACK_FAST_WORDS;
                    double v[E*W];long m=0;
                    for(long i=base;i<end;++i,++m) mplapack_blas1::fast_mul_w(xat(i),yat(i),v+m*W);
                    if(m==E) for(long st=1;st<E;st*=2) for(long k=0;k<E-st;k+=2*st) mplapack_blas1::fast_add_w(v+k*W,v+(k+st)*W,v+k*W);
                    else for(long k=1;k<m;++k) mplapack_blas1::fast_add_w(v,v+k*W,v);
                    a[(size_t)j].re=mplapack_blas1::fast_word_real(v);
#else
                    R v[E];long m=0;
                    for(long i=base;i<end;++i,++m) v[m]=xat(i)*yat(i);
                    if(m==E) for(long st=1;st<E;st*=2) for(long k=0;k<E-st;k+=2*st) v[k]=v[k]+v[k+st];
                    else for(long k=1;k<m;++k) v[0]=v[0]+v[k];
                    a[(size_t)j].re=v[0];
#endif
                }
            } else if constexpr(OP==asum) {
#ifdef MPLAPACK_FAST_WORDS
                constexpr long W=MPLAPACK_FAST_WORDS;
                double v[E*W];long m=0;
                if constexpr(C) for(long i=base;i<end;++i,++m) { mplapack_blas1::fast_abs_w(xat(i).real(),v+m*W); double t[W]; mplapack_blas1::fast_abs_w(xat(i).imag(),t); mplapack_blas1::fast_add_w(v+m*W,t,v+m*W); }
                else for(long i=base;i<end;++i,++m) mplapack_blas1::fast_abs_w(xat(i),v+m*W);
                if(m==E) for(long st=1;st<E;st*=2) for(long k=0;k<E-st;k+=2*st) mplapack_blas1::fast_add_w(v+k*W,v+(k+st)*W,v+k*W);
                else for(long k=1;k<m;++k) mplapack_blas1::fast_add_w(v,v+k*W,v);
                a[(size_t)j].re=mplapack_blas1::fast_word_real(v);
#else
                R v[E];long m=0;
                if constexpr(C) for(long i=base;i<end;++i,++m) v[m]=abs(xat(i).real())+abs(xat(i).imag());
                else for(long i=base;i<end;++i,++m) v[m]=abs(xat(i));
                if(m==E) for(long st=1;st<E;st*=2) for(long k=0;k<E-st;k+=2*st) v[k]=v[k]+v[k+st];
                else for(long k=1;k<m;++k) v[0]=v[0]+v[k];
                a[(size_t)j].re=v[0];
#endif
            } else {
                R scale(0.0),sumsq(0.0);
                for(long i=base;i<end;++i) {
                    const R b=abs(xat(i));
                    if(scale==R(0.0)) { scale=b;sumsq=R(1.0); }
                    else if(b!=R(0.0)) {
                        if(scale<b) { const R t=scale/b; sumsq=R(1.0)+sumsq*(t*t);scale=b; }
                        else { const R t=b/scale; sumsq=sumsq+R(1.0)*(t*t); }
                    }
                }
                a[(size_t)j].re=sumsq;a[(size_t)j].scale=scale;
            }
#endif
        }
        for(long stride=1;stride<fast_sub_lanes;stride*=2) {
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
            for(long q=0;q<(subs-1)/stride/2+1;++q) {
                const long i=q*stride*2;
                if(i<subs-stride) a[(size_t)i]=merge(OP,a[(size_t)i],a[(size_t)(i+stride)]);
            }
        }
    }
    const long chunks=(count+fast_chunk_length-1)/fast_chunk_length;
    if(chunks==1) return a[0];
    std::vector<state<R>> b(static_cast<size_t>(chunks));
    parallel(chunks,[&](long c) { b[(size_t)c]=a[(size_t)(c*fast_sub_lanes)]; });
    for(long stride=1;stride<chunks;) {
        const long groups=(chunks-1)/stride/2+1;
        parallel(groups,[&](long q) {
            const long i=q*stride*2;
            if(i<chunks-stride) b[(size_t)i]=merge(OP,b[(size_t)i],b[(size_t)(i+stride)]);
        });
        if(stride>chunks/2) break;
        stride*=2;
    }
    return b[0];
}
template<class R,class FX,class FY> state<R> two_level(long count,operation op,FX xat,FY yat) {
    switch(op) {
        case dot: return two_level_impl<R,dot>(count,xat,yat);
        case dotc: return two_level_impl<R,dotc>(count,xat,yat);
        case dotu: return two_level_impl<R,dotu>(count,xat,yat);
        case asum: return two_level_impl<R,asum>(count,xat,yat);
        default: return two_level_impl<R,norm>(count,xat,yat);
    }
}

template<class R,class X,bool Complex> bool fast(operation op,long n,X *x,long incx,X *y,long incy,R &real_result,X &dot_result) {
    if(!fast_enabled() || !stride_ok(n,incx) || !x) return false;
    const bool product=op==dot || op==dotc || op==dotu;
    if(product && (!stride_ok(n,incy) || !y)) return false;
    // ASUM excludes nonpositive increments; NRM2, like this tree's reference,
    // permits negative increments. Zero increments retain reference behavior.
    if(op==asum && incx<0) return false;
    if(Complex && op==norm && n>std::numeric_limits<long>::max()/2) return false;
    const long ix=incx<0?-(n-1)*incx:0,iy=product && incy<0?-(n-1)*incy:0;
    const long count=Complex && op==norm?2*n:n;
    try {
        state<R> s;
        if(op==iamax) {
            s=tree<R>(Complex && op==norm?2*n:n,op,[&](long i) {
                if constexpr(Complex) {
                    if(op==norm) return leaf<R>(norm,i%2?x[ix+(i/2)*incx].imag():x[ix+(i/2)*incx].real(),R(0.0),i);
                    return complex_leaf<R>(op,x[ix+i*incx],product?y[iy+i*incy]:X(0.0),i);
                } else return leaf<R>(op,x[ix+i*incx],product?y[iy+i*incy]:R(0.0),i);
            });
        } else if(op==norm) {
            auto xat=[&](long i)->R {
                if constexpr(Complex) return i%2?x[ix+(i/2)*incx].imag():x[ix+(i/2)*incx].real();
                else return x[ix+i*incx];
            };
            s=two_level<R>(count,op,xat,xat);
        } else {
            auto xat=[&](long i)->X { return x[ix+i*incx]; };
            auto yat=[&](long i)->X { return product?y[iy+i*incy]:X(0.0); };
            s=two_level<R>(count,op,xat,yat);
        }
        if(product) {
            if constexpr(Complex) dot_result=X(s.re,s.im);
            else dot_result=s.re;
        } else real_result=op==norm?s.scale*sqrt(s.re):s.re;
        return true;
    } catch(const std::bad_alloc &) { return false; }
      catch(const std::length_error &) { return false; }
}
}
#endif
