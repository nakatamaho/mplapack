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
template<class R,class X,bool Complex> bool fast(operation op,long n,X *x,long incx,X *y,long incy,R &real_result,X &dot_result) {
    if(!fast_enabled() || !stride_ok(n,incx) || !x) return false;
    const bool product=op==dot || op==dotc || op==dotu;
    if(product && (!stride_ok(n,incy) || !y)) return false;
    // ASUM excludes nonpositive increments; NRM2, like this tree's reference,
    // permits negative increments. Zero increments retain reference behavior.
    if(op==asum && incx<0) return false;
    if(Complex && op==norm && n>std::numeric_limits<long>::max()/2) return false;
    const long ix=incx<0?-(n-1)*incx:0,iy=product && incy<0?-(n-1)*incy:0;
    try {
        auto s=tree<R>(Complex && op==norm?2*n:n,op,[&](long i) {
            if constexpr(Complex) {
                if(op==norm) return leaf<R>(norm,i%2?x[ix+(i/2)*incx].imag():x[ix+(i/2)*incx].real(),R(0.0),i);
                return complex_leaf<R>(op,x[ix+i*incx],product?y[iy+i*incy]:X(0.0),i);
            } else return leaf<R>(op,x[ix+i*incx],product?y[iy+i*incy]:R(0.0),i);
        });
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
