// Deterministic BLAS1 reduction algebra, shared by CPU and device trees.
// No scheduling or environment policy here. Adjacent left/right subtrees are
// always combined in that order; there are no floating-point atomics.
#ifndef MPLAPACK_BLAS1_REDUCTION_COMMON_H
#define MPLAPACK_BLAS1_REDUCTION_COMMON_H
#ifndef MPLAPACK_BLAS1_HD
#define MPLAPACK_BLAS1_HD
#endif
namespace mplapack_blas1 __attribute__((visibility("hidden"))) {
enum operation { dot, dotc, dotu, asum, norm, iamax, axpy, copy };
template<class R> struct state {
    R re, im, scale;
    long index;
    MPLAPACK_BLAS1_HD state():re(0.0),im(0.0),scale(0.0),index(0) {}
};
template<class R> MPLAPACK_BLAS1_HD state<R> leaf(operation op,const R &x,const R &y,long index) {
    state<R> s;
    if(op==dot) s.re=x*y;
    else if(op==norm) {
        s.scale=abs(x);
        if(s.scale!=R(0.0)) s.re=R(1.0);else s.scale=R(0.0);
    }
    else s.re=abs(x);
    s.index=index;
    return s;
}
template<class R,class C> MPLAPACK_BLAS1_HD state<R> complex_leaf(operation op,const C &x,const C &y,long index) {
    state<R> s;
    if(op==dotc || op==dotu) {
        C p;
        if(op==dotc) p=conj(x)*y;else p=x*y;
        s.re=p.real();s.im=p.imag();
    } else s.re=abs(x.real())+abs(x.imag());
    s.index=index;
    return s;
}
template<class R> MPLAPACK_BLAS1_HD state<R> merge(operation op,state<R> a,const state<R> &b) {
    if(op==iamax) {
        // Left subtree precedes right subtree; strict comparison keeps first tie.
        // Entries reject NaNs before launch (CPU iamax keeps its existing path).
        if(b.re>a.re) return b;
    } else if(op==norm) {
        if(a.scale==R(0.0)) return b;
        if(b.scale==R(0.0)) return a;
        if(a.scale<b.scale) {
            R t=a.scale/b.scale;
            a.re=b.re+a.re*(t*t);a.scale=b.scale;
        } else {
            R t=b.scale/a.scale;
            a.re=a.re+b.re*(t*t);
        }
    } else { a.re+=b.re;a.im+=b.im; }
    return a;
}
}
#endif
