// Deterministic BLAS1 reduction algebra, shared by CPU and device trees.
// No scheduling or environment policy here. Adjacent left/right subtrees are
// always combined in that order; there are no floating-point atomics.
#ifndef MPLAPACK_BLAS1_REDUCTION_COMMON_H
#define MPLAPACK_BLAS1_REDUCTION_COMMON_H
#ifndef MPLAPACK_BLAS1_HD
#define MPLAPACK_BLAS1_HD
#endif
namespace mplapack_blas1 __attribute__((visibility("hidden"))) {
// Goal 13 fast-mode two-level order, fixed by the element count alone: a chunk
// of fast_chunk_length = fast_sub_lanes * fast_sub_elements elements is split
// into fast_sub_lanes sub-lanes of fast_sub_elements sequential elements
// (sub-lane s takes chunk*fast_chunk_length + s + fast_sub_lanes*k, a
// coalesced pattern on the GPU), each sub-lane is reduced by a fixed pairwise
// mini-tree, the sub-lane states of a chunk are combined by the adjacent-pair
// tree with ascending strides, and the chunk results are combined by the same
// adjacent-pair tree.  The CPU and the GPU run this exact order.
// The sub-lane length is 32 for the BF backends (their values are register
// sized) and 8 for MPFR (the checked device arithmetic has a bounded per-lane
// scratch; 32 live complex values overflow the launch budget).  The CPU and
// the GPU of one backend always share the same value.
#ifdef MPLAPACK_BUILD_WITH_MPFR
constexpr long fast_sub_elements=8;
#else
constexpr long fast_sub_elements=32;
#endif
constexpr long fast_sub_lanes=32;
constexpr long fast_chunk_length=fast_sub_elements*fast_sub_lanes;
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
