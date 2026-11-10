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
// Goal 13 fast-mode sub-lane job: the same fixed mini-tree as the CPU
// two_level_impl (Rblas1_fast_common.h), so the results agree bit for bit.
template<class R,class X,bool Complex,operation OP> struct sublane {
    long count;const X *x,*y;state<R> *out;
    MPLAPACK_BLAS1_HD void operator()(long j) const {
        constexpr long E=fast_sub_elements;
        constexpr bool product=OP==dot || OP==dotc || OP==dotu;
        const long base=j*E,end=base+E<count?base+E:count;
#ifdef MPLAPACK_BUILD_WITH_MPFR
        // The checked MPFR value arithmetic runs the goal 11 state chain: the
        // scalar mini-tree overflows the bounded per-lane device scratch.
        state<R> s;bool first=true;
        for(long i=base;i<end;++i) {
            state<R> t;
            if constexpr(OP==norm) {
                if constexpr(Complex) t=leaf<R>(norm,i%2?x[i/2].imag():x[i/2].real(),R(0.0),i);
                else t=leaf<R>(norm,x[i],R(0.0),i);
            }
            else if constexpr(OP==asum) {
                if constexpr(Complex) t=complex_leaf<R,X>(OP,x[i],X(0.0),i);
                else t=leaf<R>(OP,x[i],R(0.0),i);
            }
            else if constexpr(Complex) t=complex_leaf<R,X>(OP,x[i],y[i],i);
            else t=leaf<R>(OP,x[i],y[i],i);
            if(first) { s=t;first=false; }
            else s=merge(OP,s,t);
        }
        out[j]=s;
#else
        if constexpr(OP==dot || OP==dotc || OP==dotu) {
            if constexpr(Complex) {
                X v[E];long m=0;
                for(long i=base;i<end;++i,++m) v[m]=OP==dotc?conj(x[i])*y[i]:x[i]*y[i];
                if(m==E) for(long st=1;st<E;st*=2) for(long k=0;k<E-st;k+=2*st) v[k]=v[k]+v[k+st];
                else for(long k=1;k<m;++k) v[0]=v[0]+v[k];
                out[j].re=v[0].real();out[j].im=v[0].imag();
            } else {
                R v[E];long m=0;
                for(long i=base;i<end;++i,++m) v[m]=x[i]*y[i];
                if(m==E) for(long st=1;st<E;st*=2) for(long k=0;k<E-st;k+=2*st) v[k]=v[k]+v[k+st];
                else for(long k=1;k<m;++k) v[0]=v[0]+v[k];
                out[j].re=v[0];
            }
        } else if constexpr(OP==asum) {
            R v[E];long m=0;
            if constexpr(Complex) for(long i=base;i<end;++i,++m) v[m]=abs(x[i].real())+abs(x[i].imag());
            else for(long i=base;i<end;++i,++m) v[m]=abs(x[i]);
            if(m==E) for(long st=1;st<E;st*=2) for(long k=0;k<E-st;k+=2*st) v[k]=v[k]+v[k+st];
            else for(long k=1;k<m;++k) v[0]=v[0]+v[k];
            out[j].re=v[0];
        } else {
            R scale(0.0),sumsq(0.0);
            for(long i=base;i<end;++i) {
                R b;
                if constexpr(Complex) b=abs(i%2?x[i/2].imag():x[i/2].real());
                else b=abs(x[i]);
                if(scale==R(0.0)) { scale=b;sumsq=R(1.0); }
                else if(b!=R(0.0)) {
                    if(scale<b) { const R t=scale/b; sumsq=R(1.0)+sumsq*(t*t);scale=b; }
                    else { const R t=b/scale; sumsq=sumsq+R(1.0)*(t*t); }
                }
            }
            out[j].re=sumsq;out[j].scale=scale;
        }
#endif
    }
};
template<class R> struct compact {
    long chunks;const state<R> *a;state<R> *out;
    MPLAPACK_DEV_HD void operator()(long c) const { out[c]=a[c*32]; }
};
template<class R,class Root> struct finish {
    operation op;state<R> *c;state<R> *out;
    MPLAPACK_DEV_HD void operator()(long) const {
        if(op==norm) c[0].re=c[0].scale*Root::sqrt(c[0].re);
        *out=c[0];
    }
};
// Goal 13 fast-mode order, fixed by the element count alone (see the same
// order in Rblas1_fast_common.h): 32 sub-lanes of 8 sequential elements per
// 256-element chunk, in-chunk adjacent-pair tree with strides 16,8,4,2,1,
// then the adjacent-pair tree over the chunk results.  The CPU and the GPU
// return bit-identical results, independent of launch shape and threads.
// a holds ceil(count/8) sub-lane states followed by ceil(count/256) chunk
// states; callers allocate max(count, subs+chunks) states.
template<class R,class X,bool Complex,class Root,class Launcher>
void schedule(Launcher &ctx,operation op,long n,const X *x,const X *y,state<R> *a) {
    long count=Complex && op==norm?2*n:n;
    const long subs=(count+fast_sub_elements-1)/fast_sub_elements,chunks=(count+fast_chunk_length-1)/fast_chunk_length;
    switch(op) {
        case dot: ctx.launch(sublane<R,X,Complex,dot>{count,x,y,a},subs);break;
        case dotc: ctx.launch(sublane<R,X,Complex,dotc>{count,x,y,a},subs);break;
        case dotu: ctx.launch(sublane<R,X,Complex,dotu>{count,x,y,a},subs);break;
        case asum: ctx.launch(sublane<R,X,Complex,asum>{count,x,y,a},subs);break;
        default: ctx.launch(sublane<R,X,Complex,norm>{count,x,y,a},subs);break;
    }
    for(long stride=1;stride<32;stride*=2)
        ctx.launch(combine<R>{op,subs,stride,a},(subs-1)/stride/2+1);
    ctx.launch(compact<R>{chunks,a,a+subs},chunks);
    for(long stride=1;stride<chunks;) {
        ctx.launch(combine<R>{op,chunks,stride,a+subs},(chunks-1)/stride/2+1);
        if(stride>chunks/2) break;
        stride*=2;
    }
    ctx.launch(finish<R,Root>{op,a+subs,a},1);
}
}
#endif
