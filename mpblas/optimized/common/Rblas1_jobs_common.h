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
// Goal 13 (continued) fast-mode chunk job: one chunk of
// fast_chunk_length elements per block iteration, folded in shared memory in
// the exact two-level order of the CPU (Rblas1_fast_common.h): the elements
// are staged coalesced into shared memory, each sub-lane computes its fixed
// mini-tree (sequential state chain for the MPFR checked arithmetic), the
// fast_sub_lanes sub-lane states are combined by the adjacent-pair tree with
// ascending strides, and the chunk result (the slot-0 state) is written once.
// The slot-0 fold matches the global adjacent-pair tree over all sub-lanes
// because the chunk boundaries are multiples of fast_sub_lanes: every merge
// of the global tree with both operands inside one chunk reduces that
// chunk's slot 0 in the same order, and the tail chunk merges only its valid
// sub-lanes.  The staged layout halves the global traffic (the old schedule
// wrote ceil(n/E) sub-lane states and re-read them per tree stage); the
// per-stage kernel launches disappear.
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

// Goal 13 (continued) chunk job (below): used by the launchers that opt in
// (chunked_reductions); the checked MPFR arithmetic keeps the lane-parallel
// sub-lane schedule above (its state chain needs the pool-wide concurrency,
// and its bounded per-lane scratch forbids the staged shared-memory design).
template<class R,class X,bool Complex,operation OP> struct chunk_block {
    typedef X vtype;
    long count,chunks;const X *x,*y;state<R> *out;
    static constexpr bool product=OP==dot||OP==dotc||OP==dotu;
    // The element pre-product type: dot/dotu/dotc keep the (complex) product,
    // asum collapses to |re|+|im| up front.  norm and the MPFR checked
    // arithmetic run the serial sub-lane chains instead (their per-element
    // step is not a single value).
    static constexpr bool preproduct=product || OP==asum;
#ifdef MPLAPACK_BUILD_WITH_MPFR
    static constexpr bool serial=true; // checked arithmetic: state chain
#else
    static constexpr bool serial=!preproduct; // norm
#endif
    typedef typename std::conditional<preproduct && !Complex,R,X>::type ptype;
    // The pre-product path stages pass_elems elements per pass (whole
    // sub-lanes; several passes per chunk when the chunk's products exceed
    // the shared-memory budget — the sub-lane trees are per sub-lane, so the
    // passes do not change any rounding sequence; the in-chunk state tree
    // runs once after all passes).
    static constexpr long pass_elems() {
        // at least one sub-lane, at most a ~16 KB product array
        return (long)16384/sizeof(ptype)/fast_sub_elements*fast_sub_elements<fast_sub_elements
            ?fast_sub_elements
            :(long)16384/sizeof(ptype)/fast_sub_elements*fast_sub_elements;
    }
    static constexpr size_t smem_bytes() {
        return serial
            ?(product?2:1)*fast_chunk_length*sizeof(X)+fast_sub_lanes*sizeof(state<R>)
            :pass_elems()*sizeof(ptype)+fast_sub_lanes*sizeof(state<R>);
    }
    static constexpr size_t sublane_smem_bytes() { return fast_sub_lanes*sizeof(state<R>); }
    // Stream slot i -> the x element holding it.  Only the complex norm is
    // component-wise (count=2n slots over n complex elements); every other
    // form is element-wise.
    MPLAPACK_BLAS1_HD X xat(long i) const {
        if constexpr(Complex && OP==norm) return x[i/2];
        else return x[i];
    }
    MPLAPACK_BLAS1_HD ptype pre(long i) const {
        // the per-element pre-operation of the fixed order
        if constexpr(OP==asum) {
            if constexpr(Complex) return abs(x[i].real())+abs(x[i].imag());
            else return abs(x[i]);
        } else {
            if constexpr(Complex) return OP==dotc?conj(x[i])*y[i]:x[i]*y[i];
            else return x[i]*y[i];
        }
    }
    // Serial sub-lane reduction (norm and the MPFR checked arithmetic), from
    // the staged elements when sh != 0, from global memory otherwise.
    MPLAPACK_BLAS1_HD void sublane_serial(long t,long base,long len,const X *sx,const X *sy,state<R> *shs) const {
        constexpr long E=fast_sub_elements;
        const long sb=base+t*E,end=sb+E<count?sb+E:count;
#ifdef MPLAPACK_BUILD_WITH_MPFR
        // The checked MPFR value arithmetic runs the goal 11 state chain: the
        // scalar mini-tree overflows the bounded per-lane device scratch.
        state<R> s;bool first=true;
        for(long i=sb;i<end;++i) {
            state<R> u;
            if constexpr(OP==norm) {
                if constexpr(Complex) u=leaf<R>(norm,i%2?sx[i-base].imag():sx[i-base].real(),R(0.0),i);
                else u=leaf<R>(norm,sx[i-base],R(0.0),i);
            }
            else if constexpr(OP==asum) {
                if constexpr(Complex) u=complex_leaf<R,X>(OP,sx[i-base],X(0.0),i);
                else u=leaf<R>(OP,sx[i-base],R(0.0),i);
            }
            else if constexpr(Complex) u=complex_leaf<R,X>(OP,sx[i-base],sy[i-base],i);
            else u=leaf<R>(OP,sx[i-base],sy[i-base],i);
            if(first) { s=u;first=false; }
            else s=merge(OP,s,u);
        }
        shs[t]=s;
#else
        R scale(0.0),sumsq(0.0);
        for(long i=sb;i<end;++i) {
            R b;
            if constexpr(Complex) b=abs(i%2?sx[i-base].imag():sx[i-base].real());
            else b=abs(sx[i-base]);
            if(scale==R(0.0)) { scale=b;sumsq=R(1.0); }
            else if(b!=R(0.0)) {
                if(scale<b) { const R u=scale/b; sumsq=R(1.0)+sumsq*(u*u);scale=b; }
                else { const R u=b/scale; sumsq=sumsq+R(1.0)*(u*u); }
            }
        }
        shs[t].re=sumsq;shs[t].scale=scale;
#endif
    }
    // sh holds smem_bytes() bytes (staged) or sublane_smem_bytes() bytes
    // (unstaged: the elements are read straight from global memory; the same
    // sub-lane partition and trees keep the order identical).
    MPLAPACK_BLAS1_HD void block_run(int tid,int T,long c,char *sh) const {
        const long base=c*fast_chunk_length;
        const long len=base+fast_chunk_length<count?fast_chunk_length:count-base;
        state<R> *shs;
        if constexpr(serial) {
            const X *sx,*sy;
            if(sh) {
                X *sxh=reinterpret_cast<X *>(sh);
                X *syh=product?sxh+fast_chunk_length:sxh; // unused when not a product
                shs=reinterpret_cast<state<R> *>((product?syh:sxh)+fast_chunk_length);
                for(long i=tid;i<len;i+=T) {
                    sxh[i]=xat(base+i);
                    if(product) syh[i]=y[base+i];
                }
#if defined(__CUDA_ARCH__)
                __syncthreads();
#endif
                sx=sxh;sy=syh;
            } else {
                sx=x+base;sy=product?y+base:y;
                shs=reinterpret_cast<state<R> *>(sh);
            }
            const long valid=(len+fast_sub_elements-1)/fast_sub_elements;
            for(long t=tid;t<fast_sub_lanes && t<valid;t+=T) sublane_serial(t,base,len,sx,sy,shs);
        } else {
            // pre-product path: one value per element in shared memory, then
            // the fixed mini-trees as parallel shared-memory stages — the
            // pairing of every add is the CPU's, so the result is bit-identical
            ptype *pr=reinterpret_cast<ptype *>(sh);
            shs=reinterpret_cast<state<R> *>(pr+pass_elems());
            constexpr long PE=pass_elems();
            const long passes=len<PE?1:(len+PE-1)/PE;
            for(long pass=0;pass<passes;++pass) {
                const long pb=pass*PE;
                const long plen=pb+PE<len?PE:len-pb; // elements of this pass
                for(long i=tid;i<plen;i+=T) pr[i]=pre(base+pb+i);
#if defined(__CUDA_ARCH__)
                __syncthreads();
#endif
                const long full=plen/fast_sub_elements; // full sub-lanes in the pass
                for(long stride=1;stride<fast_sub_elements;stride*=2) {
#if defined(__CUDA_ARCH__)
                    __syncthreads();
#endif
                    const long per=fast_sub_elements/((long)2*stride);
                    for(long idx=tid;idx<full*per;idx+=T) {
                        const long s=idx/per,k=(idx%per)*(2*stride);
                        pr[(size_t)s*fast_sub_elements+k]=pr[(size_t)s*fast_sub_elements+k]+pr[(size_t)s*fast_sub_elements+k+stride];
                    }
                }
#if defined(__CUDA_ARCH__)
                __syncthreads();
#endif
                // sub-lane states: full sub-lanes from the tree root; the
                // single partial sub-lane (the last of the chunk) keeps the
                // CPU's left-to-right chain
                const long scount=PE/fast_sub_elements; // sub-lanes staged per pass
                for(long t=tid;t<scount;t+=T) {
                    const long gs=pass*scount+t; // global sub-lane
                    if(gs*fast_sub_elements>=len) break;
                    if(gs*fast_sub_elements+fast_sub_elements<=len) {
                        if constexpr(Complex) { shs[gs].re=pr[(size_t)t*fast_sub_elements].real();shs[gs].im=pr[(size_t)t*fast_sub_elements].imag(); }
                        else shs[gs].re=pr[(size_t)t*fast_sub_elements];
                    } else {
                        ptype v=pr[(size_t)t*fast_sub_elements];
                        for(long k=(long)t*fast_sub_elements+1;k<plen;++k) v=v+pr[k];
                        if constexpr(Complex) { shs[gs].re=v.real();shs[gs].im=v.imag(); }
                        else shs[gs].re=v;
                    }
                }
#if defined(__CUDA_ARCH__)
                __syncthreads();
#endif
            }
        }
#if defined(__CUDA_ARCH__)
        __syncthreads();
#endif
        // the in-chunk adjacent-pair tree over the sub-lane states; the chunk
        // boundaries are multiples of fast_sub_lanes, so these are exactly the
        // global tree's intra-chunk merges
        const long valid=(len+fast_sub_elements-1)/fast_sub_elements;
        for(long stride=1;stride<fast_sub_lanes;stride*=2) {
            for(long t=tid;t+stride<valid;t+=T)
                if(t%((long)2*stride)==0) shs[t]=merge(OP,shs[t],shs[t+stride]);
#if defined(__CUDA_ARCH__)
            __syncthreads();
#endif
        }
        if(tid==0) out[c]=shs[0];
    }
};
template<class R,class Root> struct finish {
    operation op;state<R> *c;state<R> *out;
    MPLAPACK_DEV_HD void operator()(long) const {
        if(op==norm) c[0].re=c[0].scale*Root::sqrt(c[0].re);
        *out=c[0];
    }
};
// Goal 13 (continued) chunk-result tree: one block folds the chunk states in
// shared memory with the same adjacent-pair tree as the CPU (strides
// 1,2,4,... over the chunk states, in place), then applies the norm finish.
// The block_run returns through the launcher only when the chunk states fit
// the shared memory budget; otherwise the schedule keeps the per-stage
// combine launches over global memory (the order is identical).
template<class R,class Root> struct chunk_tree {
    typedef state<R> vtype;
    operation op;long chunks;state<R> *a,*out;
    static constexpr size_t smem_bytes(long units) { return sizeof(state<R>)*size_t(units); }
    MPLAPACK_BLAS1_HD void block_run(int tid,int T,long,char *sh) const {
        state<R> *shs=reinterpret_cast<state<R> *>(sh);
        for(long i=tid;i<chunks;i+=T) shs[i]=a[i];
#if defined(__CUDA_ARCH__)
        __syncthreads();
#endif
        for(long stride=1;stride<chunks;) {
            for(long q=tid;q<(chunks-1)/stride/2+1;q+=T) {
                const long i=q*stride*2;
                if(i<chunks-stride) shs[i]=merge(op,shs[i],shs[i+stride]);
            }
#if defined(__CUDA_ARCH__)
            __syncthreads();
#endif
            if(stride>chunks/2) break;
            stride*=2;
        }
        if(tid==0) {
            if(op==norm) shs[0].re=shs[0].scale*Root::sqrt(shs[0].re);
            *out=shs[0];
        }
    }
};
template<class R,class X,bool Complex,class Root,class Launcher>
void schedule(Launcher &ctx,operation op,long n,const X *x,const X *y,state<R> *a) {
    long count=Complex && op==norm?2*n:n;
    const long subs=(count+fast_sub_elements-1)/fast_sub_elements,chunks=(count+fast_chunk_length-1)/fast_chunk_length;
    if(!ctx.chunked_reductions) {
        // lane-parallel sub-lane schedule (checked MPFR arithmetic): the same
        // order, staged through the sub-lane states
        switch(op) {
            case dot: ctx.launch(sublane<R,X,Complex,dot>{count,x,y,a},subs);break;
            case dotc: ctx.launch(sublane<R,X,Complex,dotc>{count,x,y,a},subs);break;
            case dotu: ctx.launch(sublane<R,X,Complex,dotu>{count,x,y,a},subs);break;
            case asum: ctx.launch(sublane<R,X,Complex,asum>{count,x,y,a},subs);break;
            default: ctx.launch(sublane<R,X,Complex,norm>{count,x,y,a},subs);break;
        }
        for(long stride=1;stride<fast_sub_lanes;stride*=2)
            ctx.launch(combine<R>{op,subs,stride,a},(subs-1)/stride/2+1);
        ctx.launch(compact<R>{chunks,a,a+subs},chunks);
        for(long stride=1;stride<chunks;) {
            ctx.launch(combine<R>{op,chunks,stride,a+subs},(chunks-1)/stride/2+1);
            if(stride>chunks/2) break;
            stride*=2;
        }
        ctx.launch(finish<R,Root>{op,a+subs,a},1);
        return;
    }
    switch(op) {
        case dot: ctx.launch_chunked(chunk_block<R,X,Complex,dot>{count,chunks,x,y,a},chunks);break;
        case dotc: ctx.launch_chunked(chunk_block<R,X,Complex,dotc>{count,chunks,x,y,a},chunks);break;
        case dotu: ctx.launch_chunked(chunk_block<R,X,Complex,dotu>{count,chunks,x,y,a},chunks);break;
        case asum: ctx.launch_chunked(chunk_block<R,X,Complex,asum>{count,chunks,x,y,a},chunks);break;
        default: ctx.launch_chunked(chunk_block<R,X,Complex,norm>{count,chunks,x,y,a},chunks);break;
    }
    // The chunk states sit directly at a[0..chunks): the chunk kernel has
    // already folded every chunk's sub-lanes, so the old sub-lane tree and
    // compact step are gone.  One shared-memory tree kernel folds the chunk
    // states when they fit; otherwise the per-stage combines run as before
    // (the same order, just with launches per stage).
    if(ctx.launch_chunk_tree(chunk_tree<R,Root>{op,chunks,a,a},chunks)==0) return;
    for(long stride=1;stride<chunks;) {
        ctx.launch(combine<R>{op,chunks,stride,a},(chunks-1)/stride/2+1);
        if(stride>chunks/2) break;
        stride*=2;
    }
    ctx.launch(finish<R,Root>{op,a,a},1);
}
}
#endif
