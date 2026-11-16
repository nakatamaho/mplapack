// Reference-order TRSM axpy forms, cooperatively scheduled within a block.
// No reductions or change to the scalar arithmetic/zero-skip rules.
#ifndef MPLAPACK_TRSM_TEAM_COMMON_H
#define MPLAPACK_TRSM_TEAM_COMMON_H
#include "Rblas_lanes_common.h"
#if defined(__CUDACC__) && !defined(MPLAPACK_MPFR_CUDA_HOST_EMULATION) && !defined(MPLAPACK_GEMM_CUDA_HOST_EMULATION)
#include <cuda_runtime.h>
#endif
namespace mplapack_blas {
// Upper transpose panels can stage each solved row for independent output
// accumulators. Every accumulator still subtracts k in increasing order.
MPLAPACK_BLAS_HD inline bool trsm_dot_shape(const blas_shape &s) {
    return s.op==L3_TRSM && s.left && s.upper && s.trans && s.m<=64;
}
MPLAPACK_BLAS_HD inline bool trsm_team_shape(const blas_shape &s) {
    return s.op==L3_TRSM && ((s.left && !s.trans) || (!s.left && s.trans) || trsm_dot_shape(s));
}
template<class E> struct trsm_team_args {
    E e; blas_shape s; blas_mats<E> x;
    const typename E::V *reciprocal=nullptr;
};
template<class E> struct trsm_team_shared {
    typename E::V pivot[16];
    int active[16];
};
template<class E, bool Conjugate> MPLAPACK_BLAS_HD typename E::V trsm_team_a(
    E &e, const typename E::M &a,long i,long j) {
    auto v=e.reg(5); e.set(v,e.at(a,i,j));
    if constexpr(Conjugate) e.conjugate(v,v);
    return v;
}
// The lane is a column (Left) or a row (Right). coord splits the other axis.
template<class E> MPLAPACK_BLAS_HD void trsm_team_init(trsm_team_args<E> &g,long lane,int coord,int stride) {
    auto &e=g.e; auto &s=g.s; auto &x=g.x;
    if(!s.left || lane>=s.n || s.alpha_one) return;
    auto y=e.reg(0);
    for(long i=coord;i<s.m;i+=stride) {
        e.mul(y,e.at(x.al,0,0),e.at(x.c,i,lane)); e.put(x.c,i,lane,y);
    }
}
template<class E,bool Conjugate> MPLAPACK_BLAS_HD void trsm_team_pivot(
    trsm_team_args<E> &g,trsm_team_shared<E> &sh,long lane,int slot,long k) {
    auto &e=g.e; auto &s=g.s; auto &x=g.x;
    sh.active[slot]=0;
    if(lane>=(s.left?s.n:s.m)) return;
    auto y=e.reg(0),t=e.reg(1),o=e.reg(2);
    if(s.left) {
        e.set(y,e.at(x.c,k,lane));
        // Keep the test BEFORE division: underflow must not suppress updates.
        if(e.iszero(y)) return;
        sh.active[slot]=1;
        if(s.nounit) { e.div(y,y,e.at(x.a,k,k)); e.put(x.c,k,lane,y); }
    } else {
        sh.active[slot]=1;
        e.set(y,e.at(x.c,lane,k));
        if(s.nounit) {
            if(g.reciprocal) e.set(t,g.reciprocal[k]);
            else { e.one(o); e.div(t,o,trsm_team_a<E,Conjugate>(e,x.a,k,k)); }
            e.mul(y,t,y); e.put(x.c,lane,k,y);
        }
    }
    e.set(sh.pivot[slot],y);
}
template<class E,bool Conjugate> MPLAPACK_BLAS_HD void trsm_team_update(
    trsm_team_args<E> &g,const trsm_team_shared<E> &sh,long lane,int slot,long k,int coord,int stride) {
    auto &e=g.e; auto &s=g.s; auto &x=g.x;
    if(!sh.active[slot]) return;
    const long end=s.left?s.m:s.n;
    const long first=s.upper?0:k+1, last=s.upper?k:end;
    auto y=e.reg(0),p=e.reg(1),v=e.reg(2);
    for(long u=first+coord;u<last;u+=stride) {
        if(s.left) {
            e.mul(p,sh.pivot[slot],e.at(x.a,u,k));
            e.sub(y,e.at(x.c,u,lane),p); e.put(x.c,u,lane,y);
        } else {
            // The reference tests the un-conjugated A, then multiplies A*B.
            if(e.iszero(e.at(x.a,u,k))) continue;
            e.set(v,trsm_team_a<E,Conjugate>(e,x.a,u,k));
            e.mul(p,v,sh.pivot[slot]);
            e.sub(y,e.at(x.c,lane,u),p); e.put(x.c,lane,u,y);
        }
    }
}
template<class E> MPLAPACK_BLAS_HD void trsm_team_finish(trsm_team_args<E> &g,long lane,long k) {
    auto &e=g.e; auto &s=g.s; auto &x=g.x;
    if(s.left || lane>=s.m || s.alpha_one) return;
    auto y=e.reg(0);
    e.mul(y,e.at(x.al,0,0),e.at(x.c,lane,k)); e.put(x.c,lane,k,y);
}
#if defined(__CUDACC__) && !defined(MPLAPACK_MPFR_CUDA_HOST_EMULATION) && !defined(MPLAPACK_GEMM_CUDA_HOST_EMULATION)
template<class E> MPLAPACK_BLAS_HD constexpr int trsm_team_width(bool left) {
    return sizeof(typename E::V)==16?(left?1:4):sizeof(typename E::V)<=64?4:(left?4:16);
}
static __device__ inline void trsm_team_barrier() {
    if(blockDim.x*blockDim.y==32) __syncwarp(); else __syncthreads();
}
template<class E,bool Conjugate=false,bool Checked=false,bool Local=false> __device__ void trsm_team_body(trsm_team_args<E> g,int *failed) {
    __shared__ __align__(16) unsigned char storage[sizeof(trsm_team_shared<E>)];
    auto &sh=*reinterpret_cast<trsm_team_shared<E>*>(storage);
    // Split independent columns (Left) / rows (Right) into small teams.
    // BF uses more teams for p64 occupancy; large MPFR values keep fewer teams.
    const int slot=g.s.left?threadIdx.y:threadIdx.x;
    const int coord=g.s.left?threadIdx.x:threadIdx.y;
    const int width=g.s.left?blockDim.y:blockDim.x;
    const int stride=g.s.left?blockDim.x:blockDim.y;
    long lane=(long)blockIdx.x*width+slot;
    const long global_lane=lane;
    const long lanes=g.s.left?g.s.n:g.s.m;
    const long order=g.s.left?g.s.m:g.s.n;
    auto original=g.x.c;
    if constexpr(Local) {
        constexpr int W=trsm_team_width<E>(false);
        __shared__ __align__(16) unsigned char matrix_bytes[64*W*sizeof(typename E::V)];
        auto *matrix=reinterpret_cast<typename E::V*>(matrix_bytes);
        if(global_lane<lanes) for(long u=coord;u<order;u+=stride)
            g.e.set(matrix[g.s.left?u+slot*order:slot+u*width],g.s.left?g.e.at(original,u,global_lane):g.e.at(original,global_lane,u));
        g.x.c.p=matrix; g.x.c.ld=g.s.left?order:width;
        const long remaining=lanes-(long)blockIdx.x*width;
        if(g.s.left) g.s.n=remaining<width?remaining:width;
        else g.s.m=remaining<width?remaining:width;
        lane=slot;
    }
    trsm_team_init(g,lane,coord,stride);
    trsm_team_barrier();
    for(long step=0;step<order;++step) {
        const long k=g.s.upper?order-1-step:step;
        if(coord==0) trsm_team_pivot<E,Conjugate>(g,sh,lane,slot,k);
        trsm_team_barrier();
        trsm_team_update<E,Conjugate>(g,sh,lane,slot,k,coord,stride);
        trsm_team_barrier();
        if(coord==0) trsm_team_finish(g,lane,k);
        // Pivot values are staged; finish writes only the completed column.
        // The next pivot does not access that column, so no third barrier.
    }
    if constexpr(Local) {
        trsm_team_barrier();
        if(global_lane<lanes) for(long u=coord;u<order;u+=stride) {
            if(g.s.left) g.e.put(original,u,global_lane,g.e.at(g.x.c,u,slot));
            else g.e.put(original,global_lane,u,g.e.at(g.x.c,slot,u));
        }
    }
    if constexpr(Checked) if(!g.e.ok) atomicExch(failed,1);
}
template<class E,bool Conjugate,bool Checked> __device__ void trsm_dot_body(trsm_team_args<E> g,int *failed) {
    using V=typename E::V;
    __shared__ __align__(16) unsigned char storage[2*sizeof(V)];
    auto *pivot=reinterpret_cast<V*>(storage);
    const long i=threadIdx.x,j=2L*blockIdx.x+threadIdx.y;
    const bool mine=i<g.s.m && j<g.s.n;
    auto value=g.e.reg(0),term=g.e.reg(1);
    // Unlike NoTrans, the reference multiplies alpha even when it equals one.
    if(mine) g.e.mul(value,g.e.at(g.x.al,0,0),g.e.at(g.x.c,i,j));
    for(long k=0;k<g.s.m;++k) {
        if(mine && i==k) {
            if(g.s.nounit) g.e.div(value,value,trsm_team_a<E,Conjugate>(g.e,g.x.a,k,k));
            g.e.set(pivot[threadIdx.y],value);
            g.e.put(g.x.c,i,j,value);
        }
        __syncthreads();
        if(mine && i>k) {
            g.e.mul(term,trsm_team_a<E,Conjugate>(g.e,g.x.a,k,i),pivot[threadIdx.y]);
            g.e.sub(value,value,term); // no zero skip in the dot-product form
        }
        __syncthreads();
    }
    if constexpr(Checked) if(!g.e.ok) atomicExch(failed,1);
}
template<class E,bool Conjugate,bool Checked> __global__ void k_trsm_dot(trsm_team_args<E> g,int *failed) {
    trsm_dot_body<E,Conjugate,Checked>(g,failed);
}
template<class E,bool Conjugate=false,bool Checked=false,bool Local=false> __global__ void k_trsm_team(trsm_team_args<E> g,int *failed) {
    trsm_team_body<E,Conjugate,Checked,Local>(g,failed);
}
template<class E,bool Conjugate=false,bool Checked=false> static int trsm_team_plain_device(
    trsm_team_args<E> g,cudaStream_t stream=nullptr,int *failed=nullptr) {
    const long lanes=g.s.left?g.s.n:g.s.m;
    if(!lanes || !g.s.m || !g.s.n) return 0;
    if(trsm_dot_shape(g.s)) {
        k_trsm_dot<E,Conjugate,Checked><<<(unsigned)((g.s.n+1)/2),dim3(64,2),0,stream>>>(g,failed);
        return cudaGetLastError();
    }
    const int width=trsm_team_width<E>(g.s.left);
    const dim3 threads=g.s.left?dim3(32,width):dim3(width,128/width);
    if constexpr(sizeof(typename E::V)<=64) {
        if((g.s.left?g.s.m:g.s.n)<=64) {
            k_trsm_team<E,Conjugate,Checked,true><<<(unsigned)((lanes+width-1)/width),threads,0,stream>>>(g,failed);
            return cudaGetLastError();
        }
    }
    k_trsm_team<E,Conjugate,Checked><<<(unsigned)((lanes+width-1)/width),threads,0,stream>>>(g,failed);
    return cudaGetLastError();
}
#else
template<class E,bool Conjugate=false,bool Checked=false> static int trsm_team_plain_device(trsm_team_args<E> g) {
    const long lanes=g.s.left?g.s.n:g.s.m,order=g.s.left?g.s.m:g.s.n;
    if(trsm_dot_shape(g.s)) {
        for(long j=0;j<g.s.n;++j) {
            typename E::V value[64];
            for(long i=0;i<g.s.m;++i) g.e.mul(value[i],g.e.at(g.x.al,0,0),g.e.at(g.x.c,i,j));
            for(long k=0;k<g.s.m;++k) {
                if(g.s.nounit) g.e.div(value[k],value[k],trsm_team_a<E,Conjugate>(g.e,g.x.a,k,k));
                g.e.put(g.x.c,k,j,value[k]);
                for(long i=k+1;i<g.s.m;++i) {
                    auto term=g.e.reg(1);
                    g.e.mul(term,trsm_team_a<E,Conjugate>(g.e,g.x.a,k,i),value[k]);
                    g.e.sub(value[i],value[i],term);
                }
            }
        }
        if constexpr(Checked) return g.e.ok?0:1;
        return 0;
    }
    const int width=g.s.left?4:16,stride=g.s.left?32:8;
    for(long base=0;base<lanes;base+=width) {
        trsm_team_shared<E> sh;
        for(int slot=0;slot<width;++slot) for(int c=0;c<stride;++c)
            trsm_team_init(g,base+slot,c,stride);
        for(long step=0;step<order;++step) {
            const long k=g.s.upper?order-1-step:step;
            for(int slot=0;slot<width;++slot) trsm_team_pivot<E,Conjugate>(g,sh,base+slot,slot,k);
            for(int slot=0;slot<width;++slot) for(int c=0;c<stride;++c)
                trsm_team_update<E,Conjugate>(g,sh,base+slot,slot,k,c,stride);
            for(int slot=0;slot<width;++slot) trsm_team_finish(g,base+slot,k);
        }
    }
    if constexpr(Checked) return g.e.ok?0:1;
    return 0;
}
#endif
}
#include "Rtrsm_block_common.h"
#endif
