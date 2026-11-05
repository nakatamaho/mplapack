// Blocked reference-order axpy TRSM. The trailing update SUBTRACTS every
// product separately; it is not GEMM with a pre-summed dot product.
#ifndef MPLAPACK_TRSM_BLOCK_COMMON_H
#define MPLAPACK_TRSM_BLOCK_COMMON_H
#include "Rtrsm_team_common.h"
#include <vector>
namespace mplapack_blas {
template<class E> struct trsm_tile_size {
    static constexpr int value=sizeof(typename E::V)<=64?16:sizeof(typename E::V)<=256?8:4;
};
static MPLAPACK_BLAS_HD inline bool trsm_block_eligible(const blas_shape &s) {
    // Nonunit Left needs the pre-division zero mask if division underflows.
    // Keep its already parallel team path rather than infer that mask later.
    return trsm_team_shape(s) && (!s.left || !s.nounit) &&
        (s.left?s.m>64:s.n>=64) && (s.left?s.n:s.m)>=32;
}
template<class E> struct trsm_update_shared {
    static constexpr int T=trsm_tile_size<E>::value;
    typename E::V a[T][T+1];
    typename E::V b[T][T+1];
};
template<class E> MPLAPACK_BLAS_HD void trsm_scale_element(trsm_team_args<E> &g,long i,long j) {
    auto p=g.e.reg(0);
    g.e.mul(p,g.e.at(g.x.al,0,0),g.e.at(g.x.c,i,j)); g.e.put(g.x.c,i,j,p);
}
template<class E,bool Conjugate> MPLAPACK_BLAS_HD void trsm_reciprocal_element(trsm_team_args<E> &g,typename E::V *r,long k) {
    auto o=g.e.reg(0),v=g.e.reg(1); g.e.one(o);
    g.e.div(v,o,trsm_team_a<E,Conjugate>(g.e,g.x.a,k,k)); g.e.set(r[k],v);
}
template<class E> MPLAPACK_BLAS_HD void trsm_update_load(trsm_team_args<E> &g,trsm_update_shared<E> &sh,
    long i0,long j0,long start,long end,long p0,int ti,int tj) {
    auto &e=g.e; auto &x=g.x;
    const long i=i0+ti,j=j0+ti;
    const long k=g.s.upper?end-1-p0-tj:start+p0+tj;
    // a[l][row]: A(i,k) for Left, B(i,k) for Right. Coalesced rows.
    if(p0+tj<end-start && i<g.s.m) e.set(sh.a[tj][ti],g.s.left?e.at(x.a,i,k):e.at(x.c,i,k));
    if(g.s.left) {
        const long kb=g.s.upper?end-1-p0-ti:start+p0+ti;
        const long col=j0+tj;
        if(p0+ti<end-start && col<g.s.n) e.set(sh.b[ti][tj],e.at(x.c,kb,col));
    } else if(p0+tj<end-start && j<g.s.n) e.set(sh.b[tj][ti],e.at(x.a,j,k));
}
template<class E,bool Conjugate> MPLAPACK_BLAS_HD void trsm_update_step(trsm_team_args<E> &g,
    const trsm_update_shared<E> &sh,typename E::V &acc,long count,int ti,int tj) {
    auto &e=g.e; auto p=e.reg(1),v=e.reg(2);
    for(long l=0;l<count;++l) {
        // Left: the saved, unit-diagonal pivot. Right: un-conjugated A(j,k).
        if(e.iszero(sh.b[l][tj])) continue;
        e.set(v,sh.b[l][tj]);
        if constexpr(Conjugate) if(!g.s.left) e.conjugate(v,v);
        e.mul(p,v,sh.a[l][ti]); e.sub(acc,acc,p);
    }
}
template<class E> static trsm_team_args<E> trsm_diagonal_args(trsm_team_args<E> g,long start,long count) {
    g.x.a.p+=start+start*g.x.a.ld;
    g.x.c.p+=g.s.left?start:start*g.x.c.ld;
    if(g.reciprocal) g.reciprocal+=start;
    if(g.s.left) g.s.m=count; else g.s.n=count;
    g.s.alpha_one=1; // alpha applied globally (Left) or after trailing update (Right)
    return g;
}
#if defined(__CUDACC__) && !defined(MPLAPACK_MPFR_CUDA_HOST_EMULATION) && !defined(MPLAPACK_GEMM_CUDA_HOST_EMULATION)
template<class E,bool Checked> __global__ void k_trsm_scale(trsm_team_args<E> g,long i0,long j0,long rows,long cols,int *failed) {
    long u=(long)blockIdx.x*blockDim.x+threadIdx.x;
    if(u<rows*cols) trsm_scale_element(g,i0+u%rows,j0+u/rows);
    if constexpr(Checked) if(!g.e.ok) atomicExch(failed,1);
}
template<class E,bool Conjugate,bool Checked> __global__ void k_trsm_reciprocal(trsm_team_args<E> g,typename E::V *r,int *failed) {
    // Bounded physical indices also fit the resident complex MPC arena.
    for(long k=(long)blockIdx.x*blockDim.x+threadIdx.x;k<g.s.n;k+=(long)blockDim.x*gridDim.x)
        trsm_reciprocal_element<E,Conjugate>(g,r,k);
    if constexpr(Checked) if(!g.e.ok) atomicExch(failed,1);
}
template<class E,bool Conjugate,bool Checked> __global__ void k_trsm_update(trsm_team_args<E> g,
    long ioff,long joff,long rows,long cols,long start,long end,int *failed) {
    constexpr int T=trsm_tile_size<E>::value;
    __shared__ __align__(16) unsigned char storage[sizeof(trsm_update_shared<E>)];
    auto &sh=*reinterpret_cast<trsm_update_shared<E>*>(storage);
    const int ti=threadIdx.x,tj=threadIdx.y;
    const long ri=(long)blockIdx.x*T+ti,cj=(long)blockIdx.y*T+tj;
    const long i0=ioff+(long)blockIdx.x*T,j0=joff+(long)blockIdx.y*T;
    const bool inside=ri<rows && cj<cols;
    auto acc=g.e.reg(0);
    if(inside) g.e.set(acc,g.e.at(g.x.c,ioff+ri,joff+cj));
    for(long p0=0;p0<end-start;p0+=T) {
        trsm_update_load(g,sh,i0,j0,start,end,p0,ti,tj); __syncthreads();
        const long count=end-start-p0<T?end-start-p0:T;
        if(inside) trsm_update_step<E,Conjugate>(g,sh,acc,count,ti,tj);
        __syncthreads();
    }
    if(inside) g.e.put(g.x.c,ioff+ri,joff+cj,acc);
    if constexpr(Checked) if(!g.e.ok) atomicExch(failed,1);
}
template<class E,bool Conjugate=false,bool Checked=false> static int trsm_team_device(
    trsm_team_args<E> g,cudaStream_t stream=nullptr,int *failed=nullptr) {
    // Large value types cannot fit two tiles in ordinary static shared memory.
    if constexpr(sizeof(typename E::V)>1024)
        return trsm_team_plain_device<E,Conjugate,Checked>(g,stream,failed);
    else {
        if(!trsm_block_eligible(g.s)) return trsm_team_plain_device<E,Conjugate,Checked>(g,stream,failed);
        const long order=g.s.left?g.s.m:g.s.n;
        // Large solves use 32-wide diagonals; p64 stays inside one team.
        constexpr long NB=32;
        constexpr int T=trsm_tile_size<E>::value;
        typename E::V *reciprocal=nullptr;
        int error=0;
        if(!g.s.left && g.s.nounit) {
            error=cudaMalloc(&reciprocal,size_t(order)*sizeof(typename E::V));
            if(error) return error;
            const unsigned blocks=(unsigned)((order+63)/64<128?(order+63)/64:128);
            k_trsm_reciprocal<E,Conjugate,Checked><<<blocks,64,0,stream>>>(g,reciprocal,failed);
            error=cudaGetLastError(); g.reciprocal=reciprocal;
        }
        if(!error && order<=64) {
            error=trsm_team_plain_device<E,Conjugate,Checked>(g,stream,failed);
            if(reciprocal) {
                int e=cudaStreamSynchronize(stream); if(!error) error=e;
                e=cudaFree(reciprocal); if(!error) error=e;
            }
            return error;
        }
        if(!error && g.s.left && !g.s.alpha_one) {
            k_trsm_scale<E,Checked><<<(unsigned)((g.s.m*g.s.n+255)/256),256,0,stream>>>(g,0,0,g.s.m,g.s.n,failed);
            error=cudaGetLastError();
        }
        for(long done=0;done<order && !error;done+=NB) {
            const long count=order-done<NB?order-done:NB;
            const long start=g.s.upper?order-done-count:done,end=start+count;
            error=trsm_team_plain_device<E,Conjugate,Checked>(trsm_diagonal_args(g,start,count),stream,failed);
            const long remain=g.s.upper?start:order-end,offset=g.s.upper?0:end;
            const long rows=g.s.left?remain:g.s.m,cols=g.s.left?g.s.n:remain;
            if(!error && rows && cols) {
                k_trsm_update<E,Conjugate,Checked><<<dim3((unsigned)((rows+T-1)/T),(unsigned)((cols+T-1)/T)),dim3(T,T),0,stream>>>(
                    g,g.s.left?offset:0,g.s.left?0:offset,rows,cols,start,end,failed);
                error=cudaGetLastError();
            }
            if(!error && !g.s.left && !g.s.alpha_one) {
                k_trsm_scale<E,Checked><<<(unsigned)((g.s.m*count+255)/256),256,0,stream>>>(g,0,start,g.s.m,count,failed);
                error=cudaGetLastError();
            }
        }
        if(reciprocal) {
            int e=cudaStreamSynchronize(stream); if(!error) error=e;
            e=cudaFree(reciprocal); if(!error) error=e;
        }
        return error;
    }
}
#else
template<class E,bool Conjugate=false,bool Checked=false> static int trsm_team_device(trsm_team_args<E> g) {
    if constexpr(sizeof(typename E::V)>1024) return trsm_team_plain_device<E,Conjugate,Checked>(g);
    else {
        if(!trsm_block_eligible(g.s)) return trsm_team_plain_device<E,Conjugate,Checked>(g);
        const long order=g.s.left?g.s.m:g.s.n;
        std::vector<typename E::V> reciprocal;
        if(!g.s.left && g.s.nounit) {
            reciprocal.resize(order);
            for(long k=0;k<order;++k) trsm_reciprocal_element<E,Conjugate>(g,reciprocal.data(),k);
            g.reciprocal=reciprocal.data();
        }
        if(order<=64) {
            int error=trsm_team_plain_device<E,Conjugate,Checked>(g);
            if constexpr(Checked) if(!g.e.ok) return 1;
            return error;
        }
        if(g.s.left && !g.s.alpha_one) for(long j=0;j<g.s.n;++j) for(long i=0;i<g.s.m;++i) trsm_scale_element(g,i,j);
        for(long done=0;done<order;done+=32) {
            const long count=order-done<32?order-done:32;
            const long start=g.s.upper?order-done-count:done,end=start+count;
            int error=trsm_team_plain_device<E,Conjugate,Checked>(trsm_diagonal_args(g,start,count));
            if(error) return error;
            const long remain=g.s.upper?start:order-end,offset=g.s.upper?0:end;
            const long rows=g.s.left?remain:g.s.m,cols=g.s.left?g.s.n:remain;
            constexpr int T=trsm_tile_size<E>::value;
            for(long j0=0;j0<cols;j0+=T) for(long i0=0;i0<rows;i0+=T) {
                const long ib=i0+(g.s.left?offset:0),jb=j0+(g.s.left?0:offset);
                trsm_update_shared<E> sh;
                typename E::V acc[T][T];
                for(int tj=0;tj<T && j0+tj<cols;++tj) for(int ti=0;ti<T && i0+ti<rows;++ti)
                    g.e.set(acc[tj][ti],g.e.at(g.x.c,ib+ti,jb+tj));
                for(long p0=0;p0<count;p0+=T) {
                    for(int tj=0;tj<T;++tj) for(int ti=0;ti<T;++ti) trsm_update_load(g,sh,ib,jb,start,end,p0,ti,tj);
                    for(int tj=0;tj<T && j0+tj<cols;++tj) for(int ti=0;ti<T && i0+ti<rows;++ti)
                        trsm_update_step<E,Conjugate>(g,sh,acc[tj][ti],count-p0<T?count-p0:T,ti,tj);
                }
                for(int tj=0;tj<T && j0+tj<cols;++tj) for(int ti=0;ti<T && i0+ti<rows;++ti)
                    g.e.put(g.x.c,ib+ti,jb+tj,acc[tj][ti]);
            }
            if(!g.s.left && !g.s.alpha_one) for(long k=start;k<end;++k) for(long i=0;i<g.s.m;++i) trsm_scale_element(g,i,k);
        }
        if constexpr(Checked) return g.e.ok?0:1;
        return 0;
    }
}
#endif
}
#endif
