// Value-semantic MPFR GEMM: shared operands, unchanged per-element rounding.
#ifndef MPLAPACK_MPFR_GEMM_TILED_CUDA_H
#define MPLAPACK_MPFR_GEMM_TILED_CUDA_H
#include "Rgemm_kernel_cuda.h"
#if defined(__CUDACC__)
#include <cuda_runtime.h>
#endif
namespace mplapack_mpfr_cuda {
namespace {
template<int PB> struct value_gemm_public_engine : gemm_fixed_engine<PB> {
    static constexpr bool ok=true;
};
template<class E,bool Nota,bool Packed> struct value_gemm_job {
    using V=typename E::V; using M=typename E::M;
    static constexpr int MR=8, NR=sizeof(V)<144?16:8, KR=8, threads=MR*NR;
    // At 1024 bits the unconstrained kernel remains register-limited despite
    // fewer operand loads. More resident warps outweigh the extra spills.
    static constexpr int min_blocks=sizeof(V)==144?8:1;
    struct shared { V a[KR][MR+1],b[KR][NR+1]; };
    E e; gemm_shape s; M al,be,a,b,ab,c;
    __host__ __device__ void bind() {}
    __host__ __device__ V base(long i,long j) {
        V v=e.reg(0);
        if(s.beta_zero) e.zero(v);
        else if(s.beta_one) e.set(v,e.at(c,i,j));
        else e.mul(v,e.at(be,0,0),e.at(c,i,j));
        return v;
    }
    __host__ __device__ void init(V &v,long i,long j) {
        if(i>=s.m || j>=s.n) return;
        if constexpr(Nota) v=base(i,j);
        else {
            // Preserve beta multiplication before the dot, but keep its exact
            // packed value in private C rather than live across the whole loop.
            if(!s.beta_zero && !s.beta_one) { V z=base(i,j);e.put(c,i,j,z); }
            e.zero(v);
        }
    }
    __host__ __device__ void load(shared &sh,long i0,long j0,long q,int t) {
        for(int u=t;u<KR*MR;u+=threads) {
            int x=Nota?u%MR:u/KR, l=Nota?u/MR:u%KR;
            if(i0+x<s.m && q+l<s.k)
                e.set(sh.a[l][x],Nota?e.at(a,i0+x,q+l):e.at(a,q+l,i0+x));
        }
        for(int u=t;u<KR*NR;u+=threads) {
            const bool kfast=Packed || s.notb;
            int y=kfast?u/KR:u%NR, l=kfast?u%KR:u/NR;
            if(j0+y>=s.n || q+l>=s.k) continue;
            if constexpr(Packed) e.set(sh.b[l][y],e.at(ab,q+l,j0+y));
            else {
                V v=s.notb?e.at(b,q+l,j0+y):e.at(b,j0+y,q+l);
                if constexpr(Nota) e.mul(sh.b[l][y],e.at(al,0,0),v);
                else e.set(sh.b[l][y],v);
            }
        }
    }
    __host__ __device__ void step(const shared &sh,V &v,long i,long j,int count,int x,int y) {
        if(i>=s.m || j>=s.n) return;
        for(int l=0;l<count;++l) {
            V z=e.reg(1);
            // No GEMM zero skip, even for an exactly zero original operand.
            if constexpr(Nota) e.mul(z,sh.b[l][y],sh.a[l][x]);
            else e.mul(z,sh.a[l][x],sh.b[l][y]);
            e.add(v,v,z);
        }
    }
    __host__ __device__ void finish(V &v,long i,long j) {
        if(i>=s.m || j>=s.n) return;
        if constexpr(!Nota) {
            V z=e.reg(1),c0=e.reg(2);
            e.mul(z,e.at(al,0,0),v);
            if(s.beta_zero) e.zero(c0); else e.set(c0,e.at(c,i,j));
            e.add(v,c0,z);
        }
        e.put(c,i,j,v);
    }
};
template<int PB,bool Nota> struct value_gemm_public_job : value_gemm_job<value_gemm_public_engine<PB>,Nota,Nota> {
    using V=real_t<PB>; V alpha,beta;
    __host__ __device__ void bind() { this->al={&alpha,1};this->be={&beta,1}; }
};
template<int PB,bool Nota> static value_gemm_public_job<PB,Nota> make_value_gemm_job(
    gemm_shape s,const real_t<PB> &alpha,const real_t<PB> &beta,
    const real_t<PB> *a,const real_t<PB> *b,const real_t<PB> *ab,real_t<PB> *c) {
    value_gemm_public_job<PB,Nota> j{};j.s=s;j.alpha=alpha;j.beta=beta;
    j.a={a,s.lda};j.b={b,s.ldb};j.ab={ab,s.k};j.c={c,s.ldc};return j;
}
#if defined(__CUDACC__)
template<class J> __device__ __forceinline__ void value_gemm_body(J j,int *failed) {
    j.bind();
    __shared__ __align__(16) unsigned char storage[sizeof(typename J::shared)];
    auto &sh=*reinterpret_cast<typename J::shared*>(storage);
    const int x=J::min_blocks>1?threadIdx.x%J::MR:threadIdx.x;
    const int y=J::min_blocks>1?threadIdx.x/J::MR:threadIdx.y;
    const int t=x+J::MR*y;
    const long i0=(long)blockIdx.x*J::MR,i=i0+x;
    for(long j0=(long)blockIdx.y*J::NR;j0<j.s.n;j0+=(long)gridDim.y*J::NR) {
        typename J::V v=j.e.reg(0);j.init(v,i,j0+y);
        for(long q=0;q<j.s.k;q+=J::KR) {
            j.load(sh,i0,j0,q,t);__syncthreads();
            j.step(sh,v,i,j0+y,j.s.k-q<J::KR?j.s.k-q:J::KR,x,y);__syncthreads();
        }
        j.finish(v,i,j0+y);
        __syncthreads();
    }
    if(failed && !j.e.ok) atomicExch(failed,1);
}
template<class J> __global__ void value_gemm_kernel(J j,int *failed) { value_gemm_body(j,failed); }
template<class J> __global__ __launch_bounds__(J::threads,J::min_blocks)
void value_gemm_bounded_kernel(J j,int *failed) { value_gemm_body(j,failed); }
template<class J> static int value_gemm_device(J j,cudaStream_t stream,int *failed=nullptr) {
    long gy=(j.s.n+J::NR-1)/J::NR;if(gy>65535) gy=65535;
    dim3 grid((unsigned)((j.s.m+J::MR-1)/J::MR),(unsigned)gy);
    // Even a one-block launch bound changes register-allocation heuristics.
    // Apply a bound only to the precision for which it was measured to help.
    if constexpr(J::min_blocks>1) value_gemm_bounded_kernel<<<grid,dim3(J::threads),0,stream>>>(j,failed);
    else value_gemm_kernel<<<grid,dim3(J::MR,J::NR),0,stream>>>(j,failed);
    return cudaGetLastError();
}
#else
template<class J> static int value_gemm_host(J j) {
    j.bind();
    for(long j0=0;j0<j.s.n;j0+=J::NR) for(long i0=0;i0<j.s.m;i0+=J::MR) {
        typename J::shared sh;
        typename J::V v[J::NR][J::MR];
        for(int y=0;y<J::NR;++y) for(int x=0;x<J::MR;++x) j.init(v[y][x],i0+x,j0+y);
        for(long q=0;q<j.s.k;q+=J::KR) {
            for(int t=0;t<J::threads;++t) j.load(sh,i0,j0,q,t);
            for(int y=0;y<J::NR;++y) for(int x=0;x<J::MR;++x)
                j.step(sh,v[y][x],i0+x,j0+y,j.s.k-q<J::KR?j.s.k-q:J::KR,x,y);
        }
        for(int y=0;y<J::NR;++y) for(int x=0;x<J::MR;++x) j.finish(v[y][x],i0+x,j0+y);
    }
    return j.e.ok?0:1;
}
#endif
} // private implementation
} // namespace mplapack_mpfr_cuda
#endif
