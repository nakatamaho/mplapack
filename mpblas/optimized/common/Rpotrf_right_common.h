// Value-semantic, ordered Lower Cholesky trailing updates. Each output uses
// its eventual reference block's SYRK/HERK or GEMM rule, not the current panel.
#ifndef MPLAPACK_POTRF_RIGHT_COMMON_H
#define MPLAPACK_POTRF_RIGHT_COMMON_H
#include <vector>
#if defined(__CUDACC__)
#include <cuda_runtime.h>
#define POTRF_RIGHT_HD __host__ __device__
#else
#define POTRF_RIGHT_HD
#endif
namespace mplapack_lapack_gpu {
template<class P> struct right_chol_job {
    using V=typename P::V; using M=typename P::M;
    P p; M a{},constants{};
    long n=0,nb=0,col=0,k=0,first=0,cols=0;
    static constexpr int T=sizeof(V)<=64?16:sizeof(V)<=256?8:sizeof(V)<=1024?4:2;
    struct shared { V a[T][T+1],b[T][T+1]; bool active[T][T+1]; };
    POTRF_RIGHT_HD bool mine(long i,long j) const { return i<n && j<first+cols && i>=j; }
    POTRF_RIGHT_HD void init(V &acc,long i,long j) {
        if(mine(i,j)) p.init(acc,p.at(a,i,j),i==j);
    }
    POTRF_RIGHT_HD void load(shared &sh,long i0,long j0,long q,int x,int y) {
        if(q+y>=k) return;
        if(i0+x<n) sh.a[y][x]=p.at(a,i0+x,col+q+y);
        if(j0+x<first+cols) {
            V b=p.at(a,j0+x,col+q+y);
            sh.active[y][x]=!p.iszero(b);
            p.scaled(sh.b[y][x],b,i0/nb==j0/nb,constants);
        }
    }
    POTRF_RIGHT_HD void step(const shared &sh,V &acc,long i,long j,long count,int x,int y) {
        if(!mine(i,j)) return;
        const bool rank=i/nb==j/nb;
        V product=p.reg(2);
        for(long q=0;q<count;++q) if(!rank || sh.active[q][y]) {
            p.product(product,sh.b[q][y],sh.a[q][x]);
            p.accumulate(acc,product,i==j);
        }
    }
    POTRF_RIGHT_HD void finish(const V &acc,long i,long j) { if(mine(i,j)) p.put(a,i,j,acc); }
    POTRF_RIGHT_HD void operator()(long lane) {
        const long rows=n-first,i=first+lane%rows,j=first+lane/rows;
        if(!mine(i,j)) return;
        V acc=p.reg(0),temp=p.reg(1),product=p.reg(2);
        init(acc,i,j);
        const bool rank=i/nb==j/nb;
        for(long q=0;q<k;++q) {
            V b=p.at(a,j,col+q);
            if(rank && p.iszero(b)) continue;
            p.scaled(temp,b,rank,constants);
            p.product(product,temp,p.at(a,i,col+q));
            p.accumulate(acc,product,i==j);
        }
        finish(acc,i,j);
    }
};
#if defined(__CUDACC__)
template<class J> __global__ void k_right_chol(J j,int *failed) {
    constexpr int T=J::T;
    const long i0=j.first+(long)blockIdx.x*T,j0=j.first+(long)blockIdx.y*T;
    if(i0+T<=j0) return;
    __shared__ __align__(16) unsigned char storage[sizeof(typename J::shared)];
    auto &sh=*reinterpret_cast<typename J::shared*>(storage);
    const int x=threadIdx.x,y=threadIdx.y;
    const long i=i0+x,col=j0+y;
    typename J::V acc=j.p.reg(0);
    j.init(acc,i,col);
    for(long q=0;q<j.k;q+=T) {
        j.load(sh,i0,j0,q,x,y); __syncthreads();
        j.step(sh,acc,i,col,j.k-q<T?j.k-q:T,x,y); __syncthreads();
    }
    j.finish(acc,i,col);
    if(failed && !j.p.ok()) atomicExch(failed,1);
}
template<class J> __global__ void k_right_chol_scalar(J j,int *failed) {
    const long u=(long)blockIdx.x*blockDim.x+threadIdx.x;
    if(u<(j.n-j.first)*j.cols) j(u);
    if(failed && !j.p.ok()) atomicExch(failed,1);
}
template<class P> static int right_chol_device(right_chol_job<P> j,cudaStream_t stream,int *failed=nullptr) {
    constexpr int T=right_chol_job<P>::T;
    if(j.nb%T) k_right_chol_scalar<<<(unsigned)(((j.n-j.first)*j.cols+63)/64),64,0,stream>>>(j,failed);
    else k_right_chol<<<dim3((unsigned)((j.n-j.first+T-1)/T),(unsigned)((j.cols+T-1)/T)),dim3(T,T),0,stream>>>(j,failed);
    return cudaGetLastError();
}
#else
template<class P> static int right_chol_host(right_chol_job<P> j) {
    using J=right_chol_job<P>; constexpr int T=J::T;
    if(j.nb%T) {
        for(long u=0;u<(j.n-j.first)*j.cols;++u) j(u);
    } else for(long col=j.first;col<j.first+j.cols;col+=T) for(long row=col;row<j.n;row+=T) {
        typename J::shared sh;
        typename J::V acc[T][T];
        for(int y=0;y<T;++y) for(int x=0;x<T;++x) j.init(acc[y][x],row+x,col+y);
        for(long q=0;q<j.k;q+=T) {
            for(int y=0;y<T;++y) for(int x=0;x<T;++x) j.load(sh,row,col,q,x,y);
            for(int y=0;y<T;++y) for(int x=0;x<T;++x) j.step(sh,acc[y][x],row+x,col+y,j.k-q<T?j.k-q:T,x,y);
        }
        for(int y=0;y<T;++y) for(int x=0;x<T;++x) j.finish(acc[y][x],row+x,col+y);
    }
    return j.p.ok()?0:1;
}
#endif
}
#undef POTRF_RIGHT_HD
#endif
