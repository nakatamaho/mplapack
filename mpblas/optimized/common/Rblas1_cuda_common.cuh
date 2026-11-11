// BF CUDA reduction driver: staged publication; no host numerical work.
#ifndef MPLAPACK_BLAS1_CUDA_COMMON_CUH
#define MPLAPACK_BLAS1_CUDA_COMMON_CUH
#include "Rblas_dev_types.h"
#include "Rlapack_sqrt_common.h"
#include "Rblas1_jobs_common.h"
#include <vector>
#include <algorithm>
#include <type_traits>
#if defined(__CUDACC__) && !defined(MPLAPACK_GEMM_CUDA_HOST_EMULATION)
#include <cuda_runtime.h>
#endif
namespace mplapack_blas1 __attribute__((visibility("hidden"))) {

// ---- one-pass iamax (factorization goal 13) ----------------------------
// A (score, index) partial; the winner is the maximum score with the smallest
// index on ties.  Scores are exact and comparisons are exact, so this rule is
// associative and commutative: any combination tree gives the reference scan
// result (strict > keeps the first index).
template<class R> struct iamax_pair { R score; long index; };
template<class X,class R,bool Complex> MPLAPACK_DEV_HD R iamax_score(const X &v) {
    if constexpr(Complex) return mplapack_blas_dev::RCabs1(v);
    else return abs(v);
}
template<class R> MPLAPACK_DEV_HD bool iamax_better(const iamax_pair<R> &a,const iamax_pair<R> &b) {
    return b.score>a.score || (b.score==a.score && b.index<a.index);
}
constexpr int iamax_threads=256,iamax_max_blocks=2048,iamax_elems_thread=16;
#if defined(__CUDACC__) && !defined(MPLAPACK_GEMM_CUDA_HOST_EMULATION)
template<class X,class R,bool Complex> __global__ void k_iamax_pass(const X *x,long n,iamax_pair<R> *out) {
    __shared__ iamax_pair<R> sh[iamax_threads];
    const int tid=threadIdx.x;
    iamax_pair<R> best;best.score=R(0.0);best.index=0;
    bool visited=false;
    for(long i=(long)blockIdx.x*iamax_threads+tid;i<n;i+=(long)gridDim.x*iamax_threads) {
        const R s=iamax_score<X,R,Complex>(x[i]);
        if(!visited || s>best.score) { best.score=s;best.index=i;visited=true; }
    }
    sh[tid]=best;
    __syncthreads();
    for(int stride=iamax_threads/2;stride>0;stride>>=1) {
        if(tid<stride && iamax_better<R>(sh[tid],sh[tid+stride])) sh[tid]=sh[tid+stride];
        __syncthreads();
    }
    if(tid==0) out[blockIdx.x]=sh[0];
}
template<class R> __global__ void k_iamax_final(const iamax_pair<R> *in,long nb,iamax_pair<R> *out) {
    __shared__ iamax_pair<R> sh[iamax_threads];
    const int tid=threadIdx.x;
    iamax_pair<R> best;best.score=R(0.0);best.index=0;
    bool visited=false;
    for(long i=tid;i<nb;i+=iamax_threads) {
        if(!visited || iamax_better<R>(best,in[i])) { best=in[i];visited=true; }
    }
    sh[tid]=best;
    __syncthreads();
    for(int stride=iamax_threads/2;stride>0;stride>>=1) {
        if(tid<stride && iamax_better<R>(sh[tid],sh[tid+stride])) sh[tid]=sh[tid+stride];
        __syncthreads();
    }
    if(tid==0) *out=sh[0];
}
#endif
#if defined(__CUDACC__) && !defined(MPLAPACK_GEMM_CUDA_HOST_EMULATION)
template<class J> __global__ void kernel(J job,long count) {
    for(long i=long(blockIdx.x)*blockDim.x+threadIdx.x;i<count;i+=long(blockDim.x)*gridDim.x) job(i);
}
// Goal 13 (continued) wavefront, grid version: the elements are spread over
// as many resident blocks as the device holds (cooperative launch), so the
// n^2/2 term updates are not confined to one SM; the columns still finalize
// one after another through per-element progress flags in global memory
// (2 ints per element: flag + owner-private term count, pre-zeroed).
template<class J> __global__ void k_wave_grid(J job,int *buf) {
    const long gtid=(long)blockIdx.x*blockDim.x+threadIdx.x;
    const long gT=(long)gridDim.x*blockDim.x;
    job.block_run(gtid,gT,buf,buf+job.s.n);
}
template<class V> struct buffer {
    V *p=nullptr;
    ~buffer() { cudaFree(p); }
    int allocate(size_t n) { return cudaMalloc(&p,n*sizeof(V)); }
    int up(const V *h,size_t n) { return cudaMemcpy(p,h,n*sizeof(V),cudaMemcpyHostToDevice); }
    int down(V *h,size_t n) { return cudaMemcpy(h,p,n*sizeof(V),cudaMemcpyDeviceToHost); }
};

struct launcher {
    int error=0;
    template<class J> void launch(J job,long n) {
        if(error) return;
        long blocks=(n+127)/128;if(blocks>65535) blocks=65535;
        kernel<<<unsigned(blocks),128>>>(job,n);
        error=cudaGetLastError();
    }
    template<class J> void launch_wave(J job) {
        if(error) return;
        const long n=job.s.n;
        if(n<=0) return;
        const int T=256;
        int perSM=0,sms=0;
        if(cudaOccupancyMaxActiveBlocksPerMultiprocessor(&perSM,(const void *)k_wave_grid<J>,T,0)!=cudaSuccess || perSM<=0)
            perSM=1;
        if(cudaDeviceGetAttribute(&sms,cudaDevAttrMultiProcessorCount,0)!=cudaSuccess || sms<=0)
            sms=1;
        long cap=(long)perSM*sms; // resident blocks (cooperative launch limit)
        long blocks=cap;         // fill the device: consecutive elements land on
                                 // consecutive threads across all blocks
        buffer<int> fb;
        if(fb.allocate(size_t(n)*2) || cudaMemset(fb.p,0,size_t(n)*2*sizeof(int))) { error=1;return; }
        J j2=job;
        void *args[2]={&j2,&fb.p};
        error=cudaLaunchCooperativeKernel((const void *)k_wave_grid<J>,dim3(unsigned(blocks)),dim3(T),args,0,0);
        if(!error) error=cudaDeviceSynchronize();
    }
};
#else
struct launcher {
    int error=0;
    template<class J> void launch(J job,long n) { for(long i=0;i<n;++i) job(i); }
    template<class J> void launch_wave(J job) {
        std::vector<int> f(size_t(job.s.n>0?job.s.n:0)*2,0);
        job.block_run(0,1,f.data(),f.data()+ (job.s.n>0?job.s.n:0));
    }
};
template<class V> struct buffer {
    std::vector<V> v;V *p=nullptr;
    int allocate(size_t n) { v.resize(n);p=v.data();return 0; }
    int up(const V *h,size_t n) { std::copy(h,h+n,p);return 0; }
    int down(V *h,size_t n) { std::copy(p,p+n,h);return 0; }
};
#endif
template<class Ops> struct root {
    using R=mplapack_blas_dev::dreal<Ops>;
    MPLAPACK_DEV_HD static R sqrt(const R &x) { return R(mplapack_lapack_gpu::square_root<Ops>(x.v)); }
};
template<class X,class R,bool Complex> int iamax_kernels(long n,const X *dx,iamax_pair<R> *part,iamax_pair<R> *out) {
#if defined(__CUDACC__) && !defined(MPLAPACK_GEMM_CUDA_HOST_EMULATION)
    const long blocks0=(n+iamax_elems_thread*iamax_threads-1)/(iamax_elems_thread*iamax_threads);
    const long nb=blocks0<iamax_max_blocks?blocks0:iamax_max_blocks;
    k_iamax_pass<X,R,Complex><<<unsigned(nb),iamax_threads>>>(dx,n,part);
    int error=cudaGetLastError();
    if(!error) { k_iamax_final<R><<<1,iamax_threads>>>(part,nb,out); error=cudaGetLastError(); }
    if(!error) error=cudaDeviceSynchronize();
    return error;
#else
    (void)part;
    iamax_pair<R> best;best.score=R(0.0);best.index=0;
    bool visited=false;
    for(long i=0;i<n;++i) {
        const R s=iamax_score<X,R,Complex>(dx[i]);
        if(!visited || s>best.score) { best.score=s;best.index=i;visited=true; }
    }
    *out=best;
    return 0;
#endif
}
// One-pass iamax (factorization goal 13): one kernel over x reduces
// (score, index) in registers and shared memory to one winner per block; a
// single small second kernel reduces the block winners.  No O(n) workspace,
// no per-level launches.  The min-index tie rule is order independent, so the
// result equals the reference scan bit for bit.  Emulation runs the same
// staged data with one serial scan (same exact scores, same first-tie rule).
template<class Ops,bool Complex> int run_iamax(long n,
    const typename std::conditional<Complex,typename mplapack_blas_dev::cgpu_types<Ops>::DC,typename mplapack_blas_dev::cgpu_types<Ops>::DR>::type *x,
    state<typename mplapack_blas_dev::cgpu_types<Ops>::DR> &result) {
    using R=typename mplapack_blas_dev::cgpu_types<Ops>::DR;
    using X=typename std::conditional<Complex,typename mplapack_blas_dev::cgpu_types<Ops>::DC,R>::type;
    using P=iamax_pair<R>;
#if defined(__CUDACC__) && !defined(MPLAPACK_GEMM_CUDA_HOST_EMULATION)
    size_t free=0,total=0;
    if(cudaMemGetInfo(&free,&total) || free/sizeof(X)<size_t(n)) return 1;
    const long blocks0=(n+iamax_elems_thread*iamax_threads-1)/(iamax_elems_thread*iamax_threads);
    const long nb=blocks0<iamax_max_blocks?blocks0:iamax_max_blocks;
    if(free<size_t(n)*sizeof(X)+size_t(nb+1)*sizeof(P)) return 1;
#else
    const long nb=0;
#endif
    buffer<X> dx;buffer<P> part,res;
    if(dx.allocate(n) || (nb && part.allocate(nb)) || res.allocate(1) || dx.up(x,n)) return 1;
    int error=iamax_kernels<X,R,Complex>(n,dx.p,part.p,res.p);
    if(error) return error;
    P r;
    if(res.down(&r,1)) return 1;
    result.re=r.score;result.index=r.index;
    return 0;
}
template<class Ops,bool Complex> int run_bf(operation op,long n,
    const typename std::conditional<Complex,typename mplapack_blas_dev::cgpu_types<Ops>::DC,typename mplapack_blas_dev::cgpu_types<Ops>::DR>::type *x,
    const typename std::conditional<Complex,typename mplapack_blas_dev::cgpu_types<Ops>::DC,typename mplapack_blas_dev::cgpu_types<Ops>::DR>::type *y,
    state<typename mplapack_blas_dev::cgpu_types<Ops>::DR> &result) {
    using R=typename mplapack_blas_dev::cgpu_types<Ops>::DR;
    using X=typename std::conditional<Complex,typename mplapack_blas_dev::cgpu_types<Ops>::DC,R>::type;
    if(op==mplapack_blas1::iamax) return run_iamax<Ops,Complex>(n,x,result);
    long count=Complex && op==norm?2*n:n;
    const bool product=op==dot || op==dotc || op==dotu;
    // goal 13 fast schedule workspace: sub-lane plus chunk states
    const long ws=count>0?std::max(count,(count+fast_sub_elements-1)/fast_sub_elements+(count+fast_chunk_length-1)/fast_chunk_length):0;
#if defined(__CUDACC__) && !defined(MPLAPACK_GEMM_CUDA_HOST_EMULATION)
    size_t free=0,total=0;
    if(cudaMemGetInfo(&free,&total) || free/sizeof(state<R>)<size_t(ws) ||
       free-size_t(ws)*sizeof(state<R>)<size_t(n)*sizeof(X)*(product?2:1)) return 1;
#endif
    buffer<X> dx,dy;buffer<state<R>> a;
    if(dx.allocate(n) || a.allocate(ws) || dx.up(x,n) || (product && (dy.allocate(n) || dy.up(y,n)))) return 1;
    launcher ctx;
    schedule<R,X,Complex,root<Ops>>(ctx,op,n,dx.p,dy.p,a.p);
    if(ctx.error) return ctx.error;
    return a.down(&result,1);
}
}
#endif
