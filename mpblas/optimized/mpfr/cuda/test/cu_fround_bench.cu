// Device throughput (kernel-only ns/operation, not single-thread latency).
// Build with nvcc -O3 -rdc=true -fmad=false, the patched include directory,
// libmpc_cuda.a and cudart. Run with optional pair count (default 65536).
#include <mpc_cuda.cuh>
#include <cuda_runtime.h>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <vector>
using namespace cu_fp;
#define CHECK(x) do { cudaError_t e=(x); if(e!=cudaSuccess) { \
 std::fprintf(stderr,"%s: %s\n",#x,cudaGetErrorString(e)); std::exit(1); } } while(0)
template<int P> struct pair_t { cu_freal<P> a,b; };
template<int P,int OP> __global__ void fixed_op(const pair_t<P>* in,cu_freal<P>* out,int n)
{
    for(int i=blockIdx.x*blockDim.x+threadIdx.x;i<n;i+=gridDim.x*blockDim.x) {
        auto a=in[i].a,b=in[i].b; cu_freal<P> r;
        if(OP==0) r=cu_fdiv<P>(a,b);
        if(OP==1) cu_fdiv_rn<P>(r,a,b);
        if(OP==2) r=cu_sqrt<P>(a);
        if(OP==3) cu_sqrt_rn<P>(r,a);
        out[i]=r;
    }
}
template<int P,bool SQRT> __global__ void runtime_op(const pair_t<P>* in,cu_freal<P>* out,int n)
{
    const int N=cu_freal<P>::N;
    for(int i=blockIdx.x*blockDim.x+threadIdx.x;i<n;i+=gridDim.x*blockDim.x) {
        mpc_cuda_arena_reset();
        cu_mp_limb_t am[N],bm[N],rm[N];
        for(int k=0;k<N;++k) { am[k]=in[i].a.m[k]; bm[k]=in[i].b.m[k]; }
        cu_mpfr_t a,b,r;
        cu_mpfr_custom_init_set(a,CU_MPFR_REGULAR_KIND,in[i].a.exp,P,am);
        cu_mpfr_custom_init_set(b,CU_MPFR_REGULAR_KIND,in[i].b.exp,P,bm);
        cu_mpfr_custom_init_set(r,CU_MPFR_ZERO_KIND,0,P,rm);
        if(SQRT) cu_mpfr_sqrt(r,a,CU_MPFR_RNDN);
        else cu_mpfr_div(r,a,b,CU_MPFR_RNDN);
        cu_freal<P> v; v.sign=1; v.exp=cu_mpfr_get_exp(r);
        for(int k=0;k<N;++k) v.m[k]=rm[k];
        out[i]=v;
    }
}
template<int P> void run(int n)
{
    const int blocks=128,threads=64;
    pair_t<P>* in; cu_freal<P>* out;
    CHECK(cudaMalloc(&in,n*sizeof(*in))); CHECK(cudaMalloc(&out,n*sizeof(*out)));
    std::vector<pair_t<P>> h(n); std::mt19937_64 rng(0x07a+P);
    for(auto& v:h) {
        v.a.sign=v.b.sign=1; v.a.exp=v.b.exp=1;
        for(int k=0;k<cu_freal<P>::N;++k) { v.a.m[k]=rng(); v.b.m[k]=rng(); }
        v.a.m[cu_freal<P>::N-1]|=1ULL<<63; v.b.m[cu_freal<P>::N-1]|=1ULL<<63;
    }
    CHECK(cudaMemcpy(in,h.data(),n*sizeof(*in),cudaMemcpyHostToDevice));
    cudaEvent_t start,stop; CHECK(cudaEventCreate(&start)); CHECK(cudaEventCreate(&stop));
    std::vector<double> times[6];
    for(int round=-1;round<5;++round) for(int k=0;k<6;++k) {
        int op=(round&1)?5-k:k; // alternate order, discard warmup
        CHECK(cudaEventRecord(start));
        switch(op) {
        case 0: fixed_op<P,0><<<blocks,threads>>>(in,out,n); break;
        case 1: fixed_op<P,1><<<blocks,threads>>>(in,out,n); break;
        case 2: fixed_op<P,2><<<blocks,threads>>>(in,out,n); break;
        case 3: fixed_op<P,3><<<blocks,threads>>>(in,out,n); break;
        case 4: runtime_op<P,false><<<blocks,threads>>>(in,out,n); break;
        case 5: runtime_op<P,true><<<blocks,threads>>>(in,out,n); break;
        }
        CHECK(cudaGetLastError()); CHECK(cudaEventRecord(stop)); CHECK(cudaEventSynchronize(stop));
        float ms; CHECK(cudaEventElapsedTime(&ms,start,stop));
        if(round>=0) times[op].push_back(double(ms)*1e6/n);
    }
    const char* names[]={"cu_fdiv","cu_fdiv_rn","cu_sqrt","cu_sqrt_rn","cu_mpfr_div","cu_mpfr_sqrt"};
    for(int op=0;op<6;++op) { std::sort(times[op].begin(),times[op].end());
        std::printf("PB=%d %s ns/op=%.3f median_of=5 n=%d\n",P,names[op],times[op][2],n); }
    CHECK(cudaEventDestroy(start)); CHECK(cudaEventDestroy(stop)); CHECK(cudaFree(in)); CHECK(cudaFree(out));
}
int main(int argc,char** argv)
{
    int dev=0; if(cudaGetDeviceCount(&dev)!=cudaSuccess || !dev) return 77;
    CHECK(cudaDeviceSetLimit(cudaLimitStackSize,65536));
    CHECK(cudaDeviceSetLimit(cudaLimitMallocHeapSize,64*1024*1024));
    const size_t slab=32768, pool=128*64;
    char* arena; size_t* top;
    CHECK(cudaMalloc(&arena,pool*slab)); CHECK(cudaMalloc(&top,pool*sizeof(size_t)));
    CHECK(cudaMemset(top,0,pool*sizeof(size_t)));
    mpc_cuda_arena_base=arena; mpc_cuda_arena_slab=slab; mpc_cuda_arena_top=top;
    int n=argc>1?std::atoi(argv[1]):65536;
    if(n<=0) return 1;
    run<256>(n); run<512>(n); run<1024>(n);
    mpc_cuda_arena_base=nullptr; mpc_cuda_arena_top=nullptr; mpc_cuda_arena_slab=0;
    CHECK(cudaFree(arena)); CHECK(cudaFree(top));
}
