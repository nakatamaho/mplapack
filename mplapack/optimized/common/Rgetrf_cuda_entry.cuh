// dd/td/qd LU public entries. No CPU fallback; publish only after all downloads.
#if defined(__CUDACC__)
#define MPLAPACK_MINMAX_MPLAPACKINT_DEFINED
#endif
#include <mpblas.h>
#include <mplapack.h>
#include <cstdlib>
#include <cstdio>
#include <chrono>
#include <vector>
#include <limits>
#include <cstring>
#include <stdexcept>
#if defined(MPLAPACK_BUILD_WITH_DD)
#include <mpblas_dd_cuda.h>
#elif defined(MPLAPACK_BUILD_WITH_TD)
#include <mpblas_td_cuda.h>
#else
#include <mpblas_qd_cuda.h>
#endif
#pragma GCC visibility push(hidden)
#include "Rgetrf_gpu_common.cuh"
#pragma GCC visibility pop
namespace {
#if defined(MPLAPACK_BUILD_WITH_DD)
using LuOps=mplapack_gemm::dd_ops;
#define LU_ENTRY(R) R##_dd_cuda
static constexpr const char *lu_prefix="MPLAPACK_DD";
static constexpr long lu_min=2048; // Measured real/complex resident LU crossover.
#elif defined(MPLAPACK_BUILD_WITH_TD)
using LuOps=mplapack_gemm::td_ops;
#define LU_ENTRY(R) R##_td_cuda
static constexpr const char *lu_prefix="MPLAPACK_TD";
static constexpr long lu_min=1024; // Measured for both Rgetrf and Cgetrf.
#else
using LuOps=mplapack_gemm::qd_ops;
#define LU_ENTRY(R) R##_qd_cuda
static constexpr const char *lu_prefix="MPLAPACK_QD";
static constexpr long lu_min=1024; // Measured for both Rgetrf and Cgetrf.
#endif
using LuReal=mplapack_blas_dev::dreal<LuOps>;
using LuComplex=mplapack_blas_dev::dcomplex<LuReal>;
static LuReal lu_pack(const REAL &a) {
    LuReal r;
    for(unsigned q=0;q<sizeof(r.v.x)/sizeof(double);++q) r.v.x[q]=a.x[q];
    return r;
}
static LuComplex lu_pack(const COMPLEX &a) { return LuComplex(lu_pack(a.real()),lu_pack(a.imag())); }
static void lu_unpack(const LuReal &r,REAL &a) {
    for(unsigned q=0;q<sizeof(r.v.x)/sizeof(double);++q) a.x[q]=r.v.x[q];
}
static void lu_unpack(const LuComplex &r,COMPLEX &a) {
    REAL re,im; lu_unpack(r.real(),re); lu_unpack(r.imag(),im); a=COMPLEX(re,im);
}
static const char *lu_env(const char *suffix) {
    char key[96]; std::snprintf(key,sizeof(key),"%s%s",lu_prefix,suffix); return std::getenv(key);
}
// Read on every call: LAPACK_CUDA=0 disables; LAPACK_CUDA_MIN controls
// min(m,n); CUDA_VERBOSE prints measured transfer/factor times. LAPACK_FAST=1
// currently retains the reference algorithm (no implicit sloppy arithmetic).
template<class Blas,class Host> static bool lu_call(const char *name,INTEGER m,INTEGER n,Host *a,INTEGER lda,INTEGER *ipiv,INTEGER &info) {
    using namespace mplapack_lapack_gpu;
    using E=typename Blas::E;
    if(m<=0 || n<=0 || lda<m || !a || !ipiv) return false;
    const char *enable=lu_env("_LAPACK_CUDA");
    if(enable && !std::strcmp(enable,"0")) return false;
    long minimum=lu_min;
    if(const char *threshold=lu_env("_LAPACK_CUDA_MIN")) {
        char *end=nullptr; long v=std::strtol(threshold,&end,10);
        if(end!=threshold && !*end && v>=0) minimum=v;
    }
    if(std::min(m,n)<minimum || (size_t)n>std::numeric_limits<size_t>::max()/(size_t)m/sizeof(E)) return false;
    // Also reject an unrepresentable host address span before accessing A.
    if((size_t)n>std::numeric_limits<size_t>::max()/(size_t)lda/sizeof(Host)) return false;
    try {
        long nb=iMlaenv(1,name," ",m,n,-1,-1);
        long workspace=std::max(1L,std::min((long)std::min(m,n),nb));
        matrix<E> device;
        if(device.allocate(m,n)) return false;
        lu_context<Blas,LuReal> ctx(std::max(m,n),workspace,lu_pack(Rlamch("S")));
        buffer<long> piv;
        if(ctx.error || piv.allocate(std::min(m,n))) return false;
        std::vector<E> staging((size_t)m*n);
        std::vector<long> host_piv(std::min(m,n));
        using Clock=std::chrono::steady_clock;
        auto start=Clock::now();
        for(long j=0;j<n;++j) for(long i=0;i<m;++i) staging[i+j*m]=lu_pack(a[i+j*lda]);
        if(device.copy(staging.data(),m,false)) return false;
#if defined(__CUDACC__)
        // Pageable H2D copies may return after staging, before DMA completes.
        // The two nonblocking streams must not race the default-stream upload.
        if(cudaStreamSynchronize(nullptr)!=cudaSuccess) return false;
#endif
        auto uploaded=Clock::now();
        long result=lu(ctx,m,n,nb,device.data(),device.ld,piv.data());
        if(ctx.error) return false;
        auto factored=Clock::now();
        if(device.copy(staging.data(),m,true) || transfer(host_piv.data(),piv.data(),host_piv.size(),true)) return false;
        for(long j=0;j<n;++j) for(long i=0;i<m;++i) lu_unpack(staging[i+j*m],a[i+j*lda]);
        for(size_t i=0;i<host_piv.size();++i) ipiv[i]=host_piv[i];
        info=result;
        auto downloaded=Clock::now();
        const char *verbose=lu_env("_CUDA_VERBOSE");
        if(verbose && std::strcmp(verbose,"0"))
            std::fprintf(stderr,"%s %s CUDA m=%ld n=%ld upload=%.9g factor=%.9g download=%.9g\n",name,lu_prefix,(long)m,(long)n,
                std::chrono::duration<double>(uploaded-start).count(),std::chrono::duration<double>(factored-uploaded).count(),
                std::chrono::duration<double>(downloaded-factored).count());
        return true;
    } catch(const std::bad_alloc &) { return false; }
      catch(const std::length_error &) { return false; }
}
}
bool LU_ENTRY(Rgetrf)(INTEGER const m,INTEGER const n,REAL *a,INTEGER const lda,INTEGER *ipiv,INTEGER &info) {
    return lu_call<mplapack_lapack_gpu::real_lu<LuOps>>("Rgetrf",m,n,a,lda,ipiv,info);
}
bool LU_ENTRY(Cgetrf)(INTEGER const m,INTEGER const n,COMPLEX *a,INTEGER const lda,INTEGER *ipiv,INTEGER &info) {
    return lu_call<mplapack_lapack_gpu::complex_lu<LuOps>>("Cgetrf",m,n,a,lda,ipiv,info);
}
#undef LU_ENTRY
