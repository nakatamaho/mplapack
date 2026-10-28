// Included by each dd/td/qd CUDA TU (or the host-emulation test).
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
#elif defined(MPLAPACK_BUILD_WITH_QD)
#include <mpblas_qd_cuda.h>
#endif
#pragma GCC visibility push(hidden)
#include "Rpotrf_gpu_common.cuh"
#pragma GCC visibility pop

namespace {
#if defined(MPLAPACK_BUILD_WITH_DD)
using PotrfOps=mplapack_gemm::dd_ops;
#define POTRF_ENTRY(R) R##_dd_cuda
static constexpr const char *potrf_prefix="MPLAPACK_DD";
static constexpr long potrf_min=4096;
static constexpr long potrf_complex_min=potrf_min;
#elif defined(MPLAPACK_BUILD_WITH_TD)
using PotrfOps=mplapack_gemm::td_ops;
#define POTRF_ENTRY(R) R##_td_cuda
static constexpr const char *potrf_prefix="MPLAPACK_TD";
static constexpr long potrf_min=2048; // Measured real/complex, both UPLOs.
static constexpr long potrf_complex_min=potrf_min;
#elif defined(MPLAPACK_BUILD_WITH_QD)
using PotrfOps=mplapack_gemm::qd_ops;
#define POTRF_ENTRY(R) R##_qd_cuda
static constexpr const char *potrf_prefix="MPLAPACK_QD";
static constexpr long potrf_min=2048; // Measured crossover after device residency.
static constexpr long potrf_complex_min=4096; // Upper loses at 2048; both UPLOs win at 4096.
#else
#error dd/td/qd only
#endif
using PotrfReal=mplapack_blas_dev::dreal<PotrfOps>;
using PotrfComplex=mplapack_blas_dev::dcomplex<PotrfReal>;
static PotrfReal potrf_pack(const REAL &a) {
    PotrfReal r;
    for(unsigned q=0;q<sizeof(r.v.x)/sizeof(double);++q) r.v.x[q]=a.x[q];
    return r;
}
static PotrfComplex potrf_pack(const COMPLEX &a) { return PotrfComplex(potrf_pack(a.real()),potrf_pack(a.imag())); }
static void potrf_unpack(const PotrfReal &r, REAL &a) {
    for(unsigned q=0;q<sizeof(r.v.x)/sizeof(double);++q) a.x[q]=r.v.x[q];
}
static void potrf_unpack(const PotrfComplex &r, COMPLEX &a) {
    REAL re,im;
    potrf_unpack(r.real(),re); potrf_unpack(r.imag(),im);
    a=COMPLEX(re,im);
}
static const char *potrf_env(const char *suffix) {
    char key[96];
    std::snprintf(key,sizeof(key),"%s%s",potrf_prefix,suffix);
    return std::getenv(key);
}
template<class Context,class Host> static bool potrf_call(const char *name,const char *uplo,INTEGER n,Host *a,INTEGER lda,INTEGER &info) {
    using namespace mplapack_lapack_gpu;
    using E=typename Context::E;
    const bool upper=uplo && (*uplo=='U' || *uplo=='u');
    if(!uplo || (!upper && *uplo!='L' && *uplo!='l') || n<=0 || lda<n || !a) return false;
    const char *enable=potrf_env("_LAPACK_CUDA");
    if(enable && std::strcmp(enable,"0")==0) return false;
    const char *threshold=potrf_env("_LAPACK_CUDA_MIN");
    long minimum=name[0]=='C'?potrf_complex_min:potrf_min;
    if(threshold) {
        char *end=nullptr;
        long v=std::strtol(threshold,&end,10);
        if(end!=threshold && !*end && v>=0) minimum=v;
    }
    if(n<minimum || (size_t)n>std::numeric_limits<size_t>::max()/(size_t)n/sizeof(E)) return false;
    try {
        const long nb=iMlaenv(1,name,uplo,n,-1,-1,-1);
        const long workspace=std::max(1L,std::min((long)n,nb));
        matrix<E> device;
        if(device.allocate(n,n)) return false;
        Context ctx(n,workspace);
        if(ctx.error) return false;
        std::vector<E> staging((size_t)n*n);
        using Clock=std::chrono::steady_clock;
        auto start=Clock::now();
        for(long j=0;j<n;++j) for(long i=0;i<n;++i)
            if(upper?i<=j:i>=j) staging[i+j*n]=potrf_pack(a[i+j*lda]);
        if(device.copy(staging.data(),n,false)) return false;
#if defined(__CUDACC__)
        // Nonblocking panel/update streams must wait for pageable H2D DMA.
        if(cudaStreamSynchronize(nullptr)!=cudaSuccess) return false;
#endif
        auto uploaded=Clock::now();
        long result=cholesky(ctx,upper,n,nb,device.data(),device.ld);
        if(ctx.error) return false;
        auto factored=Clock::now();
        // Transactional publication: CUDA errors cannot leave partial host output.
        if(device.copy(staging.data(),n,true)) return false;
        for(long j=0;j<n;++j) for(long i=0;i<n;++i)
            if(upper?i<=j:i>=j) potrf_unpack(staging[i+j*n],a[i+j*lda]);
        info=result;
        auto downloaded=Clock::now();
        const char *verbose=potrf_env("_CUDA_VERBOSE");
        if(verbose && std::strcmp(verbose,"0"))
            std::fprintf(stderr,"%s %s CUDA n=%ld upload=%.9g factor=%.9g download=%.9g\n",name,potrf_prefix,(long)n,
                         std::chrono::duration<double>(uploaded-start).count(),std::chrono::duration<double>(factored-uploaded).count(),
                         std::chrono::duration<double>(downloaded-factored).count());
        return true;
    } catch(const std::bad_alloc &) { return false; }
      catch(const std::length_error &) { return false; }
}
} // namespace
bool POTRF_ENTRY(Rpotrf)(const char *uplo,INTEGER const n,REAL *a,INTEGER const lda,INTEGER &info) {
    return potrf_call<mplapack_lapack_gpu::real_cholesky<PotrfOps>>("Rpotrf",uplo,n,a,lda,info);
}
bool POTRF_ENTRY(Cpotrf)(const char *uplo,INTEGER const n,COMPLEX *a,INTEGER const lda,INTEGER &info) {
    return potrf_call<mplapack_lapack_gpu::complex_cholesky<PotrfOps>>("Cpotrf",uplo,n,a,lda,info);
}
#undef POTRF_ENTRY
