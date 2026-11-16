// BF inverse entry points. Opt-in until CPU-library crossovers are measured.
#if defined(__CUDACC__)
#define MPLAPACK_MINMAX_MPLAPACKINT_DEFINED
#endif
#include <mpblas.h>
#include <mplapack.h>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <type_traits>
#include <stdexcept>
#include <vector>
#if defined(MPLAPACK_BUILD_WITH_DD)
#include <mpblas_dd_cuda.h>
#define INVERSE_ENTRY(R) R##_dd_cuda
#define INVERSE_PREFIX "MPLAPACK_DD"
#define INVERSE_OPS mplapack_gemm::dd_ops
#elif defined(MPLAPACK_BUILD_WITH_TD)
#include <mpblas_td_cuda.h>
#define INVERSE_ENTRY(R) R##_td_cuda
#define INVERSE_PREFIX "MPLAPACK_TD"
#define INVERSE_OPS mplapack_gemm::td_ops
#else
#include <mpblas_qd_cuda.h>
#define INVERSE_ENTRY(R) R##_qd_cuda
#define INVERSE_PREFIX "MPLAPACK_QD"
#define INVERSE_OPS mplapack_gemm::qd_ops
#endif
#pragma GCC visibility push(hidden)
#include "Rinverse_gpu_common.cuh"
#pragma GCC visibility pop
namespace {
using InverseOps=INVERSE_OPS;
using IR=mplapack_blas_dev::dreal<InverseOps>;
using IC=mplapack_blas_dev::dcomplex<IR>;
static IR inverse_pack(const REAL &x) { IR v;for(unsigned q=0;q<sizeof(v.v.x)/sizeof(double);++q)v.v.x[q]=x.x[q];return v; }
static IC inverse_pack(const COMPLEX &x) { return IC(inverse_pack(x.real()),inverse_pack(x.imag())); }
static void inverse_unpack(const IR &v,REAL &x) { for(unsigned q=0;q<sizeof(v.v.x)/sizeof(double);++q)x.x[q]=v.v.x[q]; }
static void inverse_unpack(const IC &v,COMPLEX &x) { REAL r,i;inverse_unpack(v.real(),r);inverse_unpack(v.imag(),i);x=COMPLEX(r,i); }
static long inverse_setting(const char *suffix,long fallback) {
    char key[96];std::snprintf(key,sizeof(key),"%s%s",INVERSE_PREFIX,suffix);
    if(const char *value=std::getenv(key)) {char *end=nullptr;long n=std::strtol(value,&end,10);if(end!=value && !*end && n>=0)return n;}
    return fallback;
}
static char inverse_flag(const char *s) {char c=s?*s:0;return c>='a' && c<='z'?c-32:c;}
// kind: 0 TRTRI, 1 LAUUM, 2 POTRI, 3 GETRI. No CPU arithmetic/fallback here.
template<class H> static bool inverse_call(const char *name,int kind,const char *uplo,const char *diag,
    INTEGER n,H *a,INTEGER lda,INTEGER *ipiv,H *work,INTEGER lwork,INTEGER &info) {
    using E=decltype(inverse_pack(*a));
    constexpr bool complex=std::is_same<H,COMPLEX>::value;
    char u=inverse_flag(uplo),d=inverse_flag(diag);
    if(n<=0 || lda<n || !a || (kind<3 && u!='U' && u!='L') || (kind==0 && d!='U' && d!='N') ||
       (kind==3 && (!ipiv || !work || lwork<n))) return false;
    long minimum=inverse_setting("_LAPACK_CUDA_MIN",std::numeric_limits<long>::max());
    minimum=inverse_setting("_LAPACK_INVERSE_CUDA_MIN",minimum);
    if(!inverse_setting("_LAPACK_CUDA",1) || n<minimum) return false;
    size_t max=std::numeric_limits<size_t>::max();
    if(size_t(n)>max/size_t(n)/sizeof(E) || size_t(n)>max/size_t(lda)/sizeof(H))return false;
    if(kind==3)for(long i=0;i<n;++i)if(ipiv[i]<1 || ipiv[i]>n)return false;
    try {
        const bool upper=kind==3 || u=='U',unit=kind==0 && d=='U';
        const char *tri_name=complex?"Ctrtri":"Rtrtri",*lau_name=complex?"Clauum":"Rlauum";
        const long nbtri=iMlaenv(1,tri_name,kind==3?"UN":CHAR2(uplo,kind==0?diag:"N"),n,-1,-1,-1);
        const long nblau=iMlaenv(1,lau_name,uplo,n,-1,-1,-1);
        long nb=kind==3?iMlaenv(1,name," ",n,-1,-1,-1):1,nbmin=2,iws=n,optimal=1;
        if(nb>0 && n>std::numeric_limits<long>::max()/nb)return false;
        optimal=std::max(1L,long(n)*nb);
        if(kind==3 && nb>1 && nb<n) {
            iws=long(n)*nb;
            if(lwork<iws) {nb=lwork/n;nbmin=std::max(2L,long(iMlaenv(2,name," ",n,-1,-1,-1)));}
        }
        const bool blocked=nb>=nbmin && nb<n;
        const long columns=kind==3 && blocked?nb:1;
        if(size_t(columns)>max/size_t(n)/sizeof(E))return false;
        const size_t wcount=kind==3?size_t(n)*columns:0;
#if defined(__CUDACC__)
        size_t free=0,total=0;const size_t overhead=size_t(n)*sizeof(long)+1048576;
        if(cudaMemGetInfo(&free,&total)!=cudaSuccess || free<overhead || wcount>(free-overhead)/sizeof(E) ||
           size_t(n)*size_t(n)>(free-overhead)/sizeof(E)-wcount)return false;
#endif
        mplapack_lapack_gpu::matrix<E> matrix;
        mplapack_lapack_gpu::buffer<E> workspace;
        mplapack_lapack_gpu::buffer<long> pivots;
        if(matrix.allocate(n,n) || (kind==3 && (workspace.allocate(wcount) || pivots.allocate(n))))return false;
        mplapack_inverse_gpu::context<InverseOps,complex> ctx;if(ctx.error)return false;
        std::vector<E> stage(size_t(n)*n),wstage(wcount);
        std::vector<long> p(kind==3?size_t(n):0);
        using Clock=std::chrono::steady_clock;auto start=Clock::now();
        for(long j=0;j<n;++j)for(long i=0;i<n;++i)if((kind==3 || (upper?i<=j:i>=j)) && !(unit && i==j))stage[i+j*n]=inverse_pack(a[i+j*lda]);
        if(matrix.copy(stage.data(),n,false))return false;
        if(kind==3) {for(long i=0;i<n;++i)p[i]=ipiv[i];if(mplapack_lapack_gpu::transfer(pivots.data(),p.data(),n,false))return false;}
#if defined(__CUDACC__)
        // Nonblocking kernels must wait for all pageable default-stream DMA.
        if(cudaStreamSynchronize(nullptr)!=cudaSuccess)return false;
#endif
        auto uploaded=Clock::now();long result=0;
        if(kind==0 || kind==2)result=ctx.tri(upper,unit,n,nbtri,matrix.data(),n);
        if(!result && (kind==1 || kind==2))ctx.lauum(upper,n,nblau,matrix.data(),n);
        if(kind==3)result=ctx.getri(n,nbtri,nb,nbmin,matrix.data(),workspace.data(),pivots.data(),n);
        ctx.finish();if(ctx.error)return false;auto computed=Clock::now();
        if(matrix.copy(stage.data(),n,true) || (kind==3 && !result && mplapack_lapack_gpu::transfer(wstage.data(),workspace.data(),wcount,true)))return false;
        // All transfers have succeeded: publish only reference-visible storage.
        for(long j=0;j<n;++j)for(long i=0;i<n;++i)if((kind==3 || (upper?i<=j:i>=j)) && !(unit && i==j))inverse_unpack(stage[i+j*n],a[i+j*lda]);
        if(kind==3) {
            // Workspace is output-only scratch; do not read uninitialized host
            // work or overwrite its untouched cells. Every published slot was
            // written by the reference panel-copy loops (last block starts at 0).
            if(!result)for(long j=0;j<columns;++j)for(long i=j+1;i<n;++i)inverse_unpack(wstage[i+j*n],work[i+j*n]);
            work[0]=H(double(result?optimal:iws));
        }
        info=result;
        if(inverse_setting("_CUDA_VERBOSE",0))std::fprintf(stderr,"%s %s CUDA n=%ld upload=%.9g inverse=%.9g download=%.9g\n",name,INVERSE_PREFIX,long(n),
            std::chrono::duration<double>(uploaded-start).count(),std::chrono::duration<double>(computed-uploaded).count(),
            std::chrono::duration<double>(Clock::now()-computed).count());
        return true;
    } catch(const std::bad_alloc &) {return false;} catch(const std::length_error &) {return false;}
}
}
#define INVERSE_ENTRIES(H,P) \
bool INVERSE_ENTRY(P##trtri)(const char *u,const char *d,INTEGER const n,H *a,INTEGER const lda,INTEGER &info) {return inverse_call(#P "trtri",0,u,d,n,a,lda,nullptr,(H*)nullptr,0,info);} \
bool INVERSE_ENTRY(P##lauum)(const char *u,INTEGER const n,H *a,INTEGER const lda,INTEGER &info) {return inverse_call(#P "lauum",1,u,"N",n,a,lda,nullptr,(H*)nullptr,0,info);} \
bool INVERSE_ENTRY(P##potri)(const char *u,INTEGER const n,H *a,INTEGER const lda,INTEGER &info) {return inverse_call(#P "potri",2,u,"N",n,a,lda,nullptr,(H*)nullptr,0,info);} \
bool INVERSE_ENTRY(P##getri)(INTEGER const n,H *a,INTEGER const lda,INTEGER *p,H *work,INTEGER const lwork,INTEGER &info) {return inverse_call(#P "getri",3,"U","N",n,a,lda,p,work,lwork,info);}
INVERSE_ENTRIES(REAL,R)
INVERSE_ENTRIES(COMPLEX,C)
#undef INVERSE_ENTRIES
#undef INVERSE_ENTRY
#undef INVERSE_PREFIX
#undef INVERSE_OPS
