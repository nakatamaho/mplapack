// Resident solver entries; factor arithmetic is reused, not reimplemented.
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
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <vector>
#if defined(MPLAPACK_BUILD_WITH_DD)
#include <mpblas_dd_cuda.h>
#define SOLVER_ENTRY(R) R##_dd_cuda
#elif defined(MPLAPACK_BUILD_WITH_TD)
#include <mpblas_td_cuda.h>
#define SOLVER_ENTRY(R) R##_td_cuda
#else
#include <mpblas_qd_cuda.h>
#define SOLVER_ENTRY(R) R##_qd_cuda
#endif
#pragma GCC visibility push(hidden)
#include "Rgetrf_gpu_common.cuh"
#pragma GCC visibility pop
namespace {
using namespace mplapack_lapack_gpu;
#if defined(MPLAPACK_BUILD_WITH_DD)
using SolveOps=mplapack_gemm::dd_ops;
constexpr const char *solve_prefix="MPLAPACK_DD";
constexpr long solve_lu_min=2048,solve_chol_min=4096,solve_cchol_min=4096;
#elif defined(MPLAPACK_BUILD_WITH_TD)
using SolveOps=mplapack_gemm::td_ops;
constexpr const char *solve_prefix="MPLAPACK_TD";
constexpr long solve_lu_min=1024,solve_chol_min=2048,solve_cchol_min=2048;
#else
using SolveOps=mplapack_gemm::qd_ops;
constexpr const char *solve_prefix="MPLAPACK_QD";
constexpr long solve_lu_min=1024,solve_chol_min=2048,solve_cchol_min=4096;
#endif
using SR=dreal<SolveOps>; using SC=dcomplex<SR>;
static SR pack_solve(const REAL &a) {
    SR v; for(unsigned q=0;q<sizeof(v.v.x)/sizeof(double);++q) v.v.x[q]=a.x[q]; return v;
}
static SC pack_solve(const COMPLEX &a) { return SC(pack_solve(a.real()),pack_solve(a.imag())); }
static void unpack_solve(const SR &v,REAL &a) {
    for(unsigned q=0;q<sizeof(v.v.x)/sizeof(double);++q) a.x[q]=v.v.x[q];
}
static void unpack_solve(const SC &v,COMPLEX &a) {
    REAL r,i; unpack_solve(v.real(),r); unpack_solve(v.imag(),i); a=COMPLEX(r,i);
}
static long solve_setting(const char *suffix,long fallback) {
    char key[96]; std::snprintf(key,sizeof(key),"%s%s",solve_prefix,suffix);
    if(const char *v=std::getenv(key)) {
        char *end=nullptr; long x=std::strtol(v,&end,10);
        if(end!=v && !*end && x>=0) return x;
    }
    return fallback;
}
static void triangular(chol_schedule &ctx,bool upper,char trans,bool unit,long n,long rhs,SR *a,SR *b) {
    using F=SolveOps::F;
    buffer<F> one; F h=SolveOps::from_d(1);
    if(!ctx.error) ctx.error=one.allocate(1);
    if(!ctx.error) ctx.error=transfer(one.data(),&h,1,false);
#if defined(__CUDACC__)
    if(!ctx.error) ctx.error=cudaStreamSynchronize(nullptr);
#endif
    if(ctx.error) return;
    mplapack_blas::gpu_lane_args<SolveOps> g{};
    g.s.op=mplapack_blas::L3_TRSM; g.s.left=1; g.s.upper=upper;
    g.s.trans=trans!='N'; g.s.nounit=!unit; g.s.alpha_one=1; g.s.m=n; g.s.n=rhs;
    g.a=reinterpret_cast<F*>(a); g.c=reinterpret_cast<F*>(b); g.al=one.data(); g.lda=n; g.ldc=n;
#if defined(__CUDACC__)
    mplapack_blas::k_blas_lanes<SolveOps><<<(unsigned)((rhs+63)/64),64,0,ctx.stream()>>>(g,rhs);
    ctx.error=cudaGetLastError();
#else
    ctx.error=mplapack_blas::blas_device<SolveOps>(g);
#endif
    ctx.finish(); // scalar lifetime, no matrix transfer
}
static void triangular(chol_schedule &ctx,bool upper,char trans,bool unit,long n,long rhs,SC *a,SC *b) {
    job_Ctrsm<SolveOps> j{};
    j.side[0]='L'; j.uplo[0]=upper?'U':'L'; j.ta[0]=trans; j.diag[0]=unit?'U':'N';
    j.m=n; j.n=rhs; j.lda=n; j.ldb=n; j.alpha=SC(1,0); j.A=a; j.B=b;
    ctx.launch(j,j.units());
}
struct solver_piv_delete { void operator()(long *p) const { delete[] p; } };
template<class Host> static bool solver_call(const char *name,bool factorize,bool chol,const char *flag,
    INTEGER n,INTEGER rhs,Host *a,INTEGER lda,INTEGER *ipiv,Host *b,INTEGER ldb,INTEGER &info) {
    using E=decltype(pack_solve(*a));
    char f=flag?*flag:0; if(f>='a' && f<='z') f-=32;
    if(n<=0 || rhs<=0 || lda<n || ldb<n || !a || !b || (!chol && !ipiv) ||
       (chol?(f!='U' && f!='L'):(f!='N' && f!='T' && f!='C'))) return false;
    // The measured BF standalone solves did not beat the existing path:
    // opt in with a common/solve-specific minimum instead of guessing a
    // profitable crossover beyond the measured range. Drivers retain the
    // resident factorization gates (see Rsolver_cuda.md and its timing CSV).
    long minimum=factorize?(chol?(name[0]=='C'?solve_cchol_min:solve_chol_min):solve_lu_min):std::numeric_limits<long>::max();
    minimum=solve_setting("_LAPACK_CUDA_MIN",minimum);
    minimum=solve_setting(factorize?"_LAPACK_DRIVER_CUDA_MIN":"_LAPACK_SOLVE_CUDA_MIN",minimum);
    if(!solve_setting("_LAPACK_CUDA",1) || n<minimum) return false;
    if(!factorize && rhs<solve_setting("_LAPACK_SOLVE_CUDA_MIN_RHS",0)) return false;
    size_t max=std::numeric_limits<size_t>::max();
    if(size_t(n)>max/size_t(n)/sizeof(E) || size_t(rhs)>max/size_t(n)/sizeof(E)-size_t(n) ||
       size_t(n)>max/size_t(lda)/sizeof(Host) || size_t(rhs)>max/size_t(ldb)/sizeof(Host)) return false;
    if(!factorize && !chol) for(long i=0;i<n;++i) if(ipiv[i]<1 || ipiv[i]>n) return false;
    try {
        bool upper=chol && f=='U';
        long nb=iMlaenv(1,name[0]=='C'?(chol?"Cpotrf":"Cgetrf"):(chol?"Rpotrf":"Rgetrf"),chol?flag:" ",n,chol?-1:n,-1,-1);
        long work=std::max(1L,std::min(long(n),nb));
#if defined(__CUDACC__)
        size_t free=0,total=0;
        size_t count=size_t(n)*(size_t(n)+size_t(rhs));
        size_t overhead=2*size_t(n)*size_t(work)*sizeof(E)+size_t(n)*sizeof(long)+1048576;
        if(cudaMemGetInfo(&free,&total)!=cudaSuccess || free<overhead || count>(free-overhead)/sizeof(E)) return false;
#endif
        matrix<E> da,db; buffer<long> dp;
        if(da.allocate(n,n) || db.allocate(n,rhs) || dp.allocate(n)) return false;
        std::vector<E> ha(size_t(n)*n),hb(size_t(n)*rhs);
        std::unique_ptr<long[],solver_piv_delete> piv(new long[n]);
        auto start=std::chrono::steady_clock::now();
        for(long j=0;j<n;++j) for(long i=0;i<n;++i) if(!chol || (upper?i<=j:i>=j)) ha[i+j*n]=pack_solve(a[i+j*lda]);
        for(long j=0;j<rhs;++j) for(long i=0;i<n;++i) hb[i+j*n]=pack_solve(b[i+j*ldb]);
        if(da.copy(ha.data(),n,false) || db.copy(hb.data(),n,false)) return false;
        if(!factorize && !chol) {
            for(long i=0;i<n;++i) piv[i]=ipiv[i];
            if(transfer(dp.data(),piv.get(),n,false)) return false;
        }
#if defined(__CUDACC__)
        if(cudaStreamSynchronize(nullptr)!=cudaSuccess) return false;
#endif
        auto uploaded=std::chrono::steady_clock::now();
        long result=0;
        if(factorize) {
            if(chol) {
                using C=typename std::conditional<std::is_same<E,SR>::value,real_cholesky<SolveOps>,complex_cholesky<SolveOps>>::type;
                C ctx(n,work);
#if defined(__CUDACC__)
                // The context's constants are uploaded on the default stream,
                // after the matrix fence; its nonblocking streams must wait.
                if(!ctx.error) ctx.error=cudaStreamSynchronize(nullptr);
#endif
                result=cholesky(ctx,upper,n,nb,da.data(),n); ctx.finish(); if(ctx.error) return false;
            } else {
                using L=typename std::conditional<std::is_same<E,SR>::value,real_lu<SolveOps>,complex_lu<SolveOps>>::type;
                lu_context<L,SR> ctx(n,work,pack_solve(Rlamch("S")));
#if defined(__CUDACC__)
                if(!ctx.error) ctx.error=cudaStreamSynchronize(nullptr);
#endif
                result=lu(ctx,n,n,nb,da.data(),n,dp.data()); if(ctx.error) return false;
            }
        }
        if(!result) {
            chol_schedule ctx;
            if(!chol && f=='N') ctx.launch(swaps<E>{db.data(),n,1,n,1,dp.data()},rhs);
            if(chol) {
                triangular(ctx,upper,upper?'C':'N',false,n,rhs,da.data(),db.data());
                triangular(ctx,upper,upper?'N':'C',false,n,rhs,da.data(),db.data());
            } else if(f=='N') {
                triangular(ctx,false,'N',true,n,rhs,da.data(),db.data());
                triangular(ctx,true,'N',false,n,rhs,da.data(),db.data());
            } else {
                triangular(ctx,true,f,false,n,rhs,da.data(),db.data());
                triangular(ctx,false,f,true,n,rhs,da.data(),db.data());
                ctx.launch(swaps<E>{db.data(),n,1,n,-1,dp.data()},rhs);
            }
            ctx.finish(); if(ctx.error) return false;
        }
        auto solved=std::chrono::steady_clock::now();
        if((factorize && da.copy(ha.data(),n,true)) || (!result && db.copy(hb.data(),n,true)) ||
           (factorize && !chol && transfer(piv.get(),dp.data(),n,true))) return false;
        // No publication before every operation, including downloads, succeeds.
        if(factorize) for(long j=0;j<n;++j) for(long i=0;i<n;++i) if(!chol || (upper?i<=j:i>=j)) unpack_solve(ha[i+j*n],a[i+j*lda]);
        if(!result) for(long j=0;j<rhs;++j) for(long i=0;i<n;++i) unpack_solve(hb[i+j*n],b[i+j*ldb]);
        if(factorize && !chol) for(long i=0;i<n;++i) ipiv[i]=piv[i];
        info=result;
        if(solve_setting("_CUDA_VERBOSE",0)) std::fprintf(stderr,"%s %s CUDA n=%ld nrhs=%ld upload=%.9g factor_solve=%.9g download=%.9g\n",name,solve_prefix,long(n),long(rhs),
            std::chrono::duration<double>(uploaded-start).count(),std::chrono::duration<double>(solved-uploaded).count(),
            std::chrono::duration<double>(std::chrono::steady_clock::now()-solved).count());
        return true;
    } catch(const std::bad_alloc &) { return false; } catch(const std::length_error &) { return false; }
}
}
#define SOLVER_CALL solver_call
#include "Rsolver_entries.h"
