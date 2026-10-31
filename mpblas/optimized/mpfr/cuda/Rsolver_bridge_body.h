// Included only in MPFR bridge TUs; reuse packing and factor engines locally.
namespace {
static bool solve_component_ok(mpfr_srcptr x,long p) {
    if(mpfr_get_prec(x)!=p || !mpfr_number_p(x)) return false;
    if(mpfr_zero_p(x)) return true;
    auto exponent=mpfr_get_exp(x);
    return exponent>=-RT_DEVICE_EMAX && exponent<=RT_DEVICE_EMAX;
}
#ifdef MPLAPACK_SOLVER_REAL_ONLY
using SolveHost=mpfr_class;
static bool solve_value_ok(const SolveHost &x,long p) {
    return solve_component_ok(x.get_mpfr_t(),p);
}
static void solve_publish(SolveHost &a,SolveHost &b) { mpfr_swap(a.get_mpfr_t(),b.get_mpfr_t()); }
#else
using SolveHost=mpc_class;
static bool solve_value_ok(const SolveHost &x,long p) {
    auto z=x.mpc_data();
    return solve_component_ok(mpc_realref(z),p) && solve_component_ok(mpc_imagref(z),p);
}
static void solve_publish(SolveHost &a,SolveHost &b) { mpc_swap(a.mpc_data(),b.mpc_data()); }
#endif
struct solver_host_delete {
    void operator()(SolveHost *p) const { delete[] p; }
    void operator()(long *p) const { delete[] p; }
};
// Inline backend accessors/constructors into this private bridge instead of
// adding weak C++ implementation symbols to the library's public ABI.
__attribute__((flatten)) static bool solver_call(const char *name,bool factorize,bool chol,const char *flag,long n,long rhs,
    SolveHost *a,long lda,mplapackint *ipiv,SolveHost *b,long ldb,mplapackint &info) {
    char f=flag?*flag:0; if(f>='a' && f<='z') f-=32;
    if(n<=0 || rhs<=0 || lda<n || ldb<n || !a || !b || (!chol && !ipiv) ||
       (chol?(f!='U' && f!='L'):(f!='N' && f!='T' && f!='C'))) return false;
    // Keep the measured factor gates for fused drivers; standalone solves
    // need the measured 2048-by-2048 RHS workload to amortize GPU residency.
    const long default_minimum=factorize?(chol?2048:(name[0]=='C'?1024:512)):2048;
    long minimum=env_long("MPLAPACK_MPFR_LAPACK_CUDA_MIN",default_minimum);
    if(minimum<0) minimum=default_minimum;
    long specific=env_long(factorize?"MPLAPACK_MPFR_LAPACK_DRIVER_CUDA_MIN":"MPLAPACK_MPFR_LAPACK_SOLVE_CUDA_MIN",minimum);
    if(specific>=0) minimum=specific;
    if(!env_long("MPLAPACK_MPFR_CUDA",1) || !env_long("MPLAPACK_MPFR_LAPACK_CUDA",1) || n<minimum) return false;
    if(!factorize) {
        long min_rhs=env_long("MPLAPACK_MPFR_LAPACK_SOLVE_CUDA_MIN_RHS",2048);
        if(min_rhs<0) min_rhs=2048;
        if(rhs<min_rhs) return false;
    }
    if(mpfr_get_default_rounding_mode()!=MPFR_RNDN || mpfr_get_emax()!=RT_DEVICE_EMAX || mpfr_get_emin()!=-RT_DEVICE_EMAX) return false;
    long prec=mpfrxx::default_precision_bits();
    if(prec<MPFR_PREC_MIN || prec>std::numeric_limits<int>::max()-1024L) return false;
#ifndef MPLAPACK_SOLVER_REAL_ONLY
    if(prec>4096) return false;
    if(mpfrxx::default_mpc_real_rounding_mode()!=MPFR_RNDN || mpfrxx::default_mpc_imag_rounding_mode()!=MPFR_RNDN ||
       mpfrxx::default_mpc_real_precision_bits()!=prec || mpfrxx::default_mpc_imag_precision_bits()!=prec) return false;
#endif
    size_t max=std::numeric_limits<size_t>::max();
    size_t element_bytes=size_t(2*((prec+63)/64*8+32));
    if(size_t(n)>max/size_t(n)/element_bytes || size_t(rhs)>max/size_t(n)/element_bytes-size_t(n) ||
       size_t(n)>max/size_t(lda)/sizeof(SolveHost) || size_t(rhs)>max/size_t(ldb)/sizeof(SolveHost)) return false;
    if(!factorize && !chol) for(long i=0;i<n;++i) if(ipiv[i]<1 || ipiv[i]>n) return false;
    try {
        std::lock_guard<std::mutex> lock(gpu_mutex());
        bool fixed=false;
        if(!env_long("MPLAPACK_MPFR_CUDA_FORCE_RUNTIME",0)) switch(prec) {
#define FIXED(PB) case PB: fixed=true; break;
            MPLAPACK_MPFR_CUDA_FIXED_PRECISIONS(FIXED)
#undef FIXED
        }
        if(!fixed && !env_long("MPLAPACK_MPFR_CUDA_RUNTIME",1)) return false;
#ifdef MPLAPACK_SOLVER_REAL_ONLY
        if(!factor_memory_fits(n,n+rhs,prec,fixed)) return false;
#else
        if(!complex_factor_memory_fits(n,n+rhs,prec,fixed)) return false;
#endif
        // Staging also guarantees transactional publication after range/runtime errors.
        auto conversion_start=Clock::now();
        std::unique_ptr<SolveHost[],solver_host_delete> stage(new SolveHost[size_t(n)*(n+rhs)]);
        for(long j=0;j<n;++j) for(long i=0;i<n;++i) {
            // Match Rfactor/Cfactor eligibility for the whole input matrix,
            // including an unused Cholesky triangle, but never publish it.
            if(!solve_value_ok(a[i+j*lda],prec)) return false;
            if(chol && (f=='U'?i>j:i<j)) continue;
            stage[i+j*n]=a[i+j*lda];
        }
        for(long j=0;j<rhs;++j) for(long i=0;i<n;++i) {
            if(!solve_value_ok(b[i+j*ldb],prec)) return false;
            stage[i+(j+n)*n]=b[i+j*ldb];
        }
        std::unique_ptr<long[],solver_host_delete> piv(new long[n]);
        if(!factorize && !chol) for(long i=0;i<n;++i) piv[i]=ipiv[i];
        factor_shape s{chol,chol && f=='U',n,n,iMlaenv_mpfr(1,name[0]=='C'?(chol?"Cpotrf":"Cgetrf"):(chol?"Rpotrf":"Rgetrf"),chol?flag:" ",n,chol?-1:n,-1,-1)};
        s.nrhs=rhs; s.factorize=factorize; s.trans=chol?'N':f;
        long result=0; factor_times t; bool success=false;
        double staging_seconds=std::chrono::duration<double>(Clock::now()-conversion_start).count();
        if(fixed) switch(prec) {
#ifdef MPLAPACK_SOLVER_REAL_ONLY
#define DISPATCH(PB) case PB: success=fixed_call<PB>(s,stage.get(),n,piv.get(),result,t); break;
#else
#define DISPATCH(PB) case PB: success=call<factor_fixed_tag<PB>>(s,prec,stage.get(),n,piv.get(),result,t); break;
#endif
            MPLAPACK_MPFR_CUDA_FIXED_PRECISIONS(DISPATCH)
#undef DISPATCH
        } else {
#ifdef MPLAPACK_SOLVER_REAL_ONLY
            success=runtime_call(s,prec,stage.get(),n,piv.get(),result,t);
#else
#define RT(NL) if(prec<=64*NL) success=call<factor_runtime_tag<NL>>(s,prec,stage.get(),n,piv.get(),result,t); else
            MPLAPACK_MPFR_CUDA_RT_TIERS(RT)
#undef RT
            return false;
#endif
        }
        if(!success) return false;
        t.upload+=staging_seconds;
        auto publication_start=Clock::now();
        if(factorize) for(long j=0;j<n;++j) for(long i=0;i<n;++i) if(!chol || (f=='U'?i<=j:i>=j)) solve_publish(a[i+j*lda],stage[i+j*n]);
        if(!result) for(long j=0;j<rhs;++j) for(long i=0;i<n;++i) solve_publish(b[i+j*ldb],stage[i+(j+n)*n]);
        if(factorize && !chol) for(long i=0;i<n;++i) ipiv[i]=piv[i];
        info=result;
        t.download+=std::chrono::duration<double>(Clock::now()-publication_start).count();
        if(env_long("MPLAPACK_MPFR_CUDA_VERBOSE",0)) std::fprintf(stderr,"%s MPLAPACK_MPFR CUDA p=%ld n=%ld nrhs=%ld upload=%.9g factor_solve=%.9g download=%.9g\n",name,prec,n,rhs,t.upload,t.factor,t.download);
        return true;
    } catch(const std::bad_alloc &) { return false; } catch(const std::length_error &) { return false; }
}
}
using INTEGER=mplapackint;
using REAL=mpfr_class;
using COMPLEX=mpc_class;
#define SOLVER_ENTRY(R) R##_mpfr_cuda
#define SOLVER_CALL solver_call
#include "../../../../mplapack/optimized/common/Rsolver_entries.h"
