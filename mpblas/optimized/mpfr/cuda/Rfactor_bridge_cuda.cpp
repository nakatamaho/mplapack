// Public MPFR LAPACK accelerator entries. CPU fallback belongs to dispatchers.
#include <mpblas_mpfr_cuda.h>
#include <mplapack_mpfr.h>
#include <chrono>
#include <cstdio>
#include <limits>
#include <memory>
#include <stdexcept>
#include "mpfr_pack_cuda.h"
#include "Rfactor_kernel_cuda.h"
namespace {
using namespace mplapack_mpfr_cuda;
using namespace mplapack_mpfr_cuda_host;
using Clock=std::chrono::steady_clock;
// A private deleter keeps this ownership template out of the shared ABI.
struct pivot_array_delete { void operator()(long *p) const { delete[] p; } };
static bool referenced(const factor_shape &s,long i,long j) { return j>=s.n || !s.chol || (s.upper?i<=j:i>=j); }
template<int PB> bool fixed_call(const factor_shape &s,mpfr_class *a,long lda,long *piv,long &info,factor_times &t) {
    std::vector<real_t<PB>> packed;
    auto start=Clock::now();
    if(!pack_matrix<PB>(a,lda,s.m,s.n+s.nrhs,packed,RT_DEVICE_EMAX)) return false;
    // pack() intentionally drops zero signs for the old GEMM convention.
    // The factorization engine, like fixed Rtrsm, retains signed zeros.
    for(long j=0;j<s.n+s.nrhs;++j) for(long i=0;i<s.m;++i)
        if(packed[i+j*s.m].is_zero()) packed[i+j*s.m].sign=mpfr_signbit(a[i+j*lda].get_mpfr_t())?-1:1;
    real_t<PB> safe;
    mpfr_class sfmin=Rlamch_mpfr("S");
    if(!pack<PB>(sfmin.get_mpfr_t(),safe,RT_DEVICE_EMAX)) return false;
    auto packed_at=Clock::now();
    if(factor_fixed<PB>(s,packed.data(),safe,piv,info,t)) return false;
    for(const auto &v:packed) if(!fits<PB>(v)) return false;
    auto before=Clock::now();
    for(long j=0;j<s.n+s.nrhs;++j) for(long i=0;i<s.m;++i) if(referenced(s,i,j)) {
        const auto &v=packed[i+j*s.m]; mpfr_ptr x=a[i+j*lda].get_mpfr_t();
        if(v.is_zero()) mpfr_set_zero(x,v.sign); else unpack<PB>(v,x);
    }
    t.upload+=std::chrono::duration<double>(packed_at-start).count();
    t.download+=std::chrono::duration<double>(Clock::now()-before).count();
    return true;
}
static bool runtime_call(const factor_shape &s,long prec,mpfr_class *a,long lda,long *piv,long &info,factor_times &t) {
    int nl=(prec+63)/64; rt_storage packed,safe;
    auto start=Clock::now();
    if(!rt_pack_matrix(a,lda,s.m,s.n+s.nrhs,prec,nl,packed,RT_DEVICE_EMAX)) return false;
    safe.resize(1,nl); mpfr_class sfmin=Rlamch_mpfr("S");
    if(!rt_pack_value(sfmin.get_mpfr_t(),prec,nl,safe,0,RT_DEVICE_EMAX)) return false;
    auto packed_at=Clock::now();
    if(factor_runtime(s,prec,packed.pack,safe.pack,piv,info,t)) return false;
    for(size_t i=0;i<(size_t)s.m*(s.n+s.nrhs);++i) if(!rt_fits(packed,i)) return false;
    auto before=Clock::now();
    for(long j=0;j<s.n+s.nrhs;++j) for(long i=0;i<s.m;++i) if(referenced(s,i,j))
        rt_unpack(packed,i+j*s.m,prec,nl,a[i+j*lda].get_mpfr_t());
    t.upload+=std::chrono::duration<double>(packed_at-start).count();
    t.download+=std::chrono::duration<double>(Clock::now()-before).count();
    return true;
}
static bool factor_call(bool chol,const char *uplo,long m,long n,mpfr_class *a,long lda,mplapackint *ipiv,mplapackint &info) {
    if(m<=0 || n<=0 || lda<m || !a || (!chol && !ipiv)) return false;
    bool upper=chol && uplo && (*uplo=='U' || *uplo=='u');
    if(chol && !upper && (!uplo || (*uplo!='L' && *uplo!='l'))) return false;
    // Three interleaved rounds, 512-bit A100 measurements: resident LU wins
    // at 512; Cholesky loses at 512/1024 and wins at 2048/4096.
    const long default_minimum=chol?2048:512;
    long minimum=env_long("MPLAPACK_MPFR_LAPACK_CUDA_MIN",default_minimum);
    if(minimum<0) minimum=default_minimum;
    if(!env_long("MPLAPACK_MPFR_CUDA",1) || !env_long("MPLAPACK_MPFR_LAPACK_CUDA",1) ||
       std::min(m,n)<minimum) return false;
    if(mpfr_get_default_rounding_mode()!=MPFR_RNDN ||
       mpfr_get_emax()!=RT_DEVICE_EMAX || mpfr_get_emin()!=-RT_DEVICE_EMAX) return false;
    long prec=mpfrxx::default_precision_bits();
    if(prec<MPFR_PREC_MIN || prec>std::numeric_limits<int>::max()-1024L) return false;
    // Reject both device-size and host-address overflow before reading A.
    if((size_t)n>std::numeric_limits<size_t>::max()/(size_t)m/((prec+63)/64*8+16) ||
       (size_t)n>std::numeric_limits<size_t>::max()/(size_t)lda/sizeof(mpfr_class)) return false;
    try {
        std::lock_guard<std::mutex> lock(gpu_mutex());
        bool fixed=false;
        if(!env_long("MPLAPACK_MPFR_CUDA_FORCE_RUNTIME",0)) switch(prec) {
#define FIXED_CASE(PB) case PB: fixed=true; break;
            MPLAPACK_MPFR_CUDA_FIXED_PRECISIONS(FIXED_CASE)
#undef FIXED_CASE
        }
        if(!fixed && !env_long("MPLAPACK_MPFR_CUDA_RUNTIME",1)) return false;
        if(!factor_memory_fits(m,n,prec,fixed)) return false;
        factor_shape s{chol,upper,m,n,iMlaenv_mpfr(1,chol?"Rpotrf":"Rgetrf",chol?uplo:" ",m,chol?-1:n,-1,-1)};
        long npiv=std::min(m,n);
        std::unique_ptr<long[],pivot_array_delete> piv(new long[npiv]);
        long result=0; factor_times t;
        bool success=false;
        if(!env_long("MPLAPACK_MPFR_CUDA_FORCE_RUNTIME",0)) {
            switch(prec) {
#define DISPATCH(PB) case PB: success=fixed_call<PB>(s,a,lda,piv.get(),result,t); break;
                MPLAPACK_MPFR_CUDA_FIXED_PRECISIONS(DISPATCH)
#undef DISPATCH
                default:
                    if(env_long("MPLAPACK_MPFR_CUDA_RUNTIME",1)) success=runtime_call(s,prec,a,lda,piv.get(),result,t);
            }
        } else if(env_long("MPLAPACK_MPFR_CUDA_RUNTIME",1)) success=runtime_call(s,prec,a,lda,piv.get(),result,t);
        if(!success) return false;
        if(!chol) for(long i=0;i<npiv;++i) ipiv[i]=piv[i];
        info=result;
        if(env_long("MPLAPACK_MPFR_CUDA_VERBOSE",0))
            std::fprintf(stderr,"%s MPLAPACK_MPFR CUDA p=%ld m=%ld n=%ld upload=%.9g factor=%.9g download=%.9g\n",
                         chol?"Rpotrf":"Rgetrf",prec,m,n,t.upload,t.factor,t.download);
        return true;
    } catch(const std::bad_alloc &) { return false; }
      catch(const std::length_error &) { return false; }
}
}
bool Rgetrf_mpfr_cuda(mplapackint const m,mplapackint const n,mpfr_class *a,mplapackint const lda,mplapackint *ipiv,mplapackint &info) {
    return factor_call(false,nullptr,m,n,a,lda,ipiv,info);
}
bool Rpotrf_mpfr_cuda(const char *uplo,mplapackint const n,mpfr_class *a,mplapackint const lda,mplapackint &info) {
    return factor_call(true,uplo,n,n,a,lda,nullptr,info);
}
#define MPLAPACK_SOLVER_REAL_ONLY
#include "Rsolver_bridge_body.h"
