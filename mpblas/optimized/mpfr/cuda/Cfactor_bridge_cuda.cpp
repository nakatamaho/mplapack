// Public complex MPFR resident entries. Fallback exists only in dispatchers.
#include <mpblas_mpfr_cuda.h>
#include <mplapack_mpfr.h>
#include <chrono>
#include <cstdio>
#include <limits>
#include <memory>
#include <stdexcept>
#include "mpfr_pack_cuda.h"
#include "Cfactor_kernel_cuda.h"
namespace {
using namespace mplapack_mpfr_cuda;
using namespace mplapack_mpfr_cuda_host;
using namespace mplapack_blas_dev;
using Clock=std::chrono::steady_clock;
template<class Tag> struct complex_pack;
template<int PB> struct complex_pack<factor_fixed_tag<PB>> {
    using R=factor_real<PB>;
    static bool get(mpfr_srcptr x,R &v,long) {
        if(!pack<PB>(x,v.v,RT_DEVICE_EMAX)) return false;
        if(v.v.is_zero()) v.v.sign=mpfr_signbit(x)?-1:1;
        return true;
    }
    static bool fits(const R &v) { return v.valid && mplapack_mpfr_cuda_host::fits<PB>(v.v); }
    static void put(const R &v,mpfr_ptr x,long) {
        if(v.v.is_zero()) mpfr_set_zero(x,v.v.sign); else unpack<PB>(v.v,x);
    }
};
template<int NL> struct complex_pack<factor_runtime_tag<NL>> {
    using R=rreal<NL>;
    static bool get(mpfr_srcptr x,R &v,long p) {
        if(mpfr_get_prec(x)!=p || !mpfr_number_p(x)) return false;
        v.kind=mpfr_custom_get_kind(x); v.exp=0;
        for(int i=0;i<NL;++i) v.m[i]=0;
        if(mpfr_zero_p(x)) return true;
        v.exp=mpfr_get_exp(x);
        if(v.exp < -RT_DEVICE_EMAX || v.exp > RT_DEVICE_EMAX) return false;
        auto limbs=(const mp_limb_t*)mpfr_custom_get_significand(x);
        for(int i=0;i<(p+63)/64;++i) v.m[i]=limbs[i];
        return true;
    }
    static bool fits(const R &v) {
        int k=v.kind<0?-v.kind:v.kind;
        return k==MPFR_ZERO_KIND || (k==MPFR_REGULAR_KIND && v.exp>=mpfr_get_emin() && v.exp<=mpfr_get_emax());
    }
    static void put(const R &v,mpfr_ptr x,long p) {
        mpfr_t t;
        mpfr_custom_init_set(t,v.kind,v.exp,p,const_cast<mp_limb_t*>((const mp_limb_t*)v.m));
        mpfr_set(x,t,MPFR_RNDN);
    }
};
template<class Tag> bool call(const factor_shape &s,long prec,mpc_class *a,long lda,long *piv,long &info,factor_times &times) {
    using Pack=complex_pack<Tag>; using V=typename cgpu_types<Tag>::DC;
    auto start=Clock::now(); std::vector<V> packed(size_t(s.m)*(s.n+s.nrhs));
    for(long j=0;j<s.n+s.nrhs;++j) for(long i=0;i<s.m;++i) {
        auto z=a[i+j*lda].mpc_data(); auto &v=packed[i+j*s.m];
        if(!Pack::get(mpc_realref(z),v.re_,prec) || !Pack::get(mpc_imagref(z),v.im_,prec)) return false;
    }
    V safe; mpfr_class sfmin=Rlamch_mpfr("S");
    if(!Pack::get(sfmin.get_mpfr_t(),safe.re_,prec)) return false;
    auto packed_at=Clock::now();
    if(factor_complex<Tag>(s,prec,packed.data(),safe,piv,info,times)) return false;
    for(const auto &v:packed) if(!Pack::fits(v.re_) || !Pack::fits(v.im_)) return false;
    auto before=Clock::now();
    for(long j=0;j<s.n+s.nrhs;++j) for(long i=0;i<s.m;++i) if(j>=s.n || !s.chol || (s.upper?i<=j:i>=j)) {
        auto z=a[i+j*lda].mpc_data(); const auto &v=packed[i+j*s.m];
        Pack::put(v.re_,mpc_realref(z),prec); Pack::put(v.im_,mpc_imagref(z),prec);
    }
    times.upload+=std::chrono::duration<double>(packed_at-start).count();
    times.download+=std::chrono::duration<double>(Clock::now()-before).count();
    return true;
}
static bool factor_call(bool chol,const char *uplo,long m,long n,mpc_class *a,long lda,mplapackint *ipiv,mplapackint &info) {
    if(m<=0 || n<=0 || lda<m || !a || (!chol && !ipiv)) return false;
    bool upper=chol && uplo && (*uplo=='U' || *uplo=='u');
    if(chol && !upper && (!uplo || (*uplo!='L' && *uplo!='l'))) return false;
    // Three interleaved 512-bit rounds: LU loses to old CUDA at 512,
    // wins from 1024; Cholesky beats both old CUDA and CPU at 2048.
    long default_minimum=chol?2048:1024;
    long minimum=env_long("MPLAPACK_MPFR_LAPACK_CUDA_MIN",default_minimum);
    if(minimum<0) minimum=default_minimum;
    if(!env_long("MPLAPACK_MPFR_CUDA",1) || !env_long("MPLAPACK_MPFR_LAPACK_CUDA",1) || std::min(m,n)<minimum) return false;
    if(mpfr_get_default_rounding_mode()!=MPFR_RNDN || mpfrxx::default_mpc_real_rounding_mode()!=MPFR_RNDN ||
       mpfrxx::default_mpc_imag_rounding_mode()!=MPFR_RNDN || mpfr_get_emax()!=RT_DEVICE_EMAX || mpfr_get_emin()!=-RT_DEVICE_EMAX) return false;
    long prec=mpfrxx::default_precision_bits();
    if(prec<MPFR_PREC_MIN || prec>4096 || mpfrxx::default_mpc_real_precision_bits()!=prec || mpfrxx::default_mpc_imag_precision_bits()!=prec) return false;
    if(size_t(n)>std::numeric_limits<size_t>::max()/size_t(m)/(2*((prec+63)/64*8+16)) ||
       size_t(n)>std::numeric_limits<size_t>::max()/size_t(lda)/sizeof(mpc_class)) return false;
    try {
        std::lock_guard<std::mutex> lock(gpu_mutex());
        bool fixed=false;
        if(!env_long("MPLAPACK_MPFR_CUDA_FORCE_RUNTIME",0)) switch(prec) {
#define FIXED(PB) case PB: fixed=true; break;
            MPLAPACK_MPFR_CUDA_FIXED_PRECISIONS(FIXED)
#undef FIXED
        }
        if(!fixed && !env_long("MPLAPACK_MPFR_CUDA_RUNTIME",1)) return false;
        if(!complex_factor_memory_fits(m,n,prec,fixed)) return false;
        factor_shape s{chol,upper,m,n,iMlaenv_mpfr(1,chol?"Cpotrf":"Cgetrf",chol?uplo:" ",m,chol?-1:n,-1,-1)};
        std::unique_ptr<long[]> piv(new long[std::min(m,n)]); long result=0; factor_times times; bool success=false;
        if(fixed) switch(prec) {
#define DISPATCH(PB) case PB: success=call<factor_fixed_tag<PB>>(s,prec,a,lda,piv.get(),result,times); break;
            MPLAPACK_MPFR_CUDA_FIXED_PRECISIONS(DISPATCH)
#undef DISPATCH
        } else {
#define RT(NL) if(prec<=64*NL) success=call<factor_runtime_tag<NL>>(s,prec,a,lda,piv.get(),result,times); else
            MPLAPACK_MPFR_CUDA_RT_TIERS(RT)
#undef RT
            return false;
        }
        if(!success) return false;
        if(!chol) for(long i=0;i<std::min(m,n);++i) ipiv[i]=piv[i];
        info=result;
        if(env_long("MPLAPACK_MPFR_CUDA_VERBOSE",0))
            std::fprintf(stderr,"%s MPLAPACK_MPFR CUDA p=%ld m=%ld n=%ld upload=%.9g factor=%.9g download=%.9g\n",chol?"Cpotrf":"Cgetrf",prec,m,n,times.upload,times.factor,times.download);
        return true;
    } catch(const std::bad_alloc &) { return false; }
      catch(const std::length_error &) { return false; }
}
}
bool Cgetrf_mpfr_cuda(mplapackint const m,mplapackint const n,mpc_class *a,mplapackint const lda,mplapackint *ipiv,mplapackint &info) {
    return factor_call(false,nullptr,m,n,a,lda,ipiv,info);
}
bool Cpotrf_mpfr_cuda(const char *uplo,mplapackint const n,mpc_class *a,mplapackint const lda,mplapackint &info) {
    return factor_call(true,uplo,n,n,a,lda,nullptr,info);
}
#define MPLAPACK_SOLVER_COMPLEX_ONLY
#include "Rsolver_bridge_body.h"
