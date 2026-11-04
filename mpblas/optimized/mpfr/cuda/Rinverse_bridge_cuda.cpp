// Public MPFR inverse entries. No CPU computation or fallback in this TU.
#include <mpblas_mpfr_cuda.h>
#include <mplapack_mpfr.h>
#include "mpfr_pack_cuda.h"
#include "Rinverse_dynamic_cuda.h"
#include <chrono>
#include <cstdio>
#include <limits>
#include <memory>
#include <stdexcept>
namespace {
// Keep mutex implementation details out of this library's exported ABI.
// Flatten only this tiny guard, not the numerical/packing call graph.
struct inverse_lock {
    std::mutex &mutex;
    __attribute__((flatten)) explicit inverse_lock(std::mutex &m):mutex(m) {mutex.lock();}
    __attribute__((flatten)) ~inverse_lock() {mutex.unlock();}
    inverse_lock(const inverse_lock &)=delete;
    inverse_lock &operator=(const inverse_lock &)=delete;
};
using namespace mplapack_mpfr_cuda_host;
using namespace mplapack_mpfr_cuda;
using namespace mplapack_mpfr_inverse;
using namespace mplapack_blas_dev;
struct pivot_delete { void operator()(long *p) const { delete[] p; } };
template<int PB> bool pack_real(mpfr_srcptr x,factor_real<PB> &v,long) {
    if(!pack<PB>(x,v.v,RT_DEVICE_EMAX)) return false;
    if(v.v.is_zero()) v.v.sign=mpfr_signbit(x)?-1:1;
    return true;
}
template<int NL> bool pack_real(mpfr_srcptr x,rreal<NL> &v,long p) {
    if(mpfr_get_prec(x)!=p || !mpfr_number_p(x)) return false;
    v.kind=mpfr_custom_get_kind(x); v.exp=0;
    for(int i=0;i<NL;++i) v.m[i]=0;
    if(mpfr_zero_p(x)) return true;
    v.exp=mpfr_get_exp(x);
    if(v.exp < -RT_DEVICE_EMAX || v.exp > RT_DEVICE_EMAX) return false;
    const auto limbs=(const mp_limb_t*)mpfr_custom_get_significand(x);
    for(int i=0;i<(p+63)/64;++i) v.m[i]=limbs[i];
    return true;
}
template<int PB> void unpack_real(const factor_real<PB> &v,mpfr_ptr x,long) {
    if(v.v.is_zero()) mpfr_set_zero(x,v.v.sign); else unpack<PB>(v.v,x);
}
template<int NL> void unpack_real(const rreal<NL> &v,mpfr_ptr x,long p) {
    mpfr_t t; mpfr_custom_init_set(t,v.kind,v.exp,p,const_cast<mp_limb_t*>((const mp_limb_t*)v.m));
    mpfr_set(x,t,MPFR_RNDN);
}
template<class Tag> bool get(const mpfr_class &x,value<Tag,false> &v,long p) { return pack_real(x.get_mpfr_t(),v.v,p); }
template<class Tag> bool get(const mpc_class &x,value<Tag,true> &v,long p) {
    return pack_real(mpc_realref(x.mpc_data()),v.v.re_,p) && pack_real(mpc_imagref(x.mpc_data()),v.v.im_,p);
}
template<class Tag> void put(const value<Tag,false> &v,mpfr_class &x,long p) { unpack_real(v.v,x.get_mpfr_t(),p); }
template<class Tag> void put(const value<Tag,true> &v,mpc_class &x,long p) {
    unpack_real(v.v.re_,mpc_realref(x.mpc_data()),p); unpack_real(v.v.im_,mpc_imagref(x.mpc_data()),p);
}
bool precision(const mpfr_class &x,long p) { return mpfr_get_prec(x.get_mpfr_t())==p; }
bool precision(const mpc_class &x,long p) { return mpfr_get_prec(mpc_realref(x.mpc_data()))==p && mpfr_get_prec(mpc_imagref(x.mpc_data()))==p; }
void integer(mpfr_class &x,long n) { mpfr_set_si(x.get_mpfr_t(),n,MPFR_RNDN); }
void integer(mpc_class &x,long n) { mpc_set_si(x.mpc_data(),n,MPC_RNDNN); }
bool referenced(const shape &s,long i,long j) { return (s.kind==3 || (s.upper?i<=j:i>=j)) && !(s.unit && i==j); }
template<class Tag,bool Complex,class H> bool call(const shape &s,long prec,bool fixed,H *a,long lda,
    H *work,const long *piv,long optimal,long iws,long &info,factor_times &times) {
    using V=value<Tag,Complex>;
    size_t count=size_t(s.n)*s.n,wcount=s.kind==3?size_t(s.n)*s.columns:0;
    if(count>std::numeric_limits<size_t>::max()-wcount || !memory_fits(count+wcount,sizeof(V),s.n,prec,fixed)) return false;
    using Clock=std::chrono::steady_clock; auto start=Clock::now();
    std::vector<V> stage(count),scratch(wcount);
    for(long j=0;j<s.n;++j) for(long i=0;i<s.n;++i) if(referenced(s,i,j))
        if(!get(a[i+j*lda],stage[i+j*s.n],prec)) return false;
    // Work is output-only: inspect precision, never its uninitialized value.
    if(s.kind==3 && !precision(work[0],prec)) return false;
    auto packed=Clock::now();
    long result=0;
    if(run<Tag,Complex>(s,prec,stage.data(),scratch.data(),piv,result,times)) return false;
    for(const auto &v:stage) {
        if constexpr(Complex) { if(!mplapack_mpfr_inverse::fits(v.v.re_) || !mplapack_mpfr_inverse::fits(v.v.im_)) return false; }
        else if(!mplapack_mpfr_inverse::fits(v.v)) return false;
    }
    if(s.kind==3 && !result) for(long j=0;j<s.columns;++j) for(long i=j+1;i<s.n;++i) {
        const auto &v=scratch[i+j*s.n];
        if(!precision(work[i+j*s.n],prec)) return false;
        if constexpr(Complex) { if(!mplapack_mpfr_inverse::fits(v.v.re_) || !mplapack_mpfr_inverse::fits(v.v.im_)) return false; }
        else if(!mplapack_mpfr_inverse::fits(v.v)) return false;
    }
    // All kernels, transfers, state restoration and output validation succeeded.
    // These exact MPFR stores cannot fail; publish A, scratch and info together.
    auto before=Clock::now();
    for(long j=0;j<s.n;++j) for(long i=0;i<s.n;++i) if(referenced(s,i,j)) put(stage[i+j*s.n],a[i+j*lda],prec);
    if(s.kind==3) {
        if(!result) for(long j=0;j<s.columns;++j) for(long i=j+1;i<s.n;++i) put(scratch[i+j*s.n],work[i+j*s.n],prec);
        integer(work[0],result?optimal:iws);
    }
    info=result;
    times.upload+=std::chrono::duration<double>(packed-start).count();
    times.download+=std::chrono::duration<double>(Clock::now()-before).count();
    return true;
}
// Above the owned-limb tiers, real inverses use dynamically sized scalar
// slots and packed matrix views, just as real factors have no fixed bit cap.
bool dynamic_call(const shape &s,long prec,mpfr_class *a,long lda,mpfr_class *work,const long *piv,
    long optimal,long iws,long &info,factor_times &times) {
    if(!dynamic_real_memory_fits(s,prec)) return false;
    using Clock=std::chrono::steady_clock;auto start=Clock::now();
    size_t acount=size_t(s.n)*s.n,wcount=s.kind==3?size_t(s.n)*s.columns:0,count=acount+wcount;
    int nl=(prec+63)/64;rt_storage stage;stage.resize(count,nl);
    std::fill(stage.kind.begin(),stage.kind.end(),RT_ZERO_KIND);
    for(long j=0;j<s.n;++j) for(long i=0;i<s.n;++i) if(referenced(s,i,j))
        if(!rt_pack_value(a[i+j*lda].get_mpfr_t(),prec,nl,stage,i+j*s.n,RT_DEVICE_EMAX)) return false;
    if(s.kind==3 && !precision(work[0],prec)) return false;
    auto packed=Clock::now();long result=0;
    if(run_dynamic_real(s,prec,stage.pack,piv,result,times)) return false;
    for(size_t i=0;i<acount;++i) if(!rt_fits(stage,i)) return false;
    if(s.kind==3 && !result) for(long j=0;j<s.columns;++j) for(long i=j+1;i<s.n;++i)
        if(!precision(work[i+j*s.n],prec) || !rt_fits(stage,acount+i+j*s.n)) return false;
    auto before=Clock::now();
    for(long j=0;j<s.n;++j) for(long i=0;i<s.n;++i) if(referenced(s,i,j))
        rt_unpack(stage,i+j*s.n,prec,nl,a[i+j*lda].get_mpfr_t());
    if(s.kind==3) {
        if(!result) for(long j=0;j<s.columns;++j) for(long i=j+1;i<s.n;++i)
            rt_unpack(stage,acount+i+j*s.n,prec,nl,work[i+j*s.n].get_mpfr_t());
        integer(work[0],result?optimal:iws);
    }
    info=result;
    times.upload+=std::chrono::duration<double>(packed-start).count();
    times.download+=std::chrono::duration<double>(Clock::now()-before).count();
    return true;
}
char flag(const char *s) { char c=s?*s:0; return c>='a' && c<='z'?c-32:c; }
template<class H> bool inverse_call(const char *name,int kind,const char *uplo,const char *diag,
    mplapackint n,H *a,mplapackint lda,mplapackint *ipiv,H *work,mplapackint lwork,mplapackint &info) {
    constexpr bool complex=std::is_same<H,mpc_class>::value;
    char u=flag(uplo),d=flag(diag);
    if(n<=0 || lda<n || !a || (kind<3 && u!='U' && u!='L') || (kind==0 && d!='U' && d!='N') ||
       (kind==3 && (!ipiv || !work || lwork<n))) return false;
    long minimum=env_long("MPLAPACK_MPFR_LAPACK_CUDA_MIN",std::numeric_limits<long>::max());
    minimum=env_long("MPLAPACK_MPFR_LAPACK_INVERSE_CUDA_MIN",minimum);
    if(minimum<0 || n<minimum || !env_long("MPLAPACK_MPFR_CUDA",1) || !env_long("MPLAPACK_MPFR_LAPACK_CUDA",1)) return false;
    long prec=mpfrxx::default_precision_bits();
    if(prec<MPFR_PREC_MIN || prec>std::numeric_limits<int>::max()-1024L || (complex && prec>4096) || mpfr_get_default_rounding_mode()!=MPFR_RNDN ||
       mpfr_get_emax()!=RT_DEVICE_EMAX || mpfr_get_emin()!=-RT_DEVICE_EMAX) return false;
    if(complex && (mpfrxx::default_mpc_real_precision_bits()!=prec || mpfrxx::default_mpc_imag_precision_bits()!=prec ||
       mpfrxx::default_mpc_rounding_mode()!=MPC_RNDNN)) return false;
    size_t max=std::numeric_limits<size_t>::max();
    if(size_t(n)>max/size_t(n)/(2*((prec+63)/64*8+32)) || size_t(n)>max/size_t(lda)/sizeof(H)) return false;
    try {
        inverse_lock lock(gpu_mutex());
        bool fixed=false;
        if(!env_long("MPLAPACK_MPFR_CUDA_FORCE_RUNTIME",0)) switch(prec) {
#define FIXED(PB) case PB: fixed=true; break;
            MPLAPACK_MPFR_CUDA_FIXED_PRECISIONS(FIXED)
#undef FIXED
        }
        if(!fixed && !env_long("MPLAPACK_MPFR_CUDA_RUNTIME",1)) return false;
        shape s; s.kind=kind;s.upper=kind==3 || u=='U';s.unit=kind==0 && d=='U';s.n=n;
        char opts[3]={s.upper?'U':'L',s.unit?'U':'N',0};
        s.nbtri=iMlaenv_mpfr(1,complex?"Ctrtri":"Rtrtri",opts,n,-1,-1,-1);
        s.nblau=iMlaenv_mpfr(1,complex?"Clauum":"Rlauum",opts,n,-1,-1,-1);
        s.nb=kind==3?iMlaenv_mpfr(1,name," ",n,-1,-1,-1):1;
        if(s.nb>0 && n>std::numeric_limits<long>::max()/s.nb) return false;
        long optimal=std::max(1L,long(n)*s.nb),iws=n;
        if(kind==3 && s.nb>1 && s.nb<n) {
            iws=long(n)*s.nb;
            if(lwork<iws) { s.nb=lwork/n;s.nbmin=std::max(2L,long(iMlaenv_mpfr(2,name," ",n,-1,-1,-1))); }
        }
        s.columns=kind==3 && s.nb>=s.nbmin && s.nb<n?s.nb:1;
        std::unique_ptr<long[],pivot_delete> piv(kind==3?new long[n]:nullptr);
        if(kind==3) for(long i=0;i<n;++i) { if(ipiv[i]<1 || ipiv[i]>n) return false;piv[i]=ipiv[i]; }
        factor_times times; long result=0; bool success=false;
        if(prec>4096) {
            if constexpr(!complex) success=dynamic_call(s,prec,a,lda,work,piv.get(),optimal,iws,result,times);
        } else if(fixed) switch(prec) {
#define DISPATCH(PB) case PB: success=call<factor_fixed_tag<PB>,complex>(s,prec,true,a,lda,work,piv.get(),optimal,iws,result,times);break;
            MPLAPACK_MPFR_CUDA_FIXED_PRECISIONS(DISPATCH)
#undef DISPATCH
        } else {
#define RT(NL) if(prec<=64*NL) success=call<factor_runtime_tag<NL>,complex>(s,prec,false,a,lda,work,piv.get(),optimal,iws,result,times);else
            MPLAPACK_MPFR_CUDA_RT_TIERS(RT)
#undef RT
            return false;
        }
        if(!success) return false;
        info=result;
        if(env_long("MPLAPACK_MPFR_CUDA_VERBOSE",0)) std::fprintf(stderr,"%s MPLAPACK_MPFR CUDA p=%ld n=%ld upload=%.9g inverse=%.9g download=%.9g\n",name,prec,long(n),times.upload,times.factor,times.download);
        return true;
    } catch(const std::bad_alloc &) { return false; } catch(const std::length_error &) { return false; }
}
}
#define ENTRIES(H,P) \
bool P##trtri_mpfr_cuda(const char *u,const char *d,mplapackint const n,H *a,mplapackint const lda,mplapackint &info) { return inverse_call(#P "trtri",0,u,d,n,a,lda,nullptr,(H*)nullptr,0,info); } \
bool P##lauum_mpfr_cuda(const char *u,mplapackint const n,H *a,mplapackint const lda,mplapackint &info) { return inverse_call(#P "lauum",1,u,"N",n,a,lda,nullptr,(H*)nullptr,0,info); } \
bool P##potri_mpfr_cuda(const char *u,mplapackint const n,H *a,mplapackint const lda,mplapackint &info) { return inverse_call(#P "potri",2,u,"N",n,a,lda,nullptr,(H*)nullptr,0,info); } \
bool P##getri_mpfr_cuda(mplapackint const n,H *a,mplapackint const lda,mplapackint *p,H *work,mplapackint const lwork,mplapackint &info) { return inverse_call(#P "getri",3,"U","N",n,a,lda,p,work,lwork,info); }
ENTRIES(mpfr_class,R)
ENTRIES(mpc_class,C)
#undef ENTRIES
