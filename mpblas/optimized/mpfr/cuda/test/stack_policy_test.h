// Device-only, non-vacuous stack policy probes shared by MPFR CUDA tests.
#ifndef MPLAPACK_MPFR_STACK_POLICY_TEST_H
#define MPLAPACK_MPFR_STACK_POLICY_TEST_H
#include "../mpfr_cuda_stack.h"
#include <mpblas_mpfr_cuda.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <utility>
#ifdef MPLAPACK_MPFR_STACK_TRACE
#include "stack_limit_trace.h"
#endif
namespace mpfr_stack_test {
inline void require(bool ok,const char *what) {
    if(!ok) { std::fprintf(stderr,"FAIL stack policy: %s\n",what); std::exit(1); }
}
inline bool equal(mpfr_srcptr a,mpfr_srcptr b) {
    if(mpfr_get_prec(a)!=mpfr_get_prec(b) || mpfr_signbit(a)!=mpfr_signbit(b)) return false;
    if(mpfr_zero_p(a) || mpfr_zero_p(b)) return mpfr_zero_p(a) && mpfr_zero_p(b);
    if(!mpfr_regular_p(a) || !mpfr_regular_p(b)) return false;
    return mpfr_get_exp(a)==mpfr_get_exp(b) &&
        !std::memcmp(mpfr_custom_get_significand(a),mpfr_custom_get_significand(b),
                     ((mpfr_get_prec(a)+GMP_NUMB_BITS-1)/GMP_NUMB_BITS)*sizeof(mp_limb_t));
}
inline bool equal(const mpfr_class &a,const mpfr_class &b) { return equal(a.get_mpfr_t(),b.get_mpfr_t()); }
inline bool equal(const mpc_class &a,const mpc_class &b) {
    return equal(mpc_realref(a.mpc_data()),mpc_realref(b.mpc_data())) &&
           equal(mpc_imagref(a.mpc_data()),mpc_imagref(b.mpc_data()));
}
template<class T> struct result {
    bool computed=false;
    mplapackint info=-77;
    std::vector<T> a,b,work;
    std::vector<mplapackint> piv;
    bool operator==(const result &r) const {
        if(computed!=r.computed || info!=r.info || piv!=r.piv) return false;
        for(auto pair:{std::make_pair(&a,&r.a),std::make_pair(&b,&r.b),std::make_pair(&work,&r.work)}) {
            if(pair.first->size()!=pair.second->size()) return false;
            for(size_t i=0;i<pair.first->size();++i) if(!equal((*pair.first)[i],(*pair.second)[i])) return false;
        }
        return true;
    }
};
inline bool fixed(long p) {
    const char *v=std::getenv("MPLAPACK_MPFR_CUDA_FORCE_RUNTIME");
    return !(v && v[0] && v[0]!='0') && (p==256 || p==512 || p==768 || p==1024 || p==2048);
}
inline size_t factor_bytes(long p,bool complex) {
    using namespace mplapack_mpfr_cuda;
    return complex ? (fixed(p)?complex_fixed_stack_bytes(p):complex_runtime_stack_bytes(p)) :
                     (fixed(p)?fixed_stack_bytes(p):runtime_stack_bytes(p));
}
inline size_t inverse_bytes(long p) {
    using namespace mplapack_mpfr_cuda;
    return 2*(fixed(p)?complex_fixed_stack_bytes(p):complex_runtime_stack_bytes(p));
}
inline bool restore_mode() {
    const char *v=std::getenv("MPLAPACK_MPFR_CUDA_STACK_RESTORE");
    return v && !std::strcmp(v,"1");
}
inline void begin_sweep(const char *family,bool expected) {
    require(restore_mode()==expected,"probe preserves caller policy for numerical sweep");
    std::printf("begin numerical family %s STACK_RESTORE=%s\n",family,expected?"1":"default");
    std::fflush(stdout);
}
// Inspect a real numerical sweep call without changing its environment.
// Only the test fixture changes the initial limit, never the entry's policy.
struct call_limit {
    std::string name;
    size_t original=0,small=0,required;
    bool restore;
    call_limit(const std::string &label,size_t bytes):name(label),required(bytes) {
        const char *v=std::getenv("MPLAPACK_MPFR_CUDA_STACK_RESTORE");
        restore=v && !std::strcmp(v,"1");
        require(!cudaDeviceGetLimit(&original,cudaLimitStackSize),"sweep read original");
        require(!cudaDeviceSetLimit(cudaLimitStackSize,1024),"sweep set small");
        require(!cudaDeviceGetLimit(&small,cudaLimitStackSize) && small<required,"sweep small < required");
#ifdef MPLAPACK_MPFR_STACK_TRACE
        mpfr_stack_trace::begin(small);
#endif
    }
    void finish(bool computed,mplapackint info=0,bool should_compute=true) {
        size_t after=0;
        require(!cudaDeviceGetLimit(&after,cudaLimitStackSize),"sweep read after");
        std::printf("numerical sweep %s mode=%s small=%zu required=%zu after=%zu computed=%d info=%ld\n",
                    name.c_str(),restore?"restore":"raise-only",small,required,after,int(computed),long(info));
#ifdef MPLAPACK_MPFR_STACK_TRACE
        auto trace=mpfr_stack_trace::end();
        std::printf("stack API transitions peak=%zu raises=%u lowers=%u\n",trace.peak,trace.raises,trace.lowers);
        if(should_compute) {
            require(trace.peak>=required && trace.raises>0,"real entry raised stack (API trace)");
            require(restore ? trace.lowers>0 : trace.lowers==0,"real entry restore policy (API trace)");
        } else require(!trace.raises && !trace.lowers,"refused entry does not change stack (API trace)");
#endif
        std::fflush(stdout);
        const bool unchanged=restore || !should_compute;
        require(unchanged ? after==small : after>small && after>=required,
                unchanged ? "sweep restore/unchanged == small" : "sweep raise > small and >= required");
        require(computed==should_compute,"sweep computation status");
    }
    ~call_limit() { require(!cudaDeviceSetLimit(cudaLimitStackSize,original),"sweep restore fixture"); }
};
// Run the SAME actual numerical trial (not the 3x3 eligibility probe) under
// both policies. Keep all padded outputs and info for exact comparison.
template<class Call> auto numerical_policies(const char *name,size_t required,Call run,bool computed=true) {
    const char *key="MPLAPACK_MPFR_CUDA_STACK_RESTORE",*v=std::getenv(key);
    bool had=v!=nullptr; std::string saved=v?v:"";
    decltype(run()) baseline;
    for(int mode=0;mode<2;++mode) {
        if(mode) setenv(key,"1",1); else unsetenv(key);
        call_limit limit(name,required);
        auto got=run(); limit.finish(got.computed,got.info,computed);
        if(!mode) baseline=got;
        else require(got==baseline,"numerical sweep A/B/work/ipiv/info limb/sign identity");
    }
    if(had) setenv(key,saved.c_str(),1);else unsetenv(key);
    return baseline;
}
// Save the caller's environment and stack: CTest's =1 variant must still run
// its remaining numerical sweep in restore mode, not accidentally unset it.
template<class Call> void policies(const char *name,size_t required,Call run) {
    const char *key="MPLAPACK_MPFR_CUDA_STACK_RESTORE",*v=std::getenv(key);
    bool had=v!=nullptr; std::string saved=v?v:"";
    size_t original=0;
    require(!cudaDeviceGetLimit(&original,cudaLimitStackSize),"read original");
    decltype(run()) baseline;
    // Check exact parsing as well as per-call switching (no static cache).
    for(int mode=0;mode<4;++mode) {
        if(mode==0) unsetenv(key);
        else setenv(key,mode==1?"1":mode==2?"0":"other",1);
        require(!cudaDeviceSetLimit(cudaLimitStackSize,1024),"set small");
        size_t small=0,after=0;
        require(!cudaDeviceGetLimit(&small,cudaLimitStackSize) && small<required,"small < required (non-vacuous)");
        auto got=run();
        require(got.computed,"entry computed");
        require(!cudaDeviceGetLimit(&after,cudaLimitStackSize),"read after");
        std::printf("stack policy %s mode=%s small=%zu required=%zu after=%zu\n",name,
                    mode==0?"unset":mode==1?"1":mode==2?"0":"other",small,required,after);
        std::fflush(stdout);
        require(mode==1 ? after==small : after>small && after>=required,
                mode==1?"restore == small":"raise > small and >= required");
        if(mode==0) baseline=got;
        else require(got==baseline,"A/B/work/ipiv/info limb/sign identity across policies");
    }
    require(!cudaDeviceSetLimit(cudaLimitStackSize,original),"restore test's original stack");
    if(had) setenv(key,saved.c_str(),1); else unsetenv(key);
}
}
#endif
