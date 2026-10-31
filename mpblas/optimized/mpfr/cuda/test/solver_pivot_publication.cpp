// Deterministic pivot-DMA ordering test, plus an expected-failing no-fence control.
#include <mpblas_mpfr_cuda.h>
#include <cuda_runtime_api.h>
#include <cstdio>
#include <cstdlib>
#include <dlfcn.h>
#include <vector>
namespace {
bool same(const mpfr_class &a,const mpfr_class &b) {
    auto x=a.get_mpfr_t(),y=b.get_mpfr_t();
    return mpfr_equal_p(x,y) && (!mpfr_zero_p(x) || mpfr_signbit(x)==mpfr_signbit(y));
}
bool same(const mpc_class &a,const mpc_class &b) {
    auto x=a.mpc_data(),y=b.mpc_data();
    auto equal=[](mpfr_srcptr u,mpfr_srcptr v) { return mpfr_equal_p(u,v) && (!mpfr_zero_p(u) || mpfr_signbit(u)==mpfr_signbit(v)); };
    return equal(mpc_realref(x),mpc_realref(y)) && equal(mpc_imagref(x),mpc_imagref(y));
}
void negative_zero(mpfr_class &x) { mpfr_set_zero(x.get_mpfr_t(),-1); }
void negative_zero(mpc_class &x) { mpfr_set_zero(mpc_realref(x.mpc_data()),-1); }
bool entry(mpfr_class *a,mplapackint *p,mpfr_class *b,mplapackint &info) { return Rgetrs_mpfr_cuda("N",65,1,a,65,p,b,65,info); }
bool entry(mpc_class *a,mplapackint *p,mpc_class *b,mplapackint &info) { return Cgetrs_mpfr_cuda("N",65,1,a,65,p,b,65,info); }
template<class H> bool test(const char *name) {
    constexpr long n=65;
    std::vector<H> a(n*n,H(0)),b(n); std::vector<mplapackint> p(n);
    for(long i=0;i<n;++i) { a[i+i*n]=H(1); b[i]=H(i+1); p[i]=i+1; }
    negative_zero(b[0]);
    { // Load the kernels before starting the delayed-copy window.
        auto warm_a=a,warm_b=b; auto warm_p=p; mplapackint warm_info=-1;
        if(!entry(warm_a.data(),warm_p.data(),warm_b.data(),warm_info) || warm_info) return false;
    }
    p[0]=n;
    auto expected=b,original_a=a; auto original_p=p;
    expected[0]=b[n-1]; expected[n-1]=b[0];
    mplapackint info=-987; bool computed=entry(a.data(),p.data(),b.data(),info);
    long mismatches=0; bool unchanged=true;
    for(long i=0;i<n;++i) if(!same(b[i],expected[i])) ++mismatches;
    for(long i=0;i<n*n;++i) if(!same(a[i],original_a[i])) unchanged=false;
    std::printf("%s computed=%d info=%ld B_mismatches=%ld A_untouched=%d pivots_untouched=%d\n",
        name,computed?1:0,long(info),mismatches,unchanged?1:0,p==original_p?1:0);
    return computed && info==0 && !mismatches && unchanged && p==original_p;
}
int count(const char *name) { auto f=reinterpret_cast<int(*)()>(dlsym(RTLD_DEFAULT,name)); return f?f():-1; }
}
int main(int argc,char **argv) {
    int devices=0; if(cudaGetDeviceCount(&devices)!=cudaSuccess || !devices) return 77;
    if(argc!=3) return 2;
    long precision=std::atol(argv[1]); bool control=std::atoi(argv[2])!=0;
    mpfrxx::set_default_precision_bits(precision); mpfrxx::set_default_mpc_precision_bits(precision);
    setenv("MPLAPACK_MPFR_CUDA","1",1); setenv("MPLAPACK_MPFR_LAPACK_CUDA","1",1);
    setenv("MPLAPACK_MPFR_LAPACK_CUDA_MIN","0",1); setenv("MPLAPACK_MPFR_LAPACK_SOLVE_CUDA_MIN","0",1);
    setenv("MPLAPACK_MPFR_LAPACK_SOLVE_CUDA_MIN_RHS","0",1); setenv("MPLAPACK_MPFR_CUDA_FORCE_RUNTIME","0",1);
    bool real_ok=test<mpfr_class>("Rgetrs"),complex_ok=test<mpc_class>("Cgetrs");
    int copies=count("mplapack_test_pivot_injections"),waits=count("mplapack_test_pivot_waits"),skips=count("mplapack_test_pivot_skips");
    std::printf("p=%ld copies=%d fences=%d skipped=%d control=%d\n",precision,copies,waits,skips,control?1:0);
    auto release=reinterpret_cast<void(*)()>(dlsym(RTLD_DEFAULT,"mplapack_test_pivot_release")); if(release) release();
    if(control) {
        // CTest WILL_FAIL: success here means the deliberately broken ordering
        // actually produced a mismatch, not merely that the shim was absent.
        return copies==2 && skips==2 && (!real_ok || !complex_ok)?1:0;
    }
    return real_ok && complex_ok && copies==2 && waits==2 && skips==0?0:1;
}
