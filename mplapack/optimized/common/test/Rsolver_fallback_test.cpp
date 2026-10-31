// Strict public-dispatcher fallback oracle, independent of the bool API tests.
#include <mpblas.h>
#include <mplapack.h>
#include <cuda_runtime_api.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#if defined(MPLAPACK_BUILD_WITH_DD)
#define PREFIX "MPLAPACK_DD"
#elif defined(MPLAPACK_BUILD_WITH_TD)
#define PREFIX "MPLAPACK_TD"
#elif defined(MPLAPACK_BUILD_WITH_QD)
#define PREFIX "MPLAPACK_QD"
#else
#define PREFIX "MPLAPACK_MPFR"
#endif
extern decltype(Rtrsm) Rtrsm_reference_cpu;
extern decltype(Ctrsm) Ctrsm_reference_cpu;
namespace {
#define Rtrsm Rtrsm_reference_cpu
#define Ctrsm Ctrsm_reference_cpu
#define Rgetrs oracle_Rgetrs
#define Cgetrs oracle_Cgetrs
#define Rpotrs oracle_Rpotrs
#define Cpotrs oracle_Cpotrs
#define Rgesv oracle_Rgesv
#define Cgesv oracle_Cgesv
#define Rposv oracle_Rposv
#define Cposv oracle_Cposv
#include "../../../reference/Rgetrs.cpp"
#include "../../../reference/Cgetrs.cpp"
#include "../../../reference/Rpotrs.cpp"
#include "../../../reference/Cpotrs.cpp"
#include "../../../reference/Rgesv.cpp"
#include "../../../reference/Cgesv.cpp"
#include "../../../reference/Rposv.cpp"
#include "../../../reference/Cposv.cpp"
#undef Rtrsm
#undef Ctrsm
#undef Rgetrs
#undef Cgetrs
#undef Rpotrs
#undef Cpotrs
#undef Rgesv
#undef Cgesv
#undef Rposv
#undef Cposv
long cases=0,elements=0,errors=0,zero_sign_errors=0;
bool same(const REAL &a,const REAL &b) {
    ++elements;
#if defined(MPLAPACK_BUILD_WITH_MPFR)
    bool equal=mpfr_get_prec(a.get_mpfr_t())==mpfr_get_prec(b.get_mpfr_t()) && mpfr_equal_p(a.get_mpfr_t(),b.get_mpfr_t());
    bool sign=!(mpfr_zero_p(a.get_mpfr_t()) && mpfr_zero_p(b.get_mpfr_t())) || mpfr_signbit(a.get_mpfr_t())==mpfr_signbit(b.get_mpfr_t());
    if(!sign) ++zero_sign_errors;
    return equal && sign;
#else
    bool equal=!std::memcmp(&a,&b,sizeof(REAL));
    if(!equal && a==REAL(0) && b==REAL(0)) ++zero_sign_errors;
    return equal;
#endif
}
bool same(const COMPLEX &a,const COMPLEX &b) { bool r=same(REAL(a.real()),REAL(b.real())); bool i=same(REAL(a.imag()),REAL(b.imag())); return r && i; }
template<class H> H zero(int i);
template<> REAL zero<REAL>(int i) { return REAL(i&1?-0.0:0.0); }
template<> COMPLEX zero<COMPLEX>(int i) { return COMPLEX(REAL(i&1?-0.0:0.0),REAL(i&2?-0.0:0.0)); }
#define SOLVE(H,K) \
void solve(bool oracle,bool driver,bool chol,const char *f,INTEGER n,INTEGER rhs,H *a,INTEGER lda,INTEGER *p,H *b,INTEGER ldb,INTEGER &info) { \
    if(oracle) { if(driver) { if(chol) oracle_##K##posv(f,n,rhs,a,lda,b,ldb,info); else oracle_##K##gesv(n,rhs,a,lda,p,b,ldb,info); } \
        else { if(chol) oracle_##K##potrs(f,n,rhs,a,lda,b,ldb,info); else oracle_##K##getrs(f,n,rhs,a,lda,p,b,ldb,info); } } \
    else { if(driver) { if(chol) K##posv(f,n,rhs,a,lda,b,ldb,info); else K##gesv(n,rhs,a,lda,p,b,ldb,info); } \
        else { if(chol) K##potrs(f,n,rhs,a,lda,b,ldb,info); else K##getrs(f,n,rhs,a,lda,p,b,ldb,info); } } \
}
SOLVE(REAL,R)
SOLVE(COMPLEX,C)
#undef SOLVE
template<class H> void test(bool driver,bool chol,const char *f,INTEGER rhs) {
    constexpr INTEGER n=257,lda=n+2,ldb=n+3;
    std::vector<H> a(lda*n,H(19)),b(ldb*rhs,H(23)); std::vector<INTEGER> p(n+2,777);
    for(INTEGER j=0;j<n;++j) for(INTEGER i=0;i<n;++i) a[i+j*lda]=i==j?H(chol?1:-1):zero<H>(i+3*j);
    for(INTEGER j=0;j<rhs;++j) for(INTEGER i=0;i<n;++i) b[i+j*ldb]=zero<H>(i+3*j);
    for(INTEGER i=0;i<n;++i) p[i]=i+1;
    auto expected=a,expected_b=b; auto expected_p=p; INTEGER info=-777,expected_info=-777;
    setenv(PREFIX "_CUDA","0",1); // Oracle factors use the existing optimized CPU path.
    solve(true,driver,chol,f,n,rhs,expected.data(),lda,expected_p.data(),expected_b.data(),ldb,expected_info);
    setenv(PREFIX "_CUDA","1",1); // Catch accidental re-entry into legacy CUDA BLAS.
    solve(false,driver,chol,f,n,rhs,a.data(),lda,p.data(),b.data(),ldb,info);
    long before=errors;
    for(size_t i=0;i<a.size();++i) if(!same(a[i],expected[i])) ++errors;
    for(size_t i=0;i<b.size();++i) if(!same(b[i],expected_b[i])) ++errors;
    if(info!=expected_info || p!=expected_p) ++errors;
    if(before!=errors) std::fprintf(stderr,"fallback driver=%d chol=%d flag=%s rhs=%ld new_errors=%ld\n",driver,chol,f,long(rhs),errors-before);
    ++cases;
}
}
int main(int argc,char **argv) {
    int devices=0; if(cudaGetDeviceCount(&devices)!=cudaSuccess || !devices) return 77;
    setenv("MPLAPACK_LAPACK_FAST","0",1); setenv(PREFIX "_LAPACK_CUDA","0",1);
    long precision=0;
#if defined(MPLAPACK_BUILD_WITH_MPFR)
    precision=argc>1?std::strtol(argv[1],nullptr,10):512;
    mpfrxx::set_default_precision_bits(precision); mpfrxx::set_default_mpc_precision_bits(precision);
#else
    (void)argc; (void)argv;
#endif
    for(bool driver:{false,true}) for(bool chol:{false,true}) for(INTEGER rhs:{INTEGER(1),INTEGER(3)})
        for(const char *f:chol?std::vector<const char*>{"U","L"}:driver?std::vector<const char*>{"N"}:std::vector<const char*>{"N","T","C"}) {
            test<REAL>(driver,chol,f,rhs); test<COMPLEX>(driver,chol,f,rhs);
        }
    std::printf("fallback p=%ld cases=%ld elements=%ld errors=%ld zero-sign-errors=%ld\n",precision,cases,elements,errors,zero_sign_errors);
    return errors?1:0;
}
