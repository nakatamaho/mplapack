// Focused public-dispatcher tests: default identity and opt-in bounds/repeats.
#include <mpblas.h>
#include <mplapack_arithmetic_params.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <random>
#define Rdot test_Rdot_reference
#include "../../../reference/Rdot.cpp"
#undef Rdot
#define Cdotc test_Cdotc_reference
#include "../../../reference/Cdotc.cpp"
#undef Cdotc
#define Cdotu test_Cdotu_reference
#include "../../../reference/Cdotu.cpp"
#undef Cdotu
#define Rasum test_Rasum_reference
#include "../../../reference/Rasum.cpp"
#undef Rasum
#define RCasum test_RCasum_reference
#include "../../../reference/RCasum.cpp"
#undef RCasum
#define Rnrm2 test_Rnrm2_reference
#include "../../../reference/Rnrm2.cpp"
#undef Rnrm2
#define RCnrm2 test_RCnrm2_reference
#include "../../../reference/RCnrm2.cpp"
#undef RCnrm2
#define iRamax test_iRamax_reference
#include "../../../reference/iRamax.cpp"
#undef iRamax
#define iCamax test_iCamax_reference
#include "../../../reference/iCamax.cpp"
#undef iCamax
namespace {
long checks=0,failures=0;
bool same(const REAL &a,const REAL &b) {
#if defined(MPLAPACK_BUILD_WITH_MPFR)
    return mpfr_equal_p(a.get_mpfr_t(),b.get_mpfr_t()) && mpfr_signbit(a.get_mpfr_t())==mpfr_signbit(b.get_mpfr_t());
#elif defined(MPLAPACK_BUILD_WITH_GMP)
    return a==b;
#elif defined(MPLAPACK_BUILD_WITH_BINARY80)
    return a==b && std::signbit(a)==std::signbit(b);
#else
    return std::memcmp(&a,&b,sizeof a)==0;
#endif
}
bool same(const COMPLEX &a,const COMPLEX &b) { return same(a.real(),b.real()) && same(a.imag(),b.imag()); }
void check(bool ok,const char *what) { ++checks;if(!ok) { ++failures;std::printf("FAIL %s\n",what); } }
template<class F,class X> void exercise(F call,const X &ref,const REAL &bound,const char *name) {
    unsetenv("MPLAPACK_MPBLAS_FAST");check(same(call(),ref),name);
    setenv("MPLAPACK_MPBLAS_FAST","0",1);check(same(call(),ref),name);
    setenv("MPLAPACK_MPBLAS_FAST","1",1);X a=call(),b=call();
    check(same(a,b),"repeat");
    REAL absolute_floor=0.0;
#if defined(MPLAPACK_BUILD_WITH_DD) || defined(MPLAPACK_BUILD_WITH_TD) || defined(MPLAPACK_BUILD_WITH_QD) || defined(MPLAPACK_BUILD_WITH_DOUBLE)
    absolute_floor=REAL(0x1p-1074)*REAL(32.0);
#endif
    check(abs(a-ref)<=bound+absolute_floor,name);
    unsetenv("MPLAPACK_MPBLAS_FAST");check(same(call(),ref),"per-call off");
}
}
int main() {
    std::mt19937_64 rng(20261104);
    const auto &ap=mplapack::get_arithmetic_params<REAL>();
    // BF backends' conservative effective unit roundoff (4 * 2^-p), and
    // MPFR's 2^-p. Using 4*epsilon is conservative for every backend here.
    REAL u=REAL(4.0)*ap.eps;
    setenv("MPLAPACK_BLAS_OMP_MIN","0",1);
    for(long n:{0L,1L,2L,7L,63L,256L,257L,1025L}) for(long inc:{-3L,1L,2L}) {
        size_t size=std::max(1L,1+(n-1)*std::labs(inc));
        std::vector<REAL> x(size),y(size);
        std::vector<COMPLEX> z(size),w(size);
        for(size_t j=0;j<size;++j) {
            x[j]=REAL(double(int(rng()%2001)-1000))/REAL(997.0);
            y[j]=REAL(double(int(rng()%2001)-1000))/REAL(991.0);
            z[j]=COMPLEX(x[j],y[j]);w[j]=COMPLEX(y[j],-x[j]);
        }
        REAL bound=REAL(32.0)*REAL(std::max(1L,n))*REAL(std::max(1L,n))*u;
        exercise([&]{return Rdot(n,x.data(),inc,y.data(),inc);},test_Rdot_reference(n,x.data(),inc,y.data(),inc),bound,"Rdot");
        exercise([&]{return Cdotc(n,z.data(),inc,w.data(),inc);},test_Cdotc_reference(n,z.data(),inc,w.data(),inc),bound,"Cdotc");
        exercise([&]{return Cdotu(n,z.data(),inc,w.data(),inc);},test_Cdotu_reference(n,z.data(),inc,w.data(),inc),bound,"Cdotu");
        exercise([&]{return Rasum(n,x.data(),inc);},test_Rasum_reference(n,x.data(),inc),bound,"Rasum");
        exercise([&]{return RCasum(n,z.data(),inc);},test_RCasum_reference(n,z.data(),inc),bound,"RCasum");
        REAL ref=test_Rnrm2_reference(n,x.data(),inc);
        REAL tol=REAL(16+8*std::ceil(std::log2(double(std::max(2L,2*n)))))*u;
        exercise([&]{return Rnrm2(n,x.data(),inc);},ref,tol*ref,"Rnrm2");
        ref=test_RCnrm2_reference(n,z.data(),inc);
        exercise([&]{return RCnrm2(n,z.data(),inc);},ref,tol*ref,"RCnrm2");
    }
    // Scaling: squaring these numbers without scaling over/underflows binary64.
    for(double magnitude:{0.0,1e300,1e-300}) {
        REAL x[]={REAL(magnitude),REAL(magnitude),REAL(0.0)};
        COMPLEX z[]={COMPLEX(x[0],x[1]),COMPLEX(x[2],x[0])};
        REAL ref=test_Rnrm2_reference(3,x,1);
        exercise([&]{return Rnrm2(3,x,1);},ref,REAL(64.0)*u*ref,"scaled Rnrm2");
        ref=test_RCnrm2_reference(2,z,1);
        exercise([&]{return RCnrm2(2,z,1);},ref,REAL(64.0)*u*ref,"scaled RCnrm2");
    }
    std::printf("BLAS1 fast: %ld checks, %ld failures\n",checks,failures);
    return failures?1:0;
}
