// Direct entries against CPU-only library: no dispatcher can mask failure.
#if defined(__CUDACC__)
#define MPLAPACK_MINMAX_MPLAPACKINT_DEFINED
#endif
#include <mpblas.h>
#include <mplapack.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <type_traits>
#include <vector>
#if defined(MPLAPACK_BUILD_WITH_DD)
#include <mpblas_dd_cuda.h>
#define ENTRY(R) R##_dd_cuda
#define PREFIX "MPLAPACK_DD"
#elif defined(MPLAPACK_BUILD_WITH_TD)
#include <mpblas_td_cuda.h>
#define ENTRY(R) R##_td_cuda
#define PREFIX "MPLAPACK_TD"
#elif defined(MPLAPACK_BUILD_WITH_QD)
#include <mpblas_qd_cuda.h>
#define ENTRY(R) R##_qd_cuda
#define PREFIX "MPLAPACK_QD"
#else
#include <mpblas_mpfr_cuda.h>
#define ENTRY(R) R##_mpfr_cuda
#define PREFIX "MPLAPACK_MPFR"
#endif
#if !defined(MPLAPACK_BUILD_WITH_MPFR)
#include "../Rsolver_cuda_entry.cuh"
#endif
#if defined(MPLAPACK_SOLVER_DEVICE_TEST)
#include <cuda_runtime_api.h>
#endif
namespace {
long checks=0,comparisons=0,zero_sign_errors=0;
void check(bool ok,const char *what) { ++checks; if(!ok) { std::fprintf(stderr,"FAIL %s check=%ld\n",what,checks); std::exit(1); } }
#if defined(MPLAPACK_BUILD_WITH_MPFR)
bool same(const REAL &a,const REAL &b) {
    ++comparisons;
    bool signs=(a!=REAL(0) || b!=REAL(0) || mpfr_signbit(a.get_mpfr_t())==mpfr_signbit(b.get_mpfr_t()));
    if(!signs) ++zero_sign_errors;
    return a==b && signs;
}
bool same(const COMPLEX &a,const COMPLEX &b) { return same(REAL(a.real()),REAL(b.real())) && same(REAL(a.imag()),REAL(b.imag())); }
#else
bool same(const REAL &a,const REAL &b) { ++comparisons; return !std::memcmp(&a,&b,sizeof(a)); }
bool same(const COMPLEX &a,const COMPLEX &b) { return same(a.real(),b.real()) && same(a.imag(),b.imag()); }
#endif
template<class H> bool same(const std::vector<H> &a,const std::vector<H> &b) {
    bool ok=true; for(size_t i=0;i<a.size();++i) if(!same(a[i],b[i])) {
        if(ok) std::fprintf(stderr,"first differing element %zu\n",i); ok=false;
    } return ok;
}
REAL negative_zero() {
    REAL x=-0.0;
#if defined(MPLAPACK_BUILD_WITH_MPFR)
    mpfr_set_zero(x.get_mpfr_t(),-1);
#endif
    return x;
}
REAL value(long i,long j,REAL *) { return REAL((i*17+j*31+i*j*7)%127-63)/REAL(128); }
COMPLEX value(long i,long j,COMPLEX *) { return COMPLEX(value(i,j,(REAL*)nullptr),value(j+3,i+5,(REAL*)nullptr)); }
#define API(H,K) \
void cpu_factor(bool chol,const char *u,INTEGER n,H *a,INTEGER ld,INTEGER *p,INTEGER &info) { if(chol) K##potrf(u,n,a,ld,info); else K##getrf(n,n,a,ld,p,info); } \
void cpu_solve(bool driver,bool chol,const char *f,INTEGER n,INTEGER rhs,H *a,INTEGER ld,INTEGER *p,H *b,INTEGER ldb,INTEGER &info) { \
    if(driver) { if(chol) K##posv(f,n,rhs,a,ld,b,ldb,info); else K##gesv(n,rhs,a,ld,p,b,ldb,info); } \
    else { if(chol) K##potrs(f,n,rhs,a,ld,b,ldb,info); else K##getrs(f,n,rhs,a,ld,p,b,ldb,info); } } \
bool gpu_solve(bool driver,bool chol,const char *f,INTEGER n,INTEGER rhs,H *a,INTEGER ld,INTEGER *p,H *b,INTEGER ldb,INTEGER &info) { \
    if(driver) { if(chol) return ENTRY(K##posv)(f,n,rhs,a,ld,b,ldb,info); return ENTRY(K##gesv)(n,rhs,a,ld,p,b,ldb,info); } \
    if(chol) return ENTRY(K##potrs)(f,n,rhs,a,ld,b,ldb,info); return ENTRY(K##getrs)(f,n,rhs,a,ld,p,b,ldb,info); }
API(REAL,R)
API(COMPLEX,C)
#undef API
template<class H> void test(INTEGER n,INTEGER rhs,bool driver,bool chol,const char *flag,int kind) {
    INTEGER lda=n+3,ldb=n+5,ir=-1,ig=-2;
    std::vector<H> a(size_t(lda)*n,H(19)),b(size_t(ldb)*rhs,H(23));
    for(long j=0;j<n;++j) for(long i=0;i<n;++i) {
        a[i+j*lda]=chol?(i==j?H(2*n):H(0)):value(i,j,(H*)nullptr);
        if(!chol && (kind==0 || kind==4) && i==j) a[i+j*lda]+=H(2*n);
        if(!chol && kind==1) a[i+j*lda]=H(i==n-1-j?4:(i==j?1:0));
        if(!chol && kind==2) a[i+j*lda]=H(i>=j?1:-1); // ties
        if(!chol && kind==3 && j==n/2) a[i+j*lda]=H(0); // singular
        if(chol && i!=j) a[i+j*lda]=H(REAL((i+j)%7-3)/REAL(64));
        if(chol && kind==3 && i==j && i==n/2) a[i+j*lda]=H(-1);
    }
    for(long j=0;j<rhs;++j) for(long i=0;i<n;++i) b[i+j*ldb]=kind==4?H(0):value(i+11,j+7,(H*)nullptr);
    if(kind==4) b[0]=H(negative_zero());
    std::vector<INTEGER> p(n+3,-7);
    if(!driver) {
        cpu_factor(chol,flag,n,a.data(),lda,p.data(),ir);
        if(ir) std::fprintf(stderr,"oracle n=%ld chol=%d flag=%s kind=%d info=%ld\n",long(n),chol,flag,kind,long(ir));
        check(ir==0,"oracle factors");
    }
    auto ga=a,gb=b; auto q=p;
#if defined(MPLAPACK_SOLVER_DEVICE_TEST) && defined(MPLAPACK_BUILD_WITH_MPFR)
    size_t stack=0,after=0,free_before=0,total_before=0;
    check(cudaDeviceGetLimit(&stack,cudaLimitStackSize)==cudaSuccess,"stack before");
#endif
    cpu_solve(driver,chol,flag,n,rhs,a.data(),lda,p.data(),b.data(),ldb,ir);
#if defined(MPLAPACK_SOLVER_DEVICE_TEST) && defined(MPLAPACK_BUILD_WITH_MPFR)
    cudaMemGetInfo(&free_before,&total_before);
#endif
    bool computed=gpu_solve(driver,chol,flag,n,rhs,ga.data(),lda,q.data(),gb.data(),ldb,ig);
    if(!computed) std::fprintf(stderr,"entry rejected type=%s n=%ld rhs=%ld driver=%d chol=%d flag=%s kind=%d reference_info=%ld\n",
        std::is_same<H,REAL>::value?"real":"complex",long(n),long(rhs),driver,chol,flag,kind,long(ir));
#if defined(MPLAPACK_SOLVER_DEVICE_TEST) && defined(MPLAPACK_BUILD_WITH_MPFR)
    if(!computed) {
        size_t free_after=0,total_after=0; auto error=cudaMemGetInfo(&free_after,&total_after);
        std::fprintf(stderr,"device rejection free_before=%zu free_after=%zu memory_query=%s last_error=%s\n",
            free_before,free_after,cudaGetErrorString(error),cudaGetErrorString(cudaPeekAtLastError()));
    }
#endif
    check(computed,"entry computed");
    if(ir!=ig || p!=q || !same(a,ga) || !same(b,gb)) {
        std::fprintf(stderr,"n=%ld rhs=%ld driver=%d chol=%d flag=%s kind=%d info=%ld/%ld\n",long(n),long(rhs),driver,chol,flag,kind,long(ir),long(ig));
        check(false,"A/B/ipiv/info bit identity");
    }
    check(true,"A/B/ipiv/info bit identity");
#if defined(MPLAPACK_SOLVER_DEVICE_TEST) && defined(MPLAPACK_BUILD_WITH_MPFR)
    check(cudaDeviceGetLimit(&after,cudaLimitStackSize)==cudaSuccess && stack==after,"scoped stack restoration");
#endif
}
template<class H> void reject() {
    std::vector<H> a(16,H(7)),b(16,H(11)),aa=a,bb=b;
    INTEGER p[4]={1,2,3,4},info=123;
    auto untouched=[&] { check(info==123 && p[0]==1 && p[1]==2 && p[2]==3 && p[3]==4 && same(a,aa) && same(b,bb),"false leaves outputs untouched"); };
    for(bool driver:{false,true}) for(bool chol:{false,true}) {
        const char *f=chol?"U":"N";
        for(INTEGER n:{INTEGER(-1),INTEGER(0)}) { check(!gpu_solve(driver,chol,f,n,1,a.data(),4,p,b.data(),4,info),"n quick/invalid"); untouched(); }
        for(INTEGER rhs:{INTEGER(-1),INTEGER(0)}) { check(!gpu_solve(driver,chol,f,4,rhs,a.data(),4,p,b.data(),4,info),"nrhs quick/invalid"); untouched(); }
        check(!gpu_solve(driver,chol,f,4,1,a.data(),3,p,b.data(),4,info),"lda");
        check(!gpu_solve(driver,chol,f,4,1,a.data(),4,p,b.data(),3,info),"ldb");
        check(!gpu_solve(driver,chol,f,4,1,(H*)nullptr,4,p,b.data(),4,info),"null");
        INTEGER huge=std::numeric_limits<INTEGER>::max();
        check(!gpu_solve(driver,chol,f,huge,huge,a.data(),huge,p,b.data(),huge,info),"overflow"); untouched();
        if(chol || !driver) {
            check(!gpu_solve(driver,chol,"X",4,1,a.data(),4,p,b.data(),4,info),"flag"); untouched();
        }
        if(!driver && !chol) {
            p[0]=0; check(!gpu_solve(driver,chol,f,4,1,a.data(),4,p,b.data(),4,info),"pivot range"); p[0]=1; untouched();
        }
        setenv(driver?PREFIX "_LAPACK_DRIVER_CUDA_MIN":PREFIX "_LAPACK_SOLVE_CUDA_MIN","5",1);
        check(!gpu_solve(driver,chol,f,4,1,a.data(),4,p,b.data(),4,info),"solver-specific threshold"); untouched();
        setenv(driver?PREFIX "_LAPACK_DRIVER_CUDA_MIN":PREFIX "_LAPACK_SOLVE_CUDA_MIN","0",1);
        if(!driver) {
            setenv(PREFIX "_LAPACK_SOLVE_CUDA_MIN_RHS","5",1);
            check(!gpu_solve(driver,chol,f,4,1,a.data(),4,p,b.data(),4,info),"RHS threshold"); untouched();
            setenv(PREFIX "_LAPACK_SOLVE_CUDA_MIN_RHS","0",1);
        }
        setenv(PREFIX "_LAPACK_CUDA","0",1);
        check(!gpu_solve(driver,chol,f,4,1,a.data(),4,p,b.data(),4,info),"disabled"); untouched(); unsetenv(PREFIX "_LAPACK_CUDA");
    }
}
void precision(long p) {
#if defined(MPLAPACK_BUILD_WITH_MPFR)
    mpfrxx::set_default_precision_bits(p); mpfrxx::set_default_mpc_precision_bits(p);
#else
    (void)p;
#endif
    reject<REAL>(); reject<COMPLEX>();
#if defined(MPLAPACK_BUILD_WITH_MPFR)
    std::vector<REAL> a(16,REAL(7)),b(16,REAL(11)); auto aa=a,bb=b;
    INTEGER piv[4]={1,2,3,4},info=123;
    auto rejected=[&] { check(!ENTRY(Rgesv)(4,1,a.data(),4,piv,b.data(),4,info),"MPFR eligibility"); check(info==123 && piv[0]==1 && piv[1]==2 && piv[2]==3 && piv[3]==4 && same(a,aa) && same(b,bb),"MPFR no publication"); };
    mpfr_set_default_rounding_mode(MPFR_RNDD); rejected(); mpfr_set_default_rounding_mode(MPFR_RNDN);
    mpfr_set_prec(a[0].get_mpfr_t(),p+1); a[0]=7; aa=a; rejected();
    check(mpfr_get_prec(a[0].get_mpfr_t())==p+1,"rejected precision unchanged");
    mpfr_set_prec(a[0].get_mpfr_t(),p); a[0]=7; aa=a;
    mpfr_set_inf(b[0].get_mpfr_t(),1); bb=b; rejected(); b[0]=11; bb=b;
    setenv("MPLAPACK_MPFR_CUDA","0",1); rejected(); unsetenv("MPLAPACK_MPFR_CUDA");
    // Cholesky eligibility matches factor entries even in the unused triangle.
    mpfr_set_prec(a[2].get_mpfr_t(),p+1); a[2]=7; aa=a;
    for(bool driver:{false,true}) {
        check(!gpu_solve(driver,true,"U",4,1,a.data(),4,piv,b.data(),4,info),"unused real triangle precision");
        check(info==123 && same(a,aa) && same(b,bb) && mpfr_get_prec(a[2].get_mpfr_t())==p+1,"unused real triangle untouched");
    }
    mpfr_set_prec(a[2].get_mpfr_t(),p); mpfr_set_inf(a[2].get_mpfr_t(),1); aa=a;
    for(bool driver:{false,true}) {
        check(!gpu_solve(driver,true,"U",4,1,a.data(),4,piv,b.data(),4,info),"unused real triangle finite");
        check(info==123 && same(a,aa) && same(b,bb),"unused real infinity untouched");
    }
    std::vector<COMPLEX> c(16,COMPLEX(7)),d(16,COMPLEX(11)); auto cc=c,dd=d;
    mpfr_set_prec(mpc_realref(c[2].mpc_data()),p+1); mpfr_set_ui(mpc_realref(c[2].mpc_data()),7,MPFR_RNDN); cc=c;
    for(bool driver:{false,true}) {
        check(!gpu_solve(driver,true,"U",4,1,c.data(),4,piv,d.data(),4,info),"unused complex triangle precision");
        check(info==123 && same(c,cc) && same(d,dd) && mpfr_get_prec(mpc_realref(c[2].mpc_data()))==p+1,"unused complex triangle untouched");
    }
    mpfr_set_prec(mpc_realref(c[2].mpc_data()),p); c[2]=COMPLEX(7);
    mpfr_set_inf(mpc_imagref(c[2].mpc_data()),1); cc=c;
    for(bool driver:{false,true}) {
        check(!gpu_solve(driver,true,"U",4,1,c.data(),4,piv,d.data(),4,info),"unused complex triangle finite");
        check(info==123 && same(c,cc) && same(d,dd),"unused complex infinity untouched");
    }
#endif
}
}
int main(int argc,char **argv) {
#if defined(MPLAPACK_SOLVER_DEVICE_TEST)
    int devices=0; if(cudaGetDeviceCount(&devices)!=cudaSuccess || !devices) return 77;
#endif
    setenv("MPLAPACK_LAPACK_FAST","0",1); // CPU oracle is the reference-order mode.
    setenv(PREFIX "_LAPACK_CUDA","1",1);
    setenv(PREFIX "_LAPACK_CUDA_MIN","0",1);
    setenv(PREFIX "_LAPACK_SOLVE_CUDA_MIN","0",1);
    setenv(PREFIX "_LAPACK_DRIVER_CUDA_MIN","0",1);
    setenv(PREFIX "_LAPACK_SOLVE_CUDA_MIN_RHS","0",1);
#if defined(MPLAPACK_BUILD_WITH_MPFR)
    setenv("MPLAPACK_MPFR_CUDA","1",1); setenv("MPLAPACK_MPFR_CUDA_RUNTIME","1",1);
    std::vector<long> precisions={256,512,768,1024,2048,200};
#else
    std::vector<long> precisions={0};
#endif
    if(argc>1) precisions={std::strtol(argv[1],nullptr,10)};
    bool quick=argc>2 && !std::strcmp(argv[2],"--quick");
    for(long p:precisions) {
        precision(p);
        std::vector<INTEGER> sizes=quick?std::vector<INTEGER>{1,3,17,65}:std::vector<INTEGER>{1,2,3,63,64,65,131,257,599};
        for(INTEGER n:sizes) for(bool driver:{false,true}) for(bool chol:{false,true}) {
            for(INTEGER rhs:{INTEGER(1),n}) {
                if(n>131 && rhs==n) continue;
                for(const char *f:chol?std::vector<const char*>{"U","L"}:driver?std::vector<const char*>{"N"}:std::vector<const char*>{"N","T","C"}) {
                    test<REAL>(n,rhs,driver,chol,f,0); test<COMPLEX>(n,rhs,driver,chol,f,0);
                    if(n<=65) for(int kind:{1,2,3,4}) {
                        if(kind==3 && !driver) continue;
                        test<REAL>(n,rhs,driver,chol,f,kind); test<COMPLEX>(n,rhs,driver,chol,f,kind);
                    }
                }
            }
        }
        std::printf("PASS p=%ld cumulative checks=%ld elements=%ld zero-sign-errors=%ld\n",p,checks,comparisons,zero_sign_errors); std::fflush(stdout);
    }
}
