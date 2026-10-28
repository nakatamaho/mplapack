// Full entry-point contract, including untouched opposite triangle and padding.
#if defined(__CUDACC__)
#define MPLAPACK_MINMAX_MPLAPACKINT_DEFINED
#endif
#include <mpblas.h>
#include <mplapack.h>
#include <vector>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <limits>
#ifdef _OPENMP
#include <omp.h>
#endif
#if defined(MPLAPACK_BUILD_WITH_DD)
#include <mpblas_dd_cuda.h>
#define ENTRY(R) R##_dd_cuda
#define PREFIX "MPLAPACK_DD"
#elif defined(MPLAPACK_BUILD_WITH_TD)
#include <mpblas_td_cuda.h>
#define ENTRY(R) R##_td_cuda
#define PREFIX "MPLAPACK_TD"
#else
#include <mpblas_qd_cuda.h>
#define ENTRY(R) R##_qd_cuda
#define PREFIX "MPLAPACK_QD"
#endif
#ifdef MPLAPACK_GEMM_CUDA_HOST_EMULATION
#include "../Rpotrf_cuda_entry.cuh"
#endif
#ifdef __CUDACC__
#include <cuda_runtime.h>
#endif
void Rpotrf_reference_cpu(const char *,INTEGER,REAL *,INTEGER,INTEGER &);
void Cpotrf_reference_cpu(const char *,INTEGER,COMPLEX *,INTEGER,INTEGER &);
namespace {
long checks=0;
void check(bool ok,const char *what) {
    ++checks;
    if(!ok) { std::fprintf(stderr,"FAIL %s at %ld\n",what,checks); std::exit(1); }
}
bool entry(const char *u,INTEGER n,REAL *a,INTEGER ld,INTEGER &info) { return ENTRY(Rpotrf)(u,n,a,ld,info); }
bool entry(const char *u,INTEGER n,COMPLEX *a,INTEGER ld,INTEGER &info) { return ENTRY(Cpotrf)(u,n,a,ld,info); }
void reference(const char *u,INTEGER n,REAL *a,INTEGER ld,INTEGER &info) { Rpotrf_reference_cpu(u,n,a,ld,info); }
void reference(const char *u,INTEGER n,COMPLEX *a,INTEGER ld,INTEGER &info) { Cpotrf_reference_cpu(u,n,a,ld,info); }
REAL diagonal(long n,REAL *) { return REAL(n+1); }
COMPLEX diagonal(long n,COMPLEX *) { return COMPLEX(n+1,0.125); }
REAL element(long i,long j,REAL *) { return REAL(((i*17+j*13)%31)-15)/REAL(128); }
COMPLEX element(long i,long j,COMPLEX *) { return COMPLEX(element(i,j,(REAL*)nullptr),REAL(i-j)/REAL(512)); }
REAL conjugate(const REAL &x) { return x; }
COMPLEX conjugate(const COMPLEX &x) { return conj(x); }
template<class H> void residual(const char *u,INTEGER n,INTEGER ld,const std::vector<H> &a) {
    REAL norm=0,error=0;
    for(long j=0;j<n;++j) {
        REAL col=0,err=0;
        for(long i=0;i<n;++i) {
            H original=i==j?H(n+1):element(std::min(i,j),std::max(i,j),(H*)nullptr);
            if((*u=='U' && i>j) || (*u=='L' && i<j)) original=conjugate(original);
            H sum=0;
            for(long k=0;k<=std::min(i,j);++k)
                sum+=*u=='U'?conjugate(a[k+i*ld])*a[k+j*ld]:a[i+k*ld]*conjugate(a[j+k*ld]);
            col+=abs(original); err+=abs(original-sum);
        }
        norm=std::max(norm,col); error=std::max(error,err);
    }
    check(error/(REAL(n)*norm*Rlamch("E"))<=REAL(30),"backward error <=30");
}
template<class H> void test(const char *u,INTEGER n,INTEGER bad) {
    INTEGER ld=n+3;
    std::vector<H> a(ld*n,H(19)),b;
    for(long j=0;j<n;++j) for(long i=0;i<n;++i) {
        long lo=std::min(i,j),hi=std::max(i,j);
        a[i+j*ld]=i==j?diagonal(n,(H*)nullptr):element(lo,hi,(H*)nullptr);
    }
    // The diagonal can carry an imaginary part; successful sqrt must clear it.
    if(bad) a[(bad-1)+(bad-1)*ld]=H(-1);
    b=a;
    INTEGER ir=-91,ig=-92;
    reference(u,n,a.data(),ld,ir);
    check(entry(u,n,b.data(),ld,ig),"computed");
    if(ir!=ig || std::memcmp(a.data(),b.data(),a.size()*sizeof(H))) {
        std::fprintf(stderr,"%s n=%ld bad=%ld info=%ld/%ld\n",u,(long)n,(long)bad,(long)ir,(long)ig);
        for(size_t i=0;i<a.size();++i) if(std::memcmp(&a[i],&b[i],sizeof(H))) { std::fprintf(stderr,"element=%zu\n",i); break; }
        check(false,"bit identity");
    }
    check(true,"bit identity");
    if(!bad && n<=131) residual(u,n,ld,b);
}
template<class H> void reject() {
    std::vector<H> a(16,H(7)),old=a;
    INTEGER info=123;
    for(auto u:{"X","U"}) for(INTEGER n:{INTEGER(-1),INTEGER(0),INTEGER(4)}) {
        if(*u=='U' && n==4) continue;
        check(!entry(u,n,a.data(),4,info),"reject invalid/quick");
        check(info==123 && !std::memcmp(a.data(),old.data(),16*sizeof(H)),"untouched");
    }
    check(!entry("L",4,a.data(),3,info),"reject lda");
    check(!entry(nullptr,4,a.data(),4,info),"reject null flag");
    check(!entry("U",4,(H*)nullptr,4,info),"reject null matrix");
    INTEGER huge=std::numeric_limits<INTEGER>::max();
    check(!entry("U",huge,a.data(),huge,info),"reject size overflow");
    setenv(PREFIX "_LAPACK_CUDA","0",1);
    check(!entry("U",4,a.data(),4,info),"disabled");
    unsetenv(PREFIX "_LAPACK_CUDA");
    setenv(PREFIX "_LAPACK_CUDA_MIN","5",1);
    check(!entry("U",4,a.data(),4,info),"threshold");
    check(info==123 && !std::memcmp(a.data(),old.data(),16*sizeof(H)),"untouched");
    setenv(PREFIX "_LAPACK_CUDA_MIN","0",1);
}
}
int main() {
#ifdef _OPENMP
    omp_set_num_threads(1);
#endif
#ifdef __CUDACC__
    int devices=0;
    if(cudaGetDeviceCount(&devices)!=cudaSuccess || !devices) return 77;
#endif
    setenv(PREFIX "_LAPACK_CUDA_MIN","0",1);
    reject<REAL>(); reject<COMPLEX>();
    INTEGER nb=iMlaenv(1,"Rpotrf","L",600,-1,-1,-1);
    for(INTEGER n:{INTEGER(1),INTEGER(2),INTEGER(3),nb-1,nb,nb+1,2*nb+3,INTEGER(131),INTEGER(257),INTEGER(599)})
        for(auto u:{"U","L"}) for(INTEGER bad:{INTEGER(0),INTEGER(1),n/2+1,n}) {
            test<REAL>(u,n,bad); test<COMPLEX>(u,n,bad);
        }
    setenv(PREFIX "_LAPACK_FAST","1",1);
    for(auto u:{"U","L"}) { test<REAL>(u,nb+1,0); test<COMPLEX>(u,nb+1,0); }
    unsetenv(PREFIX "_LAPACK_FAST");
    std::printf("PASS %ld Cholesky entry checks\n",checks);
}
