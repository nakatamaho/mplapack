// Direct LU entry contract: whole padded A, pivots and info versus CPU oracle.
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
#include "../Rgetrf_cuda_entry.cuh"
#endif
#ifdef __CUDACC__
#include <cuda_runtime.h>
#pragma GCC visibility push(hidden)
#include "../Rgetrf_gpu_common.cuh"
#pragma GCC visibility pop
#endif
void Rgetrf2_reference_cpu(INTEGER,INTEGER,REAL *,INTEGER,INTEGER *,INTEGER &);
void Cgetrf2_reference_cpu(INTEGER,INTEGER,COMPLEX *,INTEGER,INTEGER *,INTEGER &);
void Rgetrf_reference_cpu(INTEGER,INTEGER,REAL *,INTEGER,INTEGER *,INTEGER &);
void Cgetrf_reference_cpu(INTEGER,INTEGER,COMPLEX *,INTEGER,INTEGER *,INTEGER &);
namespace {
long checks=0;
void check(bool ok,const char *what) {
    ++checks; if(!ok) { std::fprintf(stderr,"FAIL %s at %ld\n",what,checks); std::exit(1); }
}
bool entry(INTEGER m,INTEGER n,REAL *a,INTEGER ld,INTEGER *p,INTEGER &info) { return ENTRY(Rgetrf)(m,n,a,ld,p,info); }
bool entry(INTEGER m,INTEGER n,COMPLEX *a,INTEGER ld,INTEGER *p,INTEGER &info) { return ENTRY(Cgetrf)(m,n,a,ld,p,info); }
void reference(INTEGER m,INTEGER n,REAL *a,INTEGER ld,INTEGER *p,INTEGER &info) { Rgetrf_reference_cpu(m,n,a,ld,p,info); }
void reference(INTEGER m,INTEGER n,COMPLEX *a,INTEGER ld,INTEGER *p,INTEGER &info) { Cgetrf_reference_cpu(m,n,a,ld,p,info); }
REAL element(long i,long j,REAL *) { return REAL(((i*17+j*13+i*j*7)%127)-63)/REAL(128); }
COMPLEX element(long i,long j,COMPLEX *) { return COMPLEX(element(i,j,(REAL*)nullptr),element(j+5,i+3,(REAL*)nullptr)); }
template<class H> void residual(INTEGER m,INTEGER n,INTEGER ld,const std::vector<H> &original,const std::vector<H> &a,const std::vector<INTEGER> &p) {
    auto b=original; long mn=std::min(m,n);
    for(long k=0;k<mn;++k) for(long j=0;j<n;++j) std::swap(b[k+j*ld],b[p[k]-1+j*ld]);
    REAL norm=0,error=0;
    for(long j=0;j<n;++j) {
        REAL col=0,err=0;
        for(long i=0;i<m;++i) {
            H sum=0;
            for(long k=0;k<mn && k<=i && k<=j;++k) sum+=(i==k?H(1):a[i+k*ld])*a[k+j*ld];
            col+=abs(original[i+j*ld]); err+=abs(b[i+j*ld]-sum);
        }
        norm=std::max(norm,col); error=std::max(error,err);
    }
    check(norm==REAL(0) || error/(REAL(std::max(m,n))*norm*Rlamch("E"))<=REAL(30),"backward error <=30");
}
template<class H> void test(INTEGER m,INTEGER n,int kind) {
    INTEGER ld=m+3,mn=std::min(m,n);
    std::vector<H> a((size_t)ld*n,H(19)),b;
    for(long j=0;j<n;++j) for(long i=0;i<m;++i) {
        a[i+j*ld]=element(i,j,(H*)nullptr);
        if(kind==1 && j==mn/2) a[i+j*ld]=H(0); // zero column
        if(kind==2 && j==n-1 && n>1) a[i+j*ld]=a[i]; // exact duplicate
        if(kind==3) a[i+j*ld]=H(i==m-1-j?4:(i==j?1:0)); // repeated swaps
        if(kind==4) a[i+j*ld]=H(i>=j?1:-1); // tied maxima
    }
    b=a; auto original=a;
    std::vector<INTEGER> p(mn+3,-8),q=p;
    INTEGER ir=-91,ig=-92;
    reference(m,n,a.data(),ld,p.data(),ir);
    check(entry(m,n,b.data(),ld,q.data(),ig),"computed");
    if(ir!=ig || p!=q || std::memcmp(a.data(),b.data(),a.size()*sizeof(H))) {
        std::fprintf(stderr,"m=%ld n=%ld kind=%d info=%ld/%ld\n",(long)m,(long)n,kind,(long)ir,(long)ig);
        for(size_t i=0;i<a.size();++i) if(std::memcmp(&a[i],&b[i],sizeof(H))) { std::fprintf(stderr,"element=%zu\n",i); break; }
        check(false,"A/ipiv/info bit identity");
    }
    check(true,"A/ipiv/info bit identity");
    if(std::max(m,n)<=131) residual(m,n,ld,original,b,q);
}
#if defined(MPLAPACK_BUILD_WITH_DD)
using TestOps=mplapack_gemm::dd_ops;
#elif defined(MPLAPACK_BUILD_WITH_TD)
using TestOps=mplapack_gemm::td_ops;
#else
using TestOps=mplapack_gemm::qd_ops;
#endif
using TestR=mplapack_blas_dev::dreal<TestOps>;
using TestC=mplapack_blas_dev::dcomplex<TestR>;
TestR pack(const REAL &x) { TestR r; for(size_t i=0;i<sizeof(r.v.x)/sizeof(double);++i) r.v.x[i]=x.x[i]; return r; }
TestC pack(const COMPLEX &x) { return TestC(pack(x.real()),pack(x.imag())); }
void unpack(const TestR &r,REAL &x) { for(size_t i=0;i<sizeof(r.v.x)/sizeof(double);++i) x.x[i]=r.v.x[i]; }
void unpack(const TestC &r,COMPLEX &x) { REAL re,im; unpack(r.real(),re); unpack(r.imag(),im); x=COMPLEX(re,im); }
void recursive_reference(INTEGER m,INTEGER n,REAL *a,INTEGER ld,INTEGER *p,INTEGER &i) { Rgetrf2_reference_cpu(m,n,a,ld,p,i); }
void recursive_reference(INTEGER m,INTEGER n,COMPLEX *a,INTEGER ld,INTEGER *p,INTEGER &i) { Cgetrf2_reference_cpu(m,n,a,ld,p,i); }
template<class Blas,class H> void unblocked() {
    // Force nb=1 with only O(max(m,n)) scratch, including a wide matrix.
    using namespace mplapack_lapack_gpu;
    using E=typename Blas::E;
    const INTEGER m=17,n=23,ld=20;
    std::vector<H> a(ld*n,H(19)),b;
    std::vector<E> staging(m*n);
    for(long j=0;j<n;++j) for(long i=0;i<m;++i) { a[i+j*ld]=element(i,j,(H*)nullptr); staging[i+j*m]=pack(a[i+j*ld]); }
    b=a; INTEGER ir=-1; std::vector<INTEGER> p(m),q(m);
    recursive_reference(m,n,a.data(),ld,p.data(),ir);
    matrix<E> device; buffer<long> piv;
    lu_context<Blas,TestR> ctx(n,1,pack(Rlamch("S")));
    check(!ctx.error && !device.allocate(m,n) && !piv.allocate(m) && !device.copy(staging.data(),m,false),"unblocked allocation/upload");
#ifdef __CUDACC__
    check(cudaStreamSynchronize(nullptr)==cudaSuccess,"unblocked upload fence");
#endif
    long ig=lu(ctx,m,n,1,device.data(),m,piv.data());
    std::vector<long> dp(m);
    check(!ctx.error && !device.copy(staging.data(),m,true) && !transfer(dp.data(),piv.data(),m,true),"unblocked download");
    for(long j=0;j<n;++j) for(long i=0;i<m;++i) unpack(staging[i+j*m],b[i+j*ld]);
    for(long i=0;i<m;++i) q[i]=dp[i];
    check(ir==ig && p==q && !std::memcmp(a.data(),b.data(),a.size()*sizeof(H)),"unblocked bit identity / bounded scratch");
}
template<class H> void tiny() {
    // Both reciprocal-scaling and elementwise division branches, plus zero.
    for(int kind=0;kind<3;++kind) {
        INTEGER m=4,n=1,ld=7,ir=-1,ig=-2;
        REAL s=kind==0?Rlamch("S")/REAL(4):Rlamch("S")*REAL(4);
        std::vector<H> a(ld,H(19));
        for(long i=0;i<m;++i) a[i]=kind==2?H(0):H(s/REAL(i+1));
        auto b=a; INTEGER p[4]={-1,-1,-1,-1},q[4]={-1,-1,-1,-1};
        reference(m,n,a.data(),ld,p,ir);
        check(entry(m,n,b.data(),ld,q,ig),"tiny computed");
        check(ir==ig && !std::memcmp(p,q,sizeof(p)) && !std::memcmp(a.data(),b.data(),a.size()*sizeof(H)),"tiny A/ipiv/info bit identity");
    }
}
template<class H> void reject() {
    std::vector<H> a(16,H(7)),old=a;
    INTEGER p[4]={9,8,7,6},info=123;
    auto untouched=[&]() { check(info==123 && p[0]==9 && p[1]==8 && p[2]==7 && p[3]==6 && !std::memcmp(a.data(),old.data(),16*sizeof(H)),"untouched"); };
    for(INTEGER m:{INTEGER(-1),INTEGER(0),INTEGER(4)}) for(INTEGER n:{INTEGER(-1),INTEGER(0),INTEGER(4)}) {
        if(m==4 && n==4) continue;
        check(!entry(m,n,a.data(),4,p,info),"invalid/quick"); untouched();
    }
    check(!entry(4,4,a.data(),3,p,info),"lda");
    check(!entry(4,4,(H*)nullptr,4,p,info),"null matrix");
    check(!entry(4,4,a.data(),4,nullptr,info),"null pivots");
    INTEGER huge=std::numeric_limits<INTEGER>::max();
    check(!entry(huge,huge,a.data(),huge,p,info),"overflow");
    setenv(PREFIX "_LAPACK_CUDA","0",1);
    check(!entry(4,4,a.data(),4,p,info),"disabled"); unsetenv(PREFIX "_LAPACK_CUDA");
    setenv(PREFIX "_LAPACK_CUDA_MIN","5",1);
    check(!entry(4,4,a.data(),4,p,info),"threshold"); untouched();
    setenv(PREFIX "_LAPACK_CUDA_MIN","0",1);
}
}
int main() {
#ifdef __CUDACC__
    int devices=0; if(cudaGetDeviceCount(&devices)!=cudaSuccess || !devices) return 77;
#endif
    setenv(PREFIX "_LAPACK_CUDA_MIN","0",1);
    reject<REAL>(); reject<COMPLEX>();
    tiny<REAL>(); tiny<COMPLEX>();
    unblocked<mplapack_lapack_gpu::real_lu<TestOps>,REAL>();
    unblocked<mplapack_lapack_gpu::complex_lu<TestOps>,COMPLEX>();
    INTEGER nb=iMlaenv(1,"Rgetrf"," ",600,600,-1,-1);
    for(INTEGER n:{INTEGER(1),INTEGER(2),INTEGER(3),nb-1,nb,nb+1,2*nb+3,INTEGER(131),INTEGER(257),INTEGER(599)})
        for(INTEGER m:{n,std::max(INTEGER(1),n-3),n+5}) for(int kind=0;kind<5;++kind) {
            test<REAL>(m,n,kind); test<COMPLEX>(m,n,kind);
        }
    setenv(PREFIX "_LAPACK_FAST","1",1);
    test<REAL>(nb+1,nb+1,0); test<COMPLEX>(nb+1,nb+1,0);
    unsetenv(PREFIX "_LAPACK_FAST");
    std::printf("PASS %ld LU entry checks\n",checks);
}
