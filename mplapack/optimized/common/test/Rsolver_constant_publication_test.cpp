// BF drivers must wait for context/scalar uploads before nonblocking kernels.
#include <mpblas.h>
#include <cuda_runtime_api.h>
#include <cstdio>
#include <cstdlib>
#include <dlfcn.h>
#include <vector>
#if defined(MPLAPACK_BUILD_WITH_DD)
#include <mpblas_dd_cuda.h>
#define ENTRY Rposv_dd_cuda
#define PREFIX "MPLAPACK_DD"
#elif defined(MPLAPACK_BUILD_WITH_TD)
#include <mpblas_td_cuda.h>
#define ENTRY Rposv_td_cuda
#define PREFIX "MPLAPACK_TD"
#else
#include <mpblas_qd_cuda.h>
#define ENTRY Rposv_qd_cuda
#define PREFIX "MPLAPACK_QD"
#endif
namespace {
int count(const char *name) { auto f=reinterpret_cast<int(*)()>(dlsym(RTLD_DEFAULT,name)); return f?f():-1; }
}
int main(int argc,char **argv) {
    int devices=0; if(cudaGetDeviceCount(&devices)!=cudaSuccess || !devices) return 77;
    if(argc!=2) return 2; bool control=std::atoi(argv[1])!=0;
    setenv(PREFIX "_LAPACK_CUDA","1",1); setenv(PREFIX "_LAPACK_CUDA_MIN","0",1);
    setenv(PREFIX "_LAPACK_DRIVER_CUDA_MIN","0",1);
    constexpr INTEGER n=3,lda=5,ldb=6;
    const int values[3][3]={{4,2,2},{2,5,3},{2,3,6}};
    std::vector<REAL> a(lda*n,REAL(19)),b(ldb,REAL(23));
    for(long j=0;j<n;++j) for(long i=0;i<n;++i) a[i+j*lda]=REAL(values[i][j]);
    b[0]=14; b[1]=21; b[2]=26;
    setenv("MPLAPACK_TEST_DELAY_BF_CONSTANTS","0",1);
    { auto wa=a,wb=b; INTEGER info=-1; if(!ENTRY("L",n,1,wa.data(),lda,wb.data(),ldb,info) || info) return control?0:1; }
    auto expected=a,expected_b=b;
    for(long j=0;j<n;++j) for(long i=j;i<n;++i) expected[i+j*lda]=REAL(i==j?2:1);
    for(long i=0;i<n;++i) expected_b[i]=REAL(i+1);
    setenv("MPLAPACK_TEST_DELAY_BF_CONSTANTS","1",1);
    INTEGER info=-987; bool computed=ENTRY("L",n,1,a.data(),lda,b.data(),ldb,info);
    long errors=0; for(size_t i=0;i<a.size();++i) if(a[i]!=expected[i]) ++errors;
    for(size_t i=0;i<b.size();++i) if(b[i]!=expected_b[i]) ++errors;
    int copies=count("mplapack_test_pivot_injections"),waits=count("mplapack_test_pivot_waits"),skips=count("mplapack_test_pivot_skips");
    std::printf("%s computed=%d info=%ld errors=%ld copies=%d fences=%d skipped=%d control=%d\n",PREFIX,computed?1:0,long(info),errors,copies,waits,skips,control?1:0);
    auto release=reinterpret_cast<void(*)()>(dlsym(RTLD_DEFAULT,"mplapack_test_pivot_release")); if(release) release();
    bool correct=computed && info==0 && errors==0;
    if(control) return copies==1 && skips==1 && !correct?1:0; // CTest WILL_FAIL.
    return correct && copies==1 && waits==1 && skips==0?0:1;
}
