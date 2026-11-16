// Wavefront publication-order negative control (factorization goal 13).
// The flag-driven trisolve wavefront orders its own publication: the owner
// stores the divided x(j) first and publishes the per-element flag with a
// release afterwards; dependents acquire the flag before reading x(j).  This
// control breaks the publication with the MPLAPACK_TEST_WAVE_PUBLISH_ONLY
// hook (the flag is raised, the divided value is never stored — the "mark
// first, value never follows" shape) and asserts that the result is WRONG,
// i.e. the flag protocol is load-bearing.  Without the hook the same solve
// must stay bit-identical to the reference.  Exit 0 iff both hold.
#include <mpblas.h>
#include "../Rblas_cuda_entry_common.h"
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <type_traits>
#define Rtrsv ref_Rtrsv
#include "../../../reference/Rtrsv.cpp"
#undef Rtrsv
#ifdef MPLAPACK_GEMM_CUDA_HOST_EMULATION
#if defined(MPLAPACK_BUILD_WITH_DD)
#include "../../dd/cuda/Rgemm_gpu_cuda.cu"
#elif defined(MPLAPACK_BUILD_WITH_TD)
#include "../../td/cuda/Rgemm_gpu_cuda.cu"
#else
#include "../../qd/cuda/Rgemm_gpu_cuda.cu"
#endif
#elif !defined(MPLAPACK_MPFR_CUDA_HOST_EMULATION)
#include <cuda_runtime.h>
#endif
#if defined(MPLAPACK_BUILD_WITH_MPFR)
#define PREFIX "MPLAPACK_MPFR"
#elif defined(MPLAPACK_BUILD_WITH_DD)
#define PREFIX "MPLAPACK_DD"
#elif defined(MPLAPACK_BUILD_WITH_TD)
#define PREFIX "MPLAPACK_TD"
#else
#define PREFIX "MPLAPACK_QD"
#endif
namespace {
long checks=0,failures=0;
void check(bool ok,const char *name) { ++checks;if(!ok) { std::printf("FAIL %s\n",name);++failures; } }
bool same(const REAL &a,const REAL &b) {
#ifdef MPLAPACK_BUILD_WITH_MPFR
    return mpfr_equal_p(a.get_mpfr_t(),b.get_mpfr_t()) && mpfr_signbit(a.get_mpfr_t())==mpfr_signbit(b.get_mpfr_t());
#else
    return std::memcmp(&a,&b,sizeof a)==0;
#endif
}
REAL number(long i,long j) { return REAL((i*17+j*13)%19-9)/REAL(997.0); }
}
int main() {
#if !defined(MPLAPACK_GEMM_CUDA_HOST_EMULATION) && !defined(MPLAPACK_MPFR_CUDA_HOST_EMULATION)
    int devices=0;if(cudaGetDeviceCount(&devices) || !devices) return 77;
#endif
#ifdef MPLAPACK_BUILD_WITH_MPFR
    // the wavefront runs for the checked precisions only
    if(mpfrxx::default_precision_bits()>4096) return 77;
#endif
    setenv(PREFIX "_BLAS_CUDA_MIN","0",1);setenv(PREFIX "_BLAS_CUDA_L12","1",1);
    setenv(PREFIX "_CUDA_MIN_MNK","0",1);setenv(PREFIX "_CUDA_L12","1",1);
    setenv("MPLAPACK_MPFR_CUDA_MIN_MNK","0",1);setenv("MPLAPACK_MPFR_CUDA_L12","1",1);
    // upper, non-unit, NoTrans: every element consumes its own division, so a
    // missing publication is visible in every result column.
    const long n=65,lda=n+2;
    std::vector<REAL> a(lda*n,REAL(0.0)),x(n,REAL(0.0));
    for(long j=0;j<n;++j) for(long i=0;i<=j;++i) a[i+j*lda]=i==j?REAL(1.3)+number(i,j):number(i,j);
    for(long i=0;i<n;++i) x[i]=i%7?number(i,2):REAL(-0.0);
    const std::vector<REAL> a0=a,x0=x;
    std::vector<REAL> expected=x;
    ref_Rtrsv("U","N","N",n,a.data(),lda,expected.data(),1);
    // 1. the ordered run is bit-identical to the reference
    x=x0;
    if(!MPLAPACK_CUDA_ENTRY(Rtrsv)("U","N","N",n,a.data(),lda,x.data(),1)) { std::printf("FAIL ordered run not computed\n");return 1; }
    bool correct=true;for(long i=0;i<n;++i) if(!same(x[i],expected[i])) correct=false;
    check(correct,"ordered run matches the reference");
    // 2. the broken publication (flag raised, value never stored) is WRONG
    setenv("MPLAPACK_TEST_WAVE_PUBLISH_ONLY","1",1);
    x=x0;
    if(!MPLAPACK_CUDA_ENTRY(Rtrsv)("U","N","N",n,a.data(),lda,x.data(),1)) { std::printf("FAIL broken run not computed\n");return 1; }
    unsetenv("MPLAPACK_TEST_WAVE_PUBLISH_ONLY");
    bool wrong=false;for(long i=0;i<n;++i) if(!same(x[i],expected[i])) wrong=true;
    check(wrong,"broken publication detected (the result differs)");
    bool amat=true;for(size_t i=0;i<a.size();++i) if(!same(a[i],a0[i])) amat=false;
    check(amat,"matrix unchanged");
    std::printf("wave order control: %ld checks, %ld failures\n",checks,failures);
    return failures?1:0;
}
