// End-to-end real factorization benchmark; no fallback in "entry" mode.
// Build against a CPU library for "cpu", the retained baseline CUDA library
// for "old", or the new CUDA library for "entry". One fresh process per call.
#include <mpblas.h>
#include <mplapack.h>
#include <mpblas_mpfr_cuda.h>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
int main(int argc,char **argv) {
    if(argc!=5) { std::fprintf(stderr,"usage: %s precision n getrf|potrf cpu|old|entry\n",argv[0]); return 2; }
    long prec=std::atol(argv[1]),n=std::atol(argv[2]);
    if(n<1 || prec<2) return 2;
    mpfrxx::set_default_precision_bits(prec);
    if(std::strcmp(argv[3],"getrf") && std::strcmp(argv[3],"potrf")) return 2;
#ifdef MPLAPACK_FACTOR_BENCH_CPU_ONLY
    if(std::strcmp(argv[4],"cpu") && std::strcmp(argv[4],"old")) return 2;
#else
    if(std::strcmp(argv[4],"entry")) return 2;
#endif
    bool chol=!std::strcmp(argv[3],"potrf"),entry=!std::strcmp(argv[4],"entry");
    std::vector<REAL> a(n*n); std::vector<INTEGER> piv(n); INTEGER info=-1;
    gmp_randstate_t rng; gmp_randinit_default(rng); gmp_randseed_ui(rng,42);
    for(long j=0;j<n;++j) for(long i=j;i<n;++i) {
        mpfr_urandomb(a[i+j*n].get_mpfr_t(),rng);
        a[i+j*n]=(a[i+j*n]-REAL(0.5))/REAL(8);
        if(i==j) a[i+j*n]+=REAL(n+1);
        else a[j+i*n]=a[i+j*n];
    }
    gmp_randclear(rng);
    setenv("MPLAPACK_MPFR_LAPACK_CUDA_MIN","0",1);
    setenv("MPLAPACK_MPFR_CUDA_VERBOSE","1",1);
    if(!std::strcmp(argv[4],"old")) setenv("MPLAPACK_MPFR_LAPACK_CUDA","0",1);
    auto start=std::chrono::steady_clock::now(); bool ran=true;
#ifdef MPLAPACK_FACTOR_BENCH_CPU_ONLY
    (void)entry;
    if(chol) Rpotrf("L",n,a.data(),n,info); else Rgetrf(n,n,a.data(),n,piv.data(),info);
#else
    if(entry) {
        ran=chol?Rpotrf_mpfr_cuda("L",n,a.data(),n,info):Rgetrf_mpfr_cuda(n,n,a.data(),n,piv.data(),info);
    } else {
        if(chol) Rpotrf("L",n,a.data(),n,info); else Rgetrf(n,n,a.data(),n,piv.data(),info);
    }
#endif
    double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    std::printf("%ld,%ld,%s,%s,%.9f,%d,%ld\n",prec,n,argv[3],argv[4],seconds,ran?1:0,(long)info);
    return ran && !info?0:1;
}
