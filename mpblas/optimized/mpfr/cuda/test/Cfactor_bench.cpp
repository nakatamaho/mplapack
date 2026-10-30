// End-to-end complex factor benchmark. One fresh process per timed call.
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
    mpfrxx::set_default_precision_bits(prec); mpfrxx::set_default_mpc_precision_bits(prec);
    if(std::strcmp(argv[3],"getrf") && std::strcmp(argv[3],"potrf")) return 2;
#ifdef MPLAPACK_FACTOR_BENCH_CPU_ONLY
    if(std::strcmp(argv[4],"cpu") && std::strcmp(argv[4],"old")) return 2;
#else
    if(std::strcmp(argv[4],"entry")) return 2;
#endif
    bool chol=!std::strcmp(argv[3],"potrf");
    std::vector<COMPLEX> a(n*n); std::vector<INTEGER> piv(n); INTEGER info=-1;
    gmp_randstate_t rng; gmp_randinit_default(rng); gmp_randseed_ui(rng,42);
    for(long j=0;j<n;++j) for(long i=j;i<n;++i) {
        auto z=a[i+j*n].mpc_data();
        mpfr_urandomb(mpc_realref(z),rng); mpfr_urandomb(mpc_imagref(z),rng);
        a[i+j*n]=(a[i+j*n]-COMPLEX(.5,.5))/COMPLEX(8);
        if(i==j) {
            a[i+j*n]+=COMPLEX(n+1);
            mpfr_set_zero(mpc_imagref(z),1);
        } else a[j+i*n]=conj(a[i+j*n]);
    }
    gmp_randclear(rng);
    setenv("MPLAPACK_MPFR_LAPACK_CUDA_MIN","0",1);
    setenv("MPLAPACK_MPFR_CUDA_VERBOSE","1",1);
    if(!std::strcmp(argv[4],"old")) setenv("MPLAPACK_MPFR_LAPACK_CUDA","0",1);
    auto start=std::chrono::steady_clock::now(); bool ran=true;
#ifdef MPLAPACK_FACTOR_BENCH_CPU_ONLY
    if(chol) Cpotrf("L",n,a.data(),n,info); else Cgetrf(n,n,a.data(),n,piv.data(),info);
#else
    ran=chol?Cpotrf_mpfr_cuda("L",n,a.data(),n,info):Cgetrf_mpfr_cuda(n,n,a.data(),n,piv.data(),info);
#endif
    double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    std::printf("%ld,%ld,%s,%s,%.9f,%d,%ld\n",prec,n,argv[3],argv[4],seconds,ran?1:0,(long)info);
    return ran && !info?0:1;
}
