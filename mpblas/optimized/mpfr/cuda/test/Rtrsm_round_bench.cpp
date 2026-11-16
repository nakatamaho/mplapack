// End-to-end GPU entry-point timing. Select the old engine in a separate
// process with MPLAPACK_MPFR_CUDA_FORCE_RUNTIME=1 (the setting is cached).
// Alternate old/new processes five times; report median, never best-of.
#include <mpblas_mpfr_cuda.h>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <vector>
int main(int argc,char** argv)
{
    if(argc!=6) { std::fprintf(stderr,"usage: %s PB n side uplo trans\n",argv[0]); return 1; }
    int p=std::atoi(argv[1]); mplapackint n=std::atol(argv[2]);
    if(n<1 || (p!=256 && p!=512 && p!=1024 && p!=2048)) return 1;
    mpfrxx::set_default_precision_bits(p);
    setenv("MPLAPACK_MPFR_CUDA_MIN_MNK","0",1);
    std::vector<mpfr_class> a(n*n),b(n*n),c(n*n);
    unsigned long long seed=7;
    for(mplapackint j=0;j<n;++j) for(mplapackint i=0;i<n;++i) {
        seed=seed*6364136223846793005ULL+1;
        a[i+j*n]=i==j?2.0:((double)((seed>>32)%1024)-512.0)/(1024.0*n);
        seed=seed*6364136223846793005ULL+1;
        b[i+j*n]=((double)((seed>>32)%1024)-512.0)/1024.0;
        // Fill all PB bits, not just double-sized significands. Preserve the
        // small off-diagonal exponents so the triangular system stays tame.
        for(mpfr_class* v : {&a[i+j*n], &b[i+j*n]}) {
            if(v==&a[i+j*n] && i==j) continue;
            mpfr_ptr x=v->get_mpfr_t();
            if(mpfr_zero_p(x)) continue;
            auto limbs=static_cast<mp_limb_t*>(mpfr_custom_get_significand(x));
            const int count=(p+63)/64;
            for(int k=0;k<count;++k) { seed=seed*6364136223846793005ULL+1; limbs[k]=seed; }
            limbs[count-1]|=mp_limb_t(1)<<63;
        }
    }
    c=b; mpfr_class alpha=-0.5;
    if(!Rtrsm_mpfr_cuda(argv[3],argv[4],argv[5],"N",n,n,alpha,a.data(),n,c.data(),n)) return 2;
    c=b;
    auto start=std::chrono::steady_clock::now();
    bool ok=Rtrsm_mpfr_cuda(argv[3],argv[4],argv[5],"N",n,n,alpha,a.data(),n,c.data(),n);
    double t=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    const char* runtime=std::getenv("MPLAPACK_MPFR_CUDA_FORCE_RUNTIME");
    std::printf("PB=%d n=%ld %s%s%s runtime=%s seconds=%.9f ok=%d sample=%.17g\n",
        p,(long)n,argv[3],argv[4],argv[5],runtime?runtime:"0",t,ok,mpfr_get_d(c[0].get_mpfr_t(),MPFR_RNDN));
    return ok?0:2;
}
