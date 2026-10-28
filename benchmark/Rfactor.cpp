// Reproducible factorization timings; intentionally separate from the legacy
// drivers and their benchmark/results output contracts.
#include <mpblas.h>
#include <mplapack.h>
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>
#ifdef MPLAPACK_FACTOR_CUDA
#include <cuda_runtime_api.h>
#endif

namespace {
template<class T> T scalar(std::mt19937 &rng) { return T(double(int(rng()%65536)-32768)/65536.0); }
template<> COMPLEX scalar<COMPLEX>(std::mt19937 &rng) {
    REAL r=scalar<REAL>(rng), i=scalar<REAL>(rng); return COMPLEX(r,i);
}
REAL adjoint(REAL x) { return x; }
COMPLEX adjoint(COMPLEX x) { return conj(x); }
void factor(bool pot, bool upper, INTEGER n, REAL *a, INTEGER *p, INTEGER &info) {
    if (pot) Rpotrf(upper?"U":"L",n,a,n,info); else Rgetrf(n,n,a,n,p,info);
}
void factor(bool pot, bool upper, INTEGER n, COMPLEX *a, INTEGER *p, INTEGER &info) {
    if (pot) Cpotrf(upper?"U":"L",n,a,n,info); else Cgetrf(n,n,a,n,p,info);
}
template<class T> void run(INTEGER n, bool pot, bool complex, bool upper) {
    std::mt19937 rng(20261101);
    std::vector<T> a(size_t(n)*n), g(a.size());
    for (INTEGER j=0;j<n;++j) for (INTEGER i=0;i<=j;++i) {
        T x = i==j?T(2*n):scalar<T>(rng);
        a[i+j*n]=x; a[j+i*n]=adjoint(x);
    }
    std::vector<INTEGER> piv(n);
    double best=1e100;
    int repeats=n<=2000?3:1;
    for (int r=0;r<repeats;++r) {
        g=a; INTEGER info=-1;
        auto start=std::chrono::steady_clock::now();
        factor(pot,upper,n,g.data(),piv.data(),info);
#ifdef MPLAPACK_FACTOR_CUDA
        cudaError_t status = cudaDeviceSynchronize();
#endif
        double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
#ifdef MPLAPACK_FACTOR_CUDA
        if (status != cudaSuccess)
            throw std::runtime_error(std::string("CUDA runtime (timing invalid; possible CPU fallback): ")+cudaGetErrorString(status));
#endif
        if (info) throw std::runtime_error("factorization info="+std::to_string(info));
        best=std::min(best,seconds);
    }
    double flops=(pot?1.0:2.0)*double(n)*double(n)*double(n)/3*(complex?4:1);
    std::printf("%s,%ld,%d,%.9g,%.9g\n",pot?(complex?"Cpotrf":"Rpotrf"):(complex?"Cgetrf":"Rgetrf"),long(n),repeats,best,flops/best/1e9);
    std::fflush(stdout);
}
}
int main(int argc, char **argv) {
    try {
        INTEGER limit=2048; bool pot=false, complex=false, upper=false;
        std::set<INTEGER> sizes;
        for (int i=1;i<argc;++i) {
            if (!std::strcmp(argv[i],"--potrf")) pot=true;
            else if (!std::strcmp(argv[i],"--complex")) complex=true;
            else if (!std::strcmp(argv[i],"--upper")) upper=true;
            else if ((!std::strcmp(argv[i],"--max") || !std::strcmp(argv[i],"--n")) && i+1<argc) {
                bool single=!std::strcmp(argv[i],"--n");
                std::string value=argv[++i]; size_t end=0;
                long n=std::stol(value,&end);
                if (end!=value.size() || n<1 || n>100000) throw std::runtime_error("invalid size");
                if (single) sizes.insert(n); else limit=n;
            } else throw std::runtime_error("usage: Rfactor [--potrf] [--complex] [--upper] [--max N | --n N ...]");
        }
#ifdef MPLAPACK_FACTOR_CUDA
        int devices=0;
        if (cudaGetDeviceCount(&devices)!=cudaSuccess || !devices) return 77;
#endif
#if defined(MPLAPACK_BUILD_WITH_MPFR)
        mpfr_set_default_prec(512);
#endif
        if (sizes.empty()) {
            for (INTEGER n=1;n<=limit;n*=2) sizes.insert(n);
            for (INTEGER n=13;n<=limit;n+=13) sizes.insert(n);
        }
        std::puts("routine,n,repeats,seconds,gflops");
        for (INTEGER n:sizes) {
            if (complex) run<COMPLEX>(n,pot,true,upper); else run<REAL>(n,pot,false,upper);
        }
    } catch (const std::exception &e) { std::fprintf(stderr,"%s\n",e.what()); return 1; }
}
