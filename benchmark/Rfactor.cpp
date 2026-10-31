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
#if defined(MPLAPACK_BUILD_WITH_DD)
#include <mpblas_dd_cuda.h>
#define ENTRY(R) R##_dd_cuda
#elif defined(MPLAPACK_BUILD_WITH_TD)
#include <mpblas_td_cuda.h>
#define ENTRY(R) R##_td_cuda
#elif defined(MPLAPACK_BUILD_WITH_QD)
#include <mpblas_qd_cuda.h>
#define ENTRY(R) R##_qd_cuda
#elif defined(MPLAPACK_BUILD_WITH_MPFR)
#include <mpblas_mpfr_cuda.h>
#define ENTRY(R) R##_mpfr_cuda
#endif
#include "Rsolver_legacy_bench.h"
#endif

namespace {
#ifdef MPLAPACK_FACTOR_CUDA
bool legacy_path=false;
#define SOLVE_CPU(R) (legacy_path?solver_legacy_##R:R)
#else
#define SOLVE_CPU(R) R
#endif
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
#define SOLVE(H,K) \
void solve(bool pot,bool upper,INTEGER n,INTEGER rhs,H *a,INTEGER *p,H *b,INTEGER &info,bool entry,bool driver) { \
    if(entry) { \
        SOLVE_ENTRY(K) \
    } else if(driver) { if(pot) SOLVE_CPU(K##posv)(upper?"U":"L",n,rhs,a,n,b,n,info); else SOLVE_CPU(K##gesv)(n,rhs,a,n,p,b,n,info); } \
    else { if(pot) SOLVE_CPU(K##potrs)(upper?"U":"L",n,rhs,a,n,b,n,info); else SOLVE_CPU(K##getrs)("N",n,rhs,a,n,p,b,n,info); } \
}
#ifdef MPLAPACK_FACTOR_CUDA
#define SOLVE_ENTRY(K) if(!(driver?(pot?ENTRY(K##posv)(upper?"U":"L",n,rhs,a,n,b,n,info):ENTRY(K##gesv)(n,rhs,a,n,p,b,n,info)):(pot?ENTRY(K##potrs)(upper?"U":"L",n,rhs,a,n,b,n,info):ENTRY(K##getrs)("N",n,rhs,a,n,p,b,n,info)))) throw std::runtime_error("entry rejected");
#else
#define SOLVE_ENTRY(K) throw std::runtime_error("--entry requires CUDA benchmark");
#endif
SOLVE(REAL,R)
SOLVE(COMPLEX,C)
#undef SOLVE
#undef SOLVE_ENTRY
#undef SOLVE_CPU
template<class T> void run(INTEGER n, bool pot, bool complex, bool upper,INTEGER requested_rhs,bool entry,bool driver,int requested_repeats) {
    std::mt19937 rng(20261101);
    std::vector<T> a(size_t(n)*n), g(a.size());
    for (INTEGER j=0;j<n;++j) for (INTEGER i=0;i<=j;++i) {
        T x = i==j?T(2*n):scalar<T>(rng);
        a[i+j*n]=x; a[j+i*n]=adjoint(x);
    }
    INTEGER rhs=requested_rhs==0?n:requested_rhs;
    std::vector<T> original_b(rhs<0?0:size_t(n)*rhs),b;
    for(auto &v:original_b) v=scalar<T>(rng);
    std::vector<INTEGER> piv(n);
    // Standalone solves accept synthetic nonsingular triangular factors.
    // This avoids timing (or repeatedly paying for) an unrelated factorization.
    if(rhs>=0 && !driver) for(INTEGER i=0;i<n;++i) piv[i]=n;
    double best=1e100;
    int repeats=requested_repeats>0?requested_repeats:(n<=2000?3:1);
    for (int r=0;r<repeats;++r) {
        g=a; b=original_b; INTEGER info=-1;
        auto start=std::chrono::steady_clock::now();
        if(rhs<0) factor(pot,upper,n,g.data(),piv.data(),info);
        else solve(pot,upper,n,rhs,g.data(),piv.data(),b.data(),info,entry,driver);
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
    double flops=((rhs>=0 && !driver?0:(pot?1.0:2.0)*double(n)*double(n)*double(n)/3)+(rhs<0?0:2.0*n*n*rhs))*(complex?4:1);
    const char *routine=rhs<0?(pot?(complex?"Cpotrf":"Rpotrf"):(complex?"Cgetrf":"Rgetrf")):driver?(pot?(complex?"Cposv":"Rposv"):(complex?"Cgesv":"Rgesv")):(pot?(complex?"Cpotrs":"Rpotrs"):(complex?"Cgetrs":"Rgetrs"));
    if(rhs<0) std::printf("%s,%ld,%d,%.9g,%.9g\n",routine,long(n),repeats,best,flops/best/1e9);
    else std::printf("%s,%ld,%ld,%d,%.9g,%.9g\n",routine,long(n),long(rhs),repeats,best,flops/best/1e9);
    std::fflush(stdout);
}
}
int main(int argc, char **argv) {
    try {
        INTEGER limit=2048,rhs=-1; int repeats=-1; bool pot=false, complex=false, upper=false,entry=false,driver=true;
        std::set<INTEGER> sizes;
        for (int i=1;i<argc;++i) {
            if (!std::strcmp(argv[i],"--potrf")) pot=true;
            else if (!std::strcmp(argv[i],"--gesv")) { pot=false; rhs=1; }
            else if (!std::strcmp(argv[i],"--posv")) { pot=true; rhs=1; }
            else if (!std::strcmp(argv[i],"--getrs")) { pot=false; rhs=1; driver=false; }
            else if (!std::strcmp(argv[i],"--potrs")) { pot=true; rhs=1; driver=false; }
            else if (!std::strcmp(argv[i],"--entry")) entry=true;
            else if (!std::strcmp(argv[i],"--legacy")) {
#ifdef MPLAPACK_FACTOR_CUDA
                legacy_path=true;
#else
                throw std::runtime_error("--legacy requires CUDA benchmark");
#endif
            }
            else if (!std::strcmp(argv[i],"--repeat") && i+1<argc) {
                std::string v=argv[++i]; size_t end=0; long r=std::stol(v,&end);
                if(end!=v.size() || r<1 || r>100) throw std::runtime_error("invalid repeat");
                repeats=int(r);
            }
            else if (!std::strcmp(argv[i],"--nrhs") && i+1<argc) {
                std::string v=argv[++i]; size_t end=0; long r=std::stol(v,&end);
                if(end!=v.size() || r<0 || r>100000) throw std::runtime_error("invalid nrhs (0 means n)");
                rhs=r;
            }
            else if (!std::strcmp(argv[i],"--complex")) complex=true;
            else if (!std::strcmp(argv[i],"--upper")) upper=true;
            else if ((!std::strcmp(argv[i],"--max") || !std::strcmp(argv[i],"--n")) && i+1<argc) {
                bool single=!std::strcmp(argv[i],"--n");
                std::string value=argv[++i]; size_t end=0;
                long n=std::stol(value,&end);
                if (end!=value.size() || n<1 || n>100000) throw std::runtime_error("invalid size");
                if (single) sizes.insert(n); else limit=n;
            } else throw std::runtime_error("usage: Rfactor [--potrf | --gesv | --posv | --getrs | --potrs] [--nrhs N (0=n)] [--entry | --legacy] [--repeat N] [--complex] [--upper] [--max N | --n N ...]");
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
        std::puts(rhs<0?"routine,n,repeats,seconds,gflops":"routine,n,nrhs,repeats,seconds,gflops");
        for (INTEGER n:sizes) {
            if (complex) run<COMPLEX>(n,pot,true,upper,rhs,entry,driver,repeats); else run<REAL>(n,pot,false,upper,rhs,entry,driver,repeats);
        }
    } catch (const std::exception &e) { std::fprintf(stderr,"%s\n",e.what()); return 1; }
}
