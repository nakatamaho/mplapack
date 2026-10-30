// Isolated complex panel phases, not a hybrid public factorization path.
#include <mpblas.h>
#include <mplapack.h>
#include <cuda_runtime_api.h>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "../mpfr_pack_cuda.h"
#include "../Cfactor_kernel_cuda.h"
using namespace mplapack_mpfr_cuda;
using namespace mplapack_mpfr_cuda_host;
using namespace mplapack_blas_dev;
namespace {
using Fixed=factor_complex_value<512>;
bool pack_fixed(const std::vector<COMPLEX> &a,std::vector<Fixed> &p) {
    p.resize(a.size());
    for(size_t i=0;i<a.size();++i) {
        auto z=a[i].mpc_data();
        mpfr_srcptr src[2]={mpc_realref(z),mpc_imagref(z)};
        factor_real<512> *dst[2]={&p[i].re_,&p[i].im_};
        for(int k=0;k<2;++k) {
            if(!pack<512>(src[k],dst[k]->v,RT_DEVICE_EMAX)) return false;
            if(dst[k]->v.is_zero()) dst[k]->v.sign=mpfr_signbit(src[k])?-1:1;
        }
    }
    return true;
}
void unpack_fixed(const std::vector<Fixed> &p,std::vector<COMPLEX> &a) {
    for(size_t i=0;i<p.size();++i) {
        auto z=a[i].mpc_data(); mpfr_ptr dst[2]={mpc_realref(z),mpc_imagref(z)};
        const factor_real<512> *src[2]={&p[i].re_,&p[i].im_};
        for(int k=0;k<2;++k) {
            if(src[k]->v.is_zero()) mpfr_set_zero(dst[k],src[k]->v.sign);
            else unpack<512>(src[k]->v,dst[k]);
        }
    }
}
void pack_runtime(mpfr_srcptr x,rreal<8> &v) {
    v.kind=mpfr_custom_get_kind(x); v.exp=0;
    for(int i=0;i<8;++i) v.m[i]=0;
    if(!mpfr_zero_p(x)) {
        v.exp=mpfr_get_exp(x);
        auto limbs=static_cast<const mp_limb_t*>(mpfr_custom_get_significand(x));
        for(int i=0;i<8;++i) v.m[i]=limbs[i];
    }
}
}
int main(int argc,char **argv) {
    if(argc!=4) { std::fprintf(stderr,"usage: %s m getrf|potrf fixed|runtime|host\n",argv[0]); return 2; }
    if((std::strcmp(argv[2],"getrf") && std::strcmp(argv[2],"potrf")) ||
       (std::strcmp(argv[3],"fixed") && std::strcmp(argv[3],"runtime") && std::strcmp(argv[3],"host"))) return 2;
    setenv("MPLAPACK_MPFR_CUDA","0",1); // CPU-only oracle in the hybrid comparator.
    mpfrxx::set_default_precision_bits(512); mpfrxx::set_default_mpc_precision_bits(512);
    bool chol=!std::strcmp(argv[2],"potrf"); long m=std::atol(argv[1]),n=chol?m:64;
    if(m<1) return 2;
    std::vector<COMPLEX> a(m*n); std::vector<long> piv(std::min(m,n)); long info=-1;
    gmp_randstate_t rng; gmp_randinit_default(rng); gmp_randseed_ui(rng,42);
    for(long j=0;j<n;++j) for(long i=0;i<m;++i) {
        auto z=a[i+j*m].mpc_data();
        mpfr_urandomb(mpc_realref(z),rng); mpfr_urandomb(mpc_imagref(z),rng);
        a[i+j*m]=(a[i+j*m]-COMPLEX(.5,.5))/COMPLEX(8);
        if(i==j) { a[i+j*m]+=COMPLEX(m+1); mpfr_set_zero(mpc_imagref(z),1); }
    }
    if(chol) for(long j=0;j<n;++j) for(long i=j+1;i<m;++i) a[j+i*m]=conj(a[i+j*m]);
    gmp_randclear(rng);
    REAL safe=Rlamch("S"); factor_shape s{chol,false,m,n,0}; factor_times t;
    auto start=std::chrono::steady_clock::now(); int rc=0;
    if(!std::strcmp(argv[3],"runtime")) {
        std::vector<factor_runtime_complex<8>> p(a.size()); factor_runtime_complex<8> sf;
        for(size_t i=0;i<a.size();++i) {
            auto z=a[i].mpc_data(); pack_runtime(mpc_realref(z),p[i].re_); pack_runtime(mpc_imagref(z),p[i].im_);
        }
        pack_runtime(safe.get_mpfr_t(),sf.re_);
        rc=factor_complex<factor_runtime_tag<8>>(s,512,p.data(),sf,piv.data(),info,t);
    } else {
        std::vector<Fixed> p; Fixed sf;
        if(!pack_fixed(a,p) || !pack<512>(safe.get_mpfr_t(),sf.re_.v,RT_DEVICE_EMAX)) return 1;
        if(!std::strcmp(argv[3],"fixed")) rc=factor_complex<factor_fixed_tag<512>>(s,512,p.data(),sf,piv.data(),info,t);
        else {
            Fixed *device=nullptr; size_t bytes=p.size()*sizeof(Fixed);
            if(cudaMalloc(reinterpret_cast<void**>(&device),bytes) || cudaMemcpy(device,p.data(),bytes,cudaMemcpyHostToDevice) || cudaDeviceSynchronize()) return 1;
            auto begin=std::chrono::steady_clock::now();
            if(cudaMemcpy(p.data(),device,bytes,cudaMemcpyDeviceToHost)) return 1;
            unpack_fixed(p,a);
            INTEGER result=-1; std::vector<INTEGER> host_piv(std::min(m,n));
            if(chol) Cpotrf2("L",n,a.data(),m,result); else Cgetrf2(m,n,a.data(),m,host_piv.data(),result);
            info=result;
            if(!pack_fixed(a,p) || cudaMemcpy(device,p.data(),bytes,cudaMemcpyHostToDevice) || cudaDeviceSynchronize()) return 1;
            t.factor=std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();
            cudaFree(device);
        }
    }
    double total=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    std::printf("%s,%s,%ld,%ld,phase=%.9f,total=%.9f,rc=%d,info=%ld\n",argv[2],argv[3],m,n,t.factor,total,rc,info);
    return rc || info?1:0;
}
