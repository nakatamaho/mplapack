// Isolated 512-bit panel phases: device fixed, device MPFR, and a MAGMA-style
// D2H / unpack / host MPFR panel / pack / H2D phase. Trailing updates are common
// to the three alternatives and are deliberately excluded, not estimated.
#include <mpblas.h>
#include <mplapack.h>
#include <cuda_runtime_api.h>
#include <chrono>
#include <cstdio>
#include <cstring>
#include "../mpfr_pack_cuda.h"
#include "../Rfactor_kernel_cuda.h"
using namespace mplapack_mpfr_cuda;
using namespace mplapack_mpfr_cuda_host;
int main(int argc,char **argv) {
    if(argc!=4) { std::fprintf(stderr,"usage: %s m getrf|potrf fixed|runtime|host\n",argv[0]); return 2; }
    setenv("MPLAPACK_MPFR_CUDA","0",1); // Host comparator must use CPU BLAS.
    mpfrxx::set_default_precision_bits(512);
    long m=std::atol(argv[1]),n=!std::strcmp(argv[2],"potrf")?m:64;
    bool chol=!std::strcmp(argv[2],"potrf");
    if(m<1) return 2;
    std::vector<REAL> a(m*n); std::vector<long> piv(std::min(m,n)); long info=-1;
    gmp_randstate_t rng; gmp_randinit_default(rng); gmp_randseed_ui(rng,42);
    for(long j=0;j<n;++j) for(long i=0;i<m;++i) {
        mpfr_urandomb(a[i+j*m].get_mpfr_t(),rng);
        a[i+j*m]=(a[i+j*m]-REAL(0.5))/REAL(8);
        if(i==j) a[i+j*m]+=REAL(m+1);
    }
    if(chol) for(long j=0;j<n;++j) for(long i=j+1;i<m;++i) a[j+i*m]=a[i+j*m];
    gmp_randclear(rng);
    REAL safe=Rlamch("S"); factor_shape s{chol,false,m,n,0}; factor_times t;
    auto start=std::chrono::steady_clock::now(); int rc=0;
    if(!std::strcmp(argv[3],"runtime")) {
        rt_storage p,sf; sf.resize(1,8);
        if(!rt_pack_matrix(a.data(),m,m,n,512,8,p,RT_DEVICE_EMAX) ||
           !rt_pack_value(safe.get_mpfr_t(),512,8,sf,0,RT_DEVICE_EMAX)) return 1;
        rc=factor_runtime(s,512,p.pack,sf.pack,piv.data(),info,t);
    } else {
        std::vector<real_t<512>> p; real_t<512> sf;
        if(!pack_matrix<512>(a.data(),m,m,n,p,RT_DEVICE_EMAX) || !pack<512>(safe.get_mpfr_t(),sf,RT_DEVICE_EMAX)) return 1;
        if(!std::strcmp(argv[3],"fixed")) rc=factor_fixed<512>(s,p.data(),sf,piv.data(),info,t);
        else {
            real_t<512> *device=nullptr; size_t bytes=p.size()*sizeof(*device);
            if(cudaMalloc(reinterpret_cast<void**>(&device),bytes)!=cudaSuccess || cudaMemcpy(device,p.data(),bytes,cudaMemcpyHostToDevice)!=cudaSuccess || cudaDeviceSynchronize()!=cudaSuccess) return 1;
            auto begin=std::chrono::steady_clock::now();
            if(cudaMemcpy(p.data(),device,bytes,cudaMemcpyDeviceToHost)!=cudaSuccess) return 1;
            for(size_t i=0;i<p.size();++i) unpack<512>(p[i],a[i].get_mpfr_t());
            INTEGER result=-1; std::vector<INTEGER> host_piv(std::min(m,n));
            if(chol) Rpotrf2("L",n,a.data(),m,result); else Rgetrf2(m,n,a.data(),m,host_piv.data(),result);
            info=result;
            if(!pack_matrix<512>(a.data(),m,m,n,p,RT_DEVICE_EMAX)) return 1;
            if(cudaMemcpy(device,p.data(),bytes,cudaMemcpyHostToDevice)!=cudaSuccess || cudaDeviceSynchronize()!=cudaSuccess) return 1;
            t.factor=std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();
            cudaFree(device);
        }
    }
    double total=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    std::printf("%s,%s,%ld,%ld,phase=%.9f,total=%.9f,rc=%d,info=%ld\n",argv[2],argv[3],m,n,t.factor,total,rc,info);
    return rc || info?1:0;
}
