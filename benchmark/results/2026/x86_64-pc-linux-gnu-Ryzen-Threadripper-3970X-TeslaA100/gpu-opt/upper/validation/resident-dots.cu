// Bench the actual resident fixed-complex Upper update jobs (not public BLAS).
#include "Cfactor_device_body.cuh"
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <limits>
namespace mplapack_blas_dev {
long mplapack_mpfr_rt_prec_h=512;
__device__ long mplapack_mpfr_rt_prec_d=512;
}
using namespace mplapack_mpfr_cuda;
using namespace mplapack_blas_dev;
static void check(cudaError_t e) { if(e!=cudaSuccess) { std::fprintf(stderr,"CUDA: %s\n",cudaGetErrorString(e)); std::exit(1); } }
static unsigned long long state=7;
static unsigned long long random_bits() { state^=state<<13;state^=state>>7;state^=state<<17;return state; }
int main(int argc,char **argv) {
    if(argc!=3) return 2;
    const bool rank=std::atoi(argv[1])!=0;
    const long n=std::atol(argv[2]); if(n<1) return 2;
    using E=complex_factor_engine<factor_fixed_tag<512>>; using V=E::V;
    std::vector<V> h(n*n);
    for(auto &v:h) for(auto *r:{&v.re_,&v.im_}) {
        r->valid=true; r->v.sign=(random_bits()&1)?1:-1;r->v.exp=long(random_bits()%5)-2;
        for(auto &limb:r->v.m) limb=random_bits();
        r->v.m[r->v.N-1]|=1ULL<<63;
    }
    factor_buffer<V> a,b,c,original,constants;factor_buffer<int> failed;
    if(a.allocate(h.size()) || b.allocate(h.size()) || c.allocate(h.size()) || original.allocate(h.size()) || constants.allocate(2) || failed.allocate(1)) return 1;
    if(a.upload(h.data(),h.size()) || b.upload(h.data(),h.size()) || original.upload(h.data(),h.size())) return 1;
    V cs[2]={V(1),V(-1)};if(constants.upload(cs,2)) return 1;
    cudaEvent_t start,end;check(cudaEventCreate(&start));check(cudaEventCreate(&end));
    double best=std::numeric_limits<double>::infinity();
    for(int rep=-1;rep<3;++rep) {
        check(cudaMemcpy(c.p,original.p,h.size()*sizeof(V),cudaMemcpyDeviceToDevice));
        check(cudaMemset(failed.p,0,sizeof(int)));
        check(cudaEventRecord(start));
#ifdef UPPER_AFTER
        upper_update_job<E> job{E{}, {a.p,n},{rank?a.p:b.p,n},{c.p,n},{constants.p,2},n,n,n,rank};
        check(static_cast<cudaError_t>(upper_update_device(job,nullptr,failed.p)));
#else
        factor_update_job<E> job{};job.a={a.p,n};job.b={b.p,n};job.c={c.p,n};job.constants={constants.p,2};
        job.op=rank?0:1;job.g={0,1,n,n,n,n,n,n,0,1};
        job.s.op=L3_SYRK;job.s.upper=1;job.s.trans=1;job.s.n=n;job.s.k=n;job.s.beta_one=1;
        const long blocks=std::min(factor_blocks,(n*n+factor_threads-1)/factor_threads);
        factor_kernel<<<blocks,factor_threads>>>(job,n*n,failed.p,0);check(cudaGetLastError());
#endif
        check(cudaEventRecord(end));check(cudaEventSynchronize(end));
        int bad=0;if(failed.download(&bad,1) || bad) return 1;
        float ms=0;check(cudaEventElapsedTime(&ms,start,end));if(rep>=0 && ms/1000.0<best) best=ms/1000.0;
    }
    std::printf("mpfr512,%s,%s,%ld,gpu,%.12g,%.9g,3\n",rank?"cherk":"cgemm",rank?"UC":"CN",n,best,(rank?4.0:8.0)*n*n*n/best/1e9);
    if(const char *path=std::getenv("UPPER_DUMP")) {
        if(c.download(h.data(),h.size())) return 1;
        FILE *out=std::fopen(path,"w");if(!out) return 1;
        for(const auto &v:h) for(const auto *r:{&v.re_,&v.im_}) {
            std::fprintf(out,"%d %ld",r->v.sign,(long)r->v.exp);
            for(auto limb:r->v.m) std::fprintf(out," %016llx",(unsigned long long)limb);
            std::fputc('\n',out);
        }
        std::fclose(out);
    }
    check(cudaEventDestroy(start));check(cudaEventDestroy(end));
}
