// Internal device-resident Cholesky: reference arithmetic, two-stream lookahead.
#ifndef MPLAPACK_POTRF_GPU_COMMON_CUH
#define MPLAPACK_POTRF_GPU_COMMON_CUH
#include "../../../mpblas/optimized/common/Rlapack_cuda_common.cuh"

namespace mplapack_lapack_gpu {
// Only scheduling differs from the reference. A panel-ready event releases
// the far update; update-done is waited on before touching the next block row.
class chol_schedule {
  public:
    int error=0;
    bool far=false;
#if defined(__CUDACC__)
    cudaStream_t streams[2]{};
    cudaEvent_t ready{},done{};
    chol_schedule() {
        error=cudaStreamCreateWithFlags(&streams[0],cudaStreamNonBlocking);
        if(!error) error=cudaStreamCreateWithFlags(&streams[1],cudaStreamNonBlocking);
        if(!error) error=cudaEventCreateWithFlags(&ready,cudaEventDisableTiming);
        if(!error) error=cudaEventCreateWithFlags(&done,cudaEventDisableTiming);
    }
    ~chol_schedule() {
        finish();
        cudaEventDestroy(ready); cudaEventDestroy(done);
        cudaStreamDestroy(streams[0]); cudaStreamDestroy(streams[1]);
    }
    cudaStream_t stream() const { return streams[far?1:0]; }
    void join_update() {
        if(!error) error=cudaEventRecord(done,streams[1]);
        if(!error) error=cudaStreamWaitEvent(streams[0],done,0);
    }
    void release_update() {
        if(!error) error=cudaEventRecord(ready,streams[0]);
        if(!error) error=cudaStreamWaitEvent(streams[1],ready,0);
    }
    void finish() {
        for(auto s:streams) if(s) { int e=cudaStreamSynchronize(s); if(!error) error=e; }
    }
#else
    void join_update() {}
    void release_update() {}
    void finish() {}
#endif
    template<class Job> void launch(Job j,long units) {
        if(error || units<=0) return;
#if defined(__CUDACC__)
        k_elements<<<(unsigned)((units+63)/64),64,0,stream()>>>(j,units);
        error=cudaGetLastError();
#else
        for(long i=0;i<units;++i) j(i);
#endif
    }
    long feedback(long *status) {
        long result=0;
#if defined(__CUDACC__)
        if(!error) error=cudaMemcpyAsync(&result,status,sizeof(result),cudaMemcpyDeviceToHost,stream());
        if(!error) error=cudaStreamSynchronize(stream());
#else
        if(!error) result=*status;
#endif
        return result;
    }
};
template<class Context> static long chol_recursive(Context &ctx,bool upper,long n,typename Context::E *a,long ld) {
    if(ctx.error || !n) return 0;
    if(n==1) return ctx.base(a);
    long n1=n/2,n2=n-n1;
    long info=chol_recursive(ctx,upper,n1,a,ld);
    if(ctx.error || info) return info;
    auto *off=upper?a+n1*ld:a+n1;
    ctx.solve(upper,upper?n1:n2,upper?n2:n1,a,ld,off);
    ctx.rank(upper,n2,n1,off,ld,a+n1+n1*ld);
    info=chol_recursive(ctx,upper,n2,a+n1+n1*ld,ld);
    return info?info+n1:0;
}
template<class Ops> class real_cholesky : public chol_schedule {
    using F=typename Ops::F;
    buffer<F> work[2],constants;
    buffer<long> status;
    void blas(mplapack_blas::gpu_lane_args<Ops> g) {
        if(error) return;
#if defined(__CUDACC__)
        using namespace mplapack_blas;
        long lanes=blas_lanes(g.s);
        if(lanes<=0) return;
        if(syrk_tiled<Ops>(g.s)) {
            unsigned nb=(unsigned)((g.s.n+BLAS_TILE-1)/BLAS_TILE);
            k_syrk_tiled<Ops><<<dim3(nb,nb),dim3(BLAS_TILE,BLAS_TILE),0,stream()>>>(g);
        } else k_blas_lanes<Ops><<<(unsigned)((lanes+63)/64),64,0,stream()>>>(g,lanes);
        error=cudaGetLastError();
#else
        error=mplapack_blas::blas_device<Ops>(g);
#endif
    }
  public:
    using E=dreal<Ops>;
    real_cholesky(long n,long nb) {
        if(!error) error=work[0].allocate((size_t)n*nb);
        if(!error) error=work[1].allocate((size_t)n*nb);
        if(!error) error=constants.allocate(2);
        if(!error) error=status.allocate(1);
        F c[2]={Ops::from_d(1),Ops::from_d(-1)};
        if(!error) error=transfer(constants.data(),c,2,false);
    }
    ~real_cholesky() { finish(); }
    long base(E *a) {
        launch(base_chol<Ops>{a,status.data()},1);
        return feedback(status.data());
    }
    void rank(bool upper,long n,long k,E *a,long ld,E *c) {
        if(!k) return; // Reference quick return, including signed zero limbs.
        mplapack_blas::gpu_lane_args<Ops> g{};
        g.s.op=mplapack_blas::L3_SYRK; g.s.upper=upper; g.s.trans=upper;
        g.s.n=n; g.s.k=k; g.s.beta_one=1;
        g.a=reinterpret_cast<F*>(a); g.c=reinterpret_cast<F*>(c);
        g.al=constants.data()+1; g.be=constants.data(); g.lda=ld; g.ldc=ld;
        blas(g);
    }
    void solve(bool upper,long m,long n,E *a,long ld,E *b) {
        mplapack_blas::gpu_lane_args<Ops> g{};
        g.s.op=mplapack_blas::L3_TRSM; g.s.left=upper; g.s.upper=upper;
        g.s.trans=1; g.s.nounit=1; g.s.alpha_one=1; g.s.m=m; g.s.n=n;
        g.a=reinterpret_cast<F*>(a); g.c=reinterpret_cast<F*>(b);
        g.al=constants.data(); g.lda=ld; g.ldc=ld;
        blas(g);
    }
    void update(bool upper,long m,long n,long k,E *a,E *b,E *c,long ld) {
        if(error || !k || !m || !n) return;
        using namespace mplapack_gemm;
        shape s{};
        s.m=m; s.n=n; s.k=k; s.lda=ld; s.ldb=ld; s.ldc=ld;
        s.nota=!upper; s.notb=upper;
        gpu_scalars<F> sc{Ops::from_d(-1),Ops::from_d(1),false,true};
        F *aa=reinterpret_cast<F*>(a),*bb=reinterpret_cast<F*>(b),*cc=reinterpret_cast<F*>(c),*scratch=work[far?1:0].data();
#if defined(__CUDACC__)
        if(s.nota) k_alpha_opB<Ops><<<gpu_grid2(s.k,s.n),dim3(GPU_TILE,GPU_TILE),0,stream()>>>(s,sc.alpha,bb,scratch);
        if(cudaGetLastError()!=cudaSuccess) { error=1; return; }
        k_tiled<Ops><<<gpu_grid2(s.m,s.n),dim3(GPU_TILE,GPU_TILE),0,stream()>>>(s,sc,aa,s.nota?scratch:bb,cc,0);
        error=cudaGetLastError();
#else
        error=gemm_device<Ops>(s,sc,aa,bb,cc,scratch);
#endif
    }
};
template<class Ops> struct complex_chol_base {
    using R=dreal<Ops>; using C=dcomplex<R>;
    C *a; long *info;
    MPLAPACK_DEV_HD void operator()(long) const {
        R x=a[0].real();
        *info=(x<=R(0) || Ops::is_nan(x.v))?1:0;
        if(!*info) a[0]=C(R(square_root<Ops>(x.v)),R(0));
    }
};
template<class Ops> class complex_cholesky : public chol_schedule {
    buffer<long> status;
  public:
    using R=dreal<Ops>; using E=dcomplex<R>;
    complex_cholesky(long,long) { if(!error) error=status.allocate(1); }
    ~complex_cholesky() { finish(); }
    long base(E *a) {
        launch(complex_chol_base<Ops>{a,status.data()},1);
        return feedback(status.data());
    }
    void rank(bool upper,long n,long k,E *a,long ld,E *c) {
        if(!k) return;
        job_Cherk<Ops> j{};
        j.uplo[0]=upper?'U':'L'; j.tr[0]=upper?'C':'N';
        j.n=n; j.k=k; j.lda=ld; j.ldc=ld;
        j.alpha=R(-1); j.beta=R(1); j.A=a; j.C=c;
        launch(j,j.units());
    }
    void solve(bool upper,long m,long n,E *a,long ld,E *b) {
        job_Ctrsm<Ops> j{};
        j.side[0]=upper?'L':'R'; j.uplo[0]=upper?'U':'L'; j.ta[0]='C'; j.diag[0]='N';
        j.m=m; j.n=n; j.lda=ld; j.ldb=ld; j.alpha=E(1,0); j.A=a; j.B=b;
        launch(j,j.units());
    }
    void update(bool upper,long m,long n,long k,E *a,E *b,E *c,long ld) {
        if(!k) return;
        job_Cgemm<Ops> j{};
        j.ta[0]=upper?'C':'N'; j.tb[0]=upper?'N':'C';
        j.m=m; j.n=n; j.k=k; j.lda=ld; j.ldb=ld; j.ldc=ld;
        j.alpha=E(-1,0); j.beta=E(1,0); j.A=a; j.B=b; j.C=c;
        launch(j,j.units());
    }
};
template<class Context> static long cholesky(Context &ctx,bool upper,long n,long nb,typename Context::E *a,long ld) {
    if(nb<=1 || nb>=n) return chol_recursive(ctx,upper,n,a,ld);
    for(long j=0;j<n;j+=nb) {
        long jb=std::min(nb,n-j);
        auto *diag=a+j+j*ld;
        ctx.rank(upper,jb,j,upper?a+j*ld:a+j,ld,diag);
        long info=chol_recursive(ctx,upper,jb,diag,ld);
        // Complete the previous far block row even on numerical failure:
        // reference potrf has already computed that row before this diagonal.
        if(ctx.error || info) { ctx.finish(); return info?info+j:0; }
        ctx.join_update();
        // First update/solve the next panel, then release the far update.
        // The next diagonal can factor while the far update is still running.
        for(int part=0;part<2;++part) {
            long offset=j+jb+(part?std::min(nb,n-j-jb):0);
            long rest=part?n-offset:std::min(nb,n-offset);
            ctx.far=part!=0;
            if(rest) {
                auto *off=upper?a+j+offset*ld:a+offset+j*ld;
                ctx.update(upper,upper?jb:rest,upper?rest:jb,j,
                           upper?a+j*ld:a+offset,upper?a+offset*ld:a+j,off,ld);
                ctx.solve(upper,upper?jb:rest,upper?rest:jb,diag,ld,off);
            }
            if(!part) ctx.release_update();
        }
        ctx.far=false;
        if(ctx.error) return 0;
    }
    ctx.finish();
    return 0;
}
} // namespace mplapack_lapack_gpu
#endif
