// Goal 11 benchmark extension. BSD-2-Clause; see Rkernel_bench.cpp.
// Upper-case public routine names select this path; legacy timings are untouched.
#ifndef MPLAPACK_KERNEL_NEW_H
#define MPLAPACK_KERNEL_NEW_H
#if defined(KERNEL_BENCH_BF_CUDA)
#include "../mpblas/optimized/common/Rtrisolve_cuda_common.cuh"
#elif defined(KERNEL_BENCH_MPFR_RESIDENT)
#include "../mpblas/optimized/mpfr/cuda/Rblas1_device_body.cuh"
#include "../mpblas/optimized/mpfr/cuda/Rtrisolve_device_body.cuh"
#include "../mpblas/optimized/mpfr/cuda/mpfr_pack_cuda.h"
#endif
#include <fstream>

namespace kernel_new {
const std::vector<std::string> names = {
    "Raxpy","Rcopy","iRamax","iCamax","Rdot","Cdotc","Cdotu",
    "Rasum","RCasum","Rnrm2","RCnrm2","Rtrsv","Rtbsv","Rtpsv",
    "Ctrsv","Ctbsv","Ctpsv","Rtrtri","Rlauum","Rpotri","Rgetri",
    "Ctrtri","Clauum","Cpotri","Cgetri"};
bool accepts(const std::string &s) { return std::find(names.begin(),names.end(),s)!=names.end(); }
struct shape {
    int id; long n, nbtri=64,nblau=64,nb=64;
    bool complex() const { return id==3 || id==5 || id==6 || id==8 || id==10 || (id>=14 && id<=16) || id>=21; }
    bool vector() const { return id<11; }
    bool solve() const { return id>=11 && id<=16; }
    bool inverse() const { return id>=17; }
    bool reduction() const { return id>=4 && id<=10; }
    int storage() const { return (id-11)%3; }
    int inverse_kind() const { return (id-17)%4; }
    size_t acount() const { return inverse()?size_t(n)*n:solve()?(storage()==2?size_t(n)*(n+1)/2:size_t(n)*(storage()==1?65:n)):0; }
    long index(long i,long j) const { return storage()==2?j*(j+1)/2+i:storage()==1?64+i-j+65*j:i+j*n; }
    double flops() const {
        double x=n, f=complex()?4:1;
        if(inverse()) { int k=inverse_kind(); return f*(k==0?x*x*x/3+2*x/3:k==1?x*x*x/3+x*x/2+x/6:k==2?2*x*x*x/3:4*x*x*x/3); }
        if(solve()) { double k=std::min(64L,n-1); return f*(storage()==1?n*(2*k+1)-k*(k+1):x*x); }
        if(id==0 || id==4) return 2*x;
        if(id==5 || id==6) return 8*x;
        if(id==7) return x;
        if(id==8) return 2*x;
        if(id==9) return 2*x;
        if(id==10) return 4*x;
        return 0; // copy and iamax are not floating-point operations
    }
};

#if defined(KERNEL_BENCH_BF_CUDA) || defined(KERNEL_BENCH_MPFR_RESIDENT)
void check(int e,const char *what) { if(e) throw std::runtime_error(std::string(what)+": "+cudaGetErrorString(cudaError_t(e))); }
struct event {
    cudaEvent_t p{};
    event() { check(cudaEventCreate(&p),"event create"); }
    ~event() { cudaEventDestroy(p); }
};
template<class T> struct buffer {
    T *p=nullptr; size_t count=0;
    buffer()=default; buffer(const buffer&)=delete;
    ~buffer() { cudaFree(p); }
    void allocate(size_t n) { count=n; if(n) check(cudaMalloc(&p,n*sizeof(T)),"allocate"); }
    void upload(const T *h) { if(count) check(cudaMemcpy(p,h,count*sizeof(T),cudaMemcpyHostToDevice),"upload"); }
    void download(T *h,size_t n) { check(cudaMemcpy(h,p,n*sizeof(T),cudaMemcpyDeviceToHost),"download"); }
};
#if defined(KERNEL_BENCH_MPFR_RESIDENT)
using Tag=mplapack_blas_dev::factor_fixed_tag<512>;
using R=mplapack_mpfr_inverse::value<Tag,false>;
using C=mplapack_mpfr_inverse::value<Tag,true>;
using Root=mplapack_mpfr_inverse::blas1_root<Tag>;
using Launcher=mplapack_mpfr_inverse::launcher;
void pack_value(const REAL &x,R &v) {
    if(!mplapack_mpfr_cuda_host::pack<512>(x.get_mpfr_t(),v.v.v,1073741823L)) throw std::runtime_error("pack512 rejected");
    if(v.v.v.is_zero()) v.v.v.sign=mpfr_signbit(x.get_mpfr_t())?-1:1;
}
void pack_value(const COMPLEX &x,C &v) { R r,i; pack_value(REAL(x.real()),r);pack_value(REAL(x.imag()),i);v.v.re_=r.v;v.v.im_=i.v; }
void unpack_value(const R &v,REAL &x) {
    if(!mplapack_mpfr_inverse::fits(v.v)) throw std::runtime_error("invalid output");
    if(v.v.v.is_zero()) mpfr_set_zero(x.get_mpfr_t(),v.v.v.sign); else mplapack_mpfr_cuda_host::unpack<512>(v.v.v,x.get_mpfr_t());
}
void unpack_value(const C &v,COMPLEX &x) { REAL r,i; unpack_value(R(v.v.re_),r);unpack_value(R(v.v.im_),i);x=COMPLEX(r,i); }
#else
using R=DR; using C=DC; using Root=mplapack_blas1::root<Ops>;
using Launcher=mplapack_blas1::launcher;
void pack_value(const REAL &x,R &v) { pack(x,reinterpret_cast<double*>(&v)); }
void pack_value(const COMPLEX &x,C &v) { pack(x,reinterpret_cast<double*>(&v)); }
void unpack_value(const R &v,REAL &x) { unpack(reinterpret_cast<const double*>(&v),x); }
void unpack_value(const C &v,COMPLEX &x) { unpack(reinterpret_cast<const double*>(&v),x); }
#endif
mplapack_blas1::operation operation(int id) {
    using namespace mplapack_blas1;
    return id==2 || id==3?iamax:id==4?dot:id==5?dotc:id==6?dotu:id==7 || id==8?asum:mplapack_blas1::norm;
}

template<class H> struct gpu {
    static constexpr bool complex=std::is_same<H,COMPLEX>::value;
    using V=typename std::conditional<complex,C,R>::type;
    using S=mplapack_blas1::state<R>;
    shape s;
#if defined(KERNEL_BENCH_MPFR_RESIDENT)
    std::unique_lock<std::mutex> lock{mplapack_mpfr_cuda_host::gpu_mutex()};
    mplapack_mpfr_cuda::mpfr_stack_guard stack{mplapack_mpfr_inverse::inverse_stack_bytes(512,true)};
    mplapack_mpfr_cuda::factor_buffer<char> arena;
    mplapack_mpfr_cuda::factor_buffer<size_t> tops;
    std::unique_ptr<mplapack_mpfr_cuda::complex_arena_guard> guard;
    using IA=mplapack_inverse_gpu::inverse_arithmetic<R,C,complex>;
    std::unique_ptr<mplapack_inverse_gpu::inverse_schedule<IA,Launcher>> inv;
#endif
#if defined(KERNEL_BENCH_BF_CUDA)
    buffer<F> alpha;
#endif
    buffer<V> a,x,y,p,work; buffer<S> tree; buffer<long> piv;
    std::vector<V> ha,hx,hy;
    std::unique_ptr<Launcher> ctx;
    event start,stop;
    template<class X> void stage(buffer<V> &d,std::vector<V> &h,const X &in) {
        h.resize(in.size());for(size_t i=0;i<in.size();++i) pack_value(in[i],h[i]);d.allocate(h.size());d.upload(h.data());
    }
    gpu(shape sh,const std::vector<H> &aa,const std::vector<H> &xx,const std::vector<H> &yy):s(sh) {
        size_t free=0,total=0;
        check(cudaMemGetInfo(&free,&total),"memory preflight");
        size_t bytes=(aa.size()+xx.size()+yy.size()+(s.solve()?size_t(s.n)+1:0)+
                      (s.inverse()?size_t(s.n)*s.nb:0))*sizeof(V)+size_t(s.n)*sizeof(long);
        if(s.vector() && s.id>=2) bytes+=size_t(complex && s.id==10?2*s.n:s.n)*sizeof(S);
#if defined(KERNEL_BENCH_MPFR_RESIDENT)
        // Stack has already been reserved by its scoped guard (inherited limit
        // is preserved). Include the bounded production arena and context too.
        check(stack.error,"stack");
        bytes+=mplapack_mpfr_cuda::factor_pool*(16384+sizeof(size_t));
#endif
        if(bytes>free || free-bytes<1048576) throw std::runtime_error("insufficient device memory");
#if defined(KERNEL_BENCH_MPFR_RESIDENT)
        check(arena.allocate(mplapack_mpfr_cuda::factor_pool*16384),"arena");
        check(tops.allocate(mplapack_mpfr_cuda::factor_pool),"tops");
        check(cudaMemset(tops.p,0,mplapack_mpfr_cuda::factor_pool*sizeof(size_t)),"tops reset");
        guard.reset(new mplapack_mpfr_cuda::complex_arena_guard(arena.p,16384,tops.p,512));check(guard->error,"arena guard");
#endif
#if defined(KERNEL_BENCH_BF_CUDA)
        if(s.id<2) { alpha.allocate(1);F one=Ops::from_d(1);alpha.upload(&one); }
#endif
        stage(a,ha,aa);stage(x,hx,xx);stage(y,hy,yy);
        if(s.solve()) p.allocate(s.n+1);
        if(s.vector() && s.id>=2) tree.allocate(complex && s.id==10?2*s.n:s.n);
#if defined(KERNEL_BENCH_MPFR_RESIDENT)
        if(s.inverse()) {
            inv.reset(new mplapack_inverse_gpu::inverse_schedule<IA,Launcher>());check(inv->error,"inverse context");
            work.allocate(size_t(s.n)*s.nb);piv.allocate(s.n);
            std::vector<long> ip(s.n);for(long i=0;i<s.n;++i) ip[i]=i+1;piv.upload(ip.data());
            check(cudaMemset(work.p,0,work.count*sizeof(V)),"work reset");
        }
#else
        if(s.inverse()) throw std::runtime_error("BF inverses are not new goal11 work");
#endif
        ctx.reset(new Launcher());check(ctx->error,"launcher");check(cudaDeviceSynchronize(),"setup complete");
    }
    double run() {
        // Restore from the prepacked input before recording either event.
        if(s.inverse()) a.upload(ha.data());else if(s.solve()) x.upload(hx.data());else if(s.id<2) y.upload(hy.data());
        cudaStream_t stream=nullptr;
#if defined(KERNEL_BENCH_MPFR_RESIDENT)
        stream=s.inverse()?inv->stream():ctx->stream();
        check(cudaMemset((s.inverse()?inv->failed.p:ctx->failed.p),0,sizeof(int)),"status reset");
#endif
        check(cudaDeviceSynchronize(),"reset complete");
        check(cudaEventRecord(start.p,stream),"start");
        long info=0;
        if(s.solve()) {
            mplapack_trisolve::shape t{s.n,64,s.storage()==1?65:s.n,s.storage(),true,true,false,false};
            using A=typename std::conditional<complex,mplapack_trisolve::Carithmetic<V>,mplapack_trisolve::Rarithmetic<V>>::type;
            mplapack_trisolve::schedule<V,A>(*ctx,t,a.p,x.p,p.p);
        } else if(s.vector() && s.id>=2) mplapack_blas1::schedule<R,V,complex,Root>(*ctx,operation(s.id),s.n,x.p,y.p,tree.p);
        else if(s.id<2) {
#if defined(KERNEL_BENCH_MPFR_RESIDENT)
            ctx->launch(mplapack_mpfr_inverse::blas1_element_job<V,R>{s.id==0,x.p,y.p,R(1)},s.n);
#else
            if constexpr(!complex) {
                mplapack_blas::gpu_lane_args<Ops> g{};
                g.s.op=s.id==0?mplapack_blas::L1_AXPY:mplapack_blas::L1_COPY;
                g.s.n=s.n;g.lda=g.ldc=1;g.al=alpha.p;
                g.a=reinterpret_cast<F*>(x.p);g.c=reinterpret_cast<F*>(y.p);
                // Same production kernel and64-thread launch policy, without
                // blas_device's post-kernel host synchronization in the interval.
                mplapack_blas::k_blas_lanes<Ops><<<unsigned((s.n+63)/64),64>>>(g,s.n);
                check(cudaGetLastError(),"element kernel");
            }
#endif
        }
#if defined(KERNEL_BENCH_MPFR_RESIDENT)
        else if(s.inverse()) {
            int k=s.inverse_kind();
            if(k==0 || k==2) info=inv->tri(true,false,s.n,s.nbtri,a.p,s.n);
            if(!info && (k==1 || k==2)) inv->lauum(true,s.n,s.nblau,a.p,s.n);
            if(k==3) info=inv->getri(s.n,s.nbtri,s.nb,2,a.p,work.p,piv.p,s.n);
        }
#endif
        check(cudaEventRecord(stop.p,stream),"stop");check(cudaEventSynchronize(stop.p),"stop complete");
        float ms=0;check(cudaEventElapsedTime(&ms,start.p,stop.p),"elapsed");
        check(ctx->error,"launch");
#if defined(KERNEL_BENCH_MPFR_RESIDENT)
        if(s.inverse()) { inv->finish();check(inv->error,"inverse result"); } else {ctx->finish();check(ctx->error,"checked arithmetic");}
#endif
        if(info) throw std::runtime_error("inverse info="+std::to_string(info));
        return double(ms)*1e-3;
    }
    void restore_state() {
#if defined(KERNEL_BENCH_MPFR_RESIDENT)
        check(guard->restore(),"restore arena/precision");
        check(stack.restore(),"restore stack");
#endif
    }
    void result(std::vector<H> &out,REAL &r,COMPLEX &c,long &index) {
        if(s.id>=2 && s.id<=10) {
            S q;tree.download(&q,1);unpack_value(q.re,r);
            if(s.id==5 || s.id==6) {REAL im;unpack_value(q.im,im);c=COMPLEX(r,im);} index=q.index+1;
        } else {
            auto &d=s.inverse()?a:s.solve()?x:y;
            std::vector<V> h(d.count);d.download(h.data(),h.size());out.resize(h.size());
            for(size_t i=0;i<h.size();++i) unpack_value(h[i],out[i]);
        }
    }
};
#endif

// CPU programs link only libmplapack_<backend>_opt, never the CUDA dispatcher.
template<class H> void cpu(shape s,std::vector<H> &a,std::vector<H> &x,std::vector<H> &y,
                           std::vector<H> &work,std::vector<mplapackint> &piv,REAL &r,COMPLEX &c,long &idx) {
    mplapackint info=0;
    if constexpr(std::is_same<H,REAL>::value) {
        switch(s.id) {
        case 0:Raxpy(s.n,REAL(1),x.data(),1,y.data(),1);break;
        case 1:Rcopy(s.n,x.data(),1,y.data(),1);break;
        case 2:idx=iRamax(s.n,x.data(),1);break;
        case 4:r=Rdot(s.n,x.data(),1,y.data(),1);break;
        case 7:r=Rasum(s.n,x.data(),1);break;
        case 9:r=Rnrm2(s.n,x.data(),1);break;
        case 11:Rtrsv("U","N","N",s.n,a.data(),s.n,x.data(),1);break;
        case 12:Rtbsv("U","N","N",s.n,64,a.data(),65,x.data(),1);break;
        case 13:Rtpsv("U","N","N",s.n,a.data(),x.data(),1);break;
        case 17:Rtrtri("U","N",s.n,a.data(),s.n,info);break;
        case 18:Rlauum("U",s.n,a.data(),s.n,info);break;
        case 19:Rpotri("U",s.n,a.data(),s.n,info);break;
        case 20:Rgetri(s.n,a.data(),s.n,piv.data(),work.data(),work.size(),info);break;
        default:throw std::runtime_error("bad real routine");
        }
    } else {
        switch(s.id) {
        case 3:idx=iCamax(s.n,x.data(),1);break;
        case 5:c=Cdotc(s.n,x.data(),1,y.data(),1);break;
        case 6:c=Cdotu(s.n,x.data(),1,y.data(),1);break;
        case 8:r=RCasum(s.n,x.data(),1);break;
        case 10:r=RCnrm2(s.n,x.data(),1);break;
        case 14:Ctrsv("U","N","N",s.n,a.data(),s.n,x.data(),1);break;
        case 15:Ctbsv("U","N","N",s.n,64,a.data(),65,x.data(),1);break;
        case 16:Ctpsv("U","N","N",s.n,a.data(),x.data(),1);break;
        case 21:Ctrtri("U","N",s.n,a.data(),s.n,info);break;
        case 22:Clauum("U",s.n,a.data(),s.n,info);break;
        case 23:Cpotri("U",s.n,a.data(),s.n,info);break;
        case 24:Cgetri(s.n,a.data(),s.n,piv.data(),work.data(),work.size(),info);break;
        default:throw std::runtime_error("bad complex routine");
        }
    }
    if(info) throw std::runtime_error("CPU info="+std::to_string(info));
}

template<class H> void point(shape s,const std::string &mode,int reps,int warmup,const std::string &dumpfile) {
    generator g(20261104+s.n);
    std::vector<H> a(s.acount()),x(s.vector() || s.solve()?s.n:0),y(s.vector()?s.n:0);
    for(auto &v:x) v=draw<H>(g);
    for(auto &v:y) v=draw<H>(g);
    if(s.solve() || s.inverse()) for(long j=0;j<s.n;++j) {
        long first=s.solve() && s.storage()==1?std::max(0L,j-64):0;
        for(long i=first;i<s.n;++i) {
            if(s.solve() && i>j) break;
            H v=draw<H>(g)/REAL(double(s.n));
            if(i==j) v=H(REAL(1)+magnitude(v));
            a[s.solve()?s.index(i,j):i+j*s.n]=v;
        }
    }
    const std::string &name=names[s.id];
    s.nbtri=iMlaenv(1,s.complex()?"Ctrtri":"Rtrtri","UN",s.n,-1,-1,-1);
    s.nblau=iMlaenv(1,s.complex()?"Clauum":"Rlauum","U",s.n,-1,-1,-1);
    s.nb=iMlaenv(1,s.complex()?"Cgetri":"Rgetri"," ",s.n,-1,-1,-1);
    REAL r;COMPLEX c;long idx=0;std::vector<H> out;
    double best=1e300;int runs=0;
#if defined(KERNEL_BENCH_BF_CUDA) || defined(KERNEL_BENCH_MPFR_RESIDENT)
    gpu<H> device(s,a,x,y);
    auto run=[&]() { return device.run(); };
    const char *method="gpu";
#else
    std::vector<H> aa=a,xx=x,yy=y,work(s.inverse()?size_t(s.n)*s.nb:0);
    std::vector<mplapackint> piv(s.inverse()?s.n:0);for(size_t i=0;i<piv.size();++i) piv[i]=i+1;
    auto run=[&]() {
        if(s.inverse()) aa=a; else if(s.solve()) xx=x; else if(s.id<2) yy=y;
        auto t=std::chrono::steady_clock::now();cpu(s,aa,xx,yy,work,piv,r,c,idx);
        return std::chrono::duration<double>(std::chrono::steady_clock::now()-t).count();
    };
    const char *method=mode=="default"?"cpu-default":"cpu";
#endif
    for(int i=0;i<warmup;++i) {if(run()>10) reps=1;}
    for(int i=0;i<reps;++i) { double t=run();best=std::min(best,t);++runs;if(t>10) break; }
    // Download/check is outside the event interval. Dumps allow an independent
    // CPU/GPU comparison of the actual benchmark operands, including low limbs.
#if defined(KERNEL_BENCH_BF_CUDA) || defined(KERNEL_BENCH_MPFR_RESIDENT)
    device.result(out,r,c,idx);
    device.restore_state();
#else
    out=s.inverse()?aa:s.solve()?xx:yy;
#endif
    if(!dumpfile.empty()) {
        FILE *f=std::fopen(dumpfile.c_str(),"w");if(!f) throw std::runtime_error("dump open");
        if(s.id==2 || s.id==3) std::fprintf(f,"%ld\n",idx);
        else if(s.id==5 || s.id==6) { hex(f,c);std::fputc('\n',f); }
        else if(s.reduction()) { hex(f,r);std::fputc('\n',f); }
        else for(const auto &v:out) {hex(f,v);std::fputc('\n',f);}
        std::fclose(f);
    }
    std::printf("%s,%s,%s,%ld,%s,%.9g,%.9g,%d\n",default_label(),name.c_str(),s.solve()?"UNN_k64":s.inverse()?"UN":"inc1",s.n,method,best,s.flops()/best*1e-9,runs);std::fflush(stdout);
}
int main(int argc,char **argv) {
    int id=std::find(names.begin(),names.end(),argv[1])-names.begin();
    const std::string mode=argv[2];if(mode!="fast" && mode!="default") throw std::runtime_error("mode must be fast/default");
#if defined(KERNEL_BENCH_BF_CUDA) || defined(KERNEL_BENCH_MPFR_RESIDENT)
    // Load modules before event intervals, including an explicit --warmup0.
    // Otherwise CUDA lazy loading can make first-launch host latency look like
    // kernel time for the single-measurement large-size protocol.
    setenv("CUDA_MODULE_LOADING","EAGER",1);
    int devices=0;
    if(cudaGetDeviceCount(&devices)!=cudaSuccess || !devices) { std::fprintf(stderr,"no CUDA device\n");return 77; }
    if(mode=="default" && id>=4 && id<=10) throw std::runtime_error("GPU reductions require fast mode");
#endif
#if !defined(MPLAPACK_BUILD_WITH_MPFR)
    if(id>=17) throw std::runtime_error("BF inverses are existing work; use the legacy benchmark");
#endif
    // Only the seven reduction series opt in. Inverse leaves themselves call
    // DOT/ASUM; enabling FAST globally would change their reference arithmetic.
    setenv("MPLAPACK_MPBLAS_FAST",mode=="fast" && id>=4 && id<=10?"1":"0",1);
    setenv("MPLAPACK_LAPACK_FAST","0",1);
    auto sizes=list(argv[3]);if(sizes.empty()) throw std::runtime_error("empty size list");
    int reps=-1,warmup=-1;std::string file;
    for(int i=4;i<argc;++i) {
        std::string arg=argv[i];if(i+1==argc) throw std::runtime_error("missing option value");
        if(arg=="--reps") reps=std::stoi(argv[++i]);else if(arg=="--warmup") warmup=std::stoi(argv[++i]);else if(arg=="--dump") file=argv[++i];else throw std::runtime_error("unknown option "+arg);
    }
#if defined(KERNEL_BENCH_MPFR_CUDA)
    throw std::runtime_error("new MPFR kernels require Rkernel_bench.mpfr_cuda_resident (CUDA events, not CUPTI)");
#endif
#if defined(KERNEL_BENCH_MPFR_RESIDENT)
    if(mpfr_get_default_prec()!=512) throw std::runtime_error("resident performance adapter requires MPFR512");
#endif
    for(long n:sizes) {
        shape s{id,n};if(n<1 || n>(s.vector()?(1L<<24):4096)) throw std::runtime_error("size out of benchmark range");
        int nr=reps<0?(n<=(s.vector()?(1L<<21):2000)?3:1):reps;
        int nw=warmup<0?(n<=(s.vector()?(1L<<21):2000)?1:0):warmup;
#if defined(KERNEL_BENCH_BF_CUDA) || defined(KERNEL_BENCH_MPFR_RESIDENT)
        if(warmup<0) nw=1; // Always prime GPU launches; large sizes still measure once.
#endif
        if(nr<1 || nw<0) throw std::runtime_error("invalid repetitions");
        if(s.complex()) point<COMPLEX>(s,mode,nr,nw,file);else point<REAL>(s,mode,nr,nw,file);
    }
    return 0;
}
}
#endif
