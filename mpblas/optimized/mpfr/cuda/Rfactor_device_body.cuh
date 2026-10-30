// MPFR-only reference-order real LU/Cholesky scheduler and panel jobs.
// Included after the arithmetic engine; never includes host mpfr_class.
#ifndef MPLAPACK_MPFR_FACTOR_DEVICE_BODY_CUH
#define MPLAPACK_MPFR_FACTOR_DEVICE_BODY_CUH
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <vector>
#include <limits>
#include "Rfactor_kernel_cuda.h"
#include "blas3_kernel_cuda.h"
#include "mpfr_cuda_stack.h"
#ifndef MPLAPACK_MPFR_CUDA_HOST_EMULATION
#include <cuda_runtime.h>
#endif
namespace mplapack_mpfr_cuda {
namespace {
constexpr long factor_threads=64,factor_blocks=128,factor_pool=factor_threads*factor_blocks;
template<class E> struct factor_job {
    E e;
    typename E::M a{}, b{}, c{}, constants{};
    l3_shape s{};
    gemm_shape g{};
    long *piv=nullptr, *status=nullptr, *candidates=nullptr;
    long m=0,n=0,k1=0,k2=0,offset=0;
    int op=0;
    // 0 BLAS, 1 GEMM, 2 pivot chunks, 3 pivot/reciprocal, 4 scale,
    // 5 swaps, 6 adjust pivots, 7 sqrt. Pivots are local, one-based.
    L3_HD void operator()(long lane) {
        if(op==0) {
            blas_mats<E> x{};
            x.a=a; x.b=b; x.c=c;
            x.al=constants; x.be=constants;
#ifdef MPLAPACK_MPFR_FACTOR_COMPLEX
            e.trsm(s,x,lane);
#else
            l3_trsm(e,s,x,lane);
#endif
        } else if(op==2) {
            long start=lane*128, end=start+128<m?start+128:m, best=start;
            auto largest=e.reg(0), v=e.reg(1);
            e.absolute(largest,e.at(a,start,0));
            for(long i=start+1;i<end;++i) {
                e.absolute(v,e.at(a,i,0));
                if(e.compare(v,largest)>0) { best=i; e.set(largest,v); }
            }
            candidates[lane]=best;
        } else if(op==3) {
            long best=0;
            if(m>1) {
                auto largest=e.reg(0),v=e.reg(1);
                best=candidates[0]; e.absolute(largest,e.at(a,best,0));
                for(long i=1;i<(m+127)/128;++i) {
                    long p=candidates[i]; e.absolute(v,e.at(a,p,0));
                    if(e.compare(v,largest)>0) { best=p; e.set(largest,v); }
                }
            }
            piv[0]=best+1;
            status[0]=e.iszero(e.at(a,best,0))?1:0;
            if(m>1 && !status[0]) {
                auto v=e.reg(0), t=e.reg(1);
                e.set(t,e.at(a,0,0)); e.set(v,e.at(a,best,0));
                e.put(a,0,0,v); e.put(a,best,0,t);
#ifdef MPLAPACK_MPFR_FACTOR_COMPLEX
                e.magnitude(v,e.at(a,0,0));
#else
                e.absolute(v,e.at(a,0,0));
#endif
                status[1]=e.compare(v,e.at(constants,2,0))>=0;
                if(status[1]) {
                    e.div(t,e.at(constants,0,0),e.at(a,0,0));
                    e.put(constants,3,0,t);
                }
            }
        } else if(op==4) {
            auto t=e.reg(0);
            if(status[1]) e.mul(t,e.at(constants,3,0),e.at(a,lane+1,0));
            else e.div(t,e.at(a,lane+1,0),e.at(a,0,0));
            e.put(a,lane+1,0,t);
        } else if(op==5) {
            auto t=e.reg(0),v=e.reg(1);
            for(long i=k1;i<k2;++i) {
                long p=piv[i]-1;
                if(p==i) continue;
                e.set(t,e.at(a,i,lane)); e.set(v,e.at(a,p,lane));
                e.put(a,i,lane,v); e.put(a,p,lane,t);
            }
        } else if(op==6) piv[lane]+=offset;
        else if(op==7) {
            auto t=e.reg(0),z=e.reg(1); e.zero(z);
            status[0]=e.compare(e.at(a,0,0),z)<=0?1:0;
            if(!status[0]) { e.root(t,e.at(a,0,0)); e.put(a,0,0,t); }
        }
    }
};
// Inlined update arithmetic, separate from the panel's noinline Newton and
// BLAS job. This reuses the GEMM/SYRK element bodies without routing every
// multiply/add through a function call or compiling all panel branches into
// the high-throughput update kernel.
template<class E> struct factor_update_job {
    E e; typename E::M a{},b{},c{},constants{};
    l3_shape s{}; gemm_shape g{}; int op=0;
    L3_HD void operator()(long lane) {
        if(op==0) {
            blas_mats<E> x{}; x.a=a; x.c=c;
            x.al=E::shift(constants,1); x.be=constants;
#ifdef MPLAPACK_MPFR_FACTOR_COMPLEX
            e.herk(s,x,lane);
#else
            l3_syrk(e,s,x,lane%s.n,lane/s.n);
#endif
        } else {
            auto al=E::shift(constants,1),ab=a;
#ifdef MPLAPACK_MPFR_FACTOR_COMPLEX
            e.gemm(g,al,constants,a,b,c,lane%g.m,lane/g.m);
#else
            gemm_engine_element(e,g,al,constants,a,b,ab,c,lane%g.m,lane/g.m,false);
#endif
        }
    }
};
#ifndef MPLAPACK_MPFR_CUDA_HOST_EMULATION
template<class Job> __global__ void factor_kernel(Job job,long units,int *failed,long offset) {
    long tid=(long)blockIdx.x*blockDim.x+threadIdx.x;
    if(tid<offset) return;
    // libmpc_cuda indexes its arena by the physical CUDA thread index. Far
    // updates occupy a disjoint index range, without changing arena globals
    // while either stream is running. The skipped blocks do no arithmetic.
    job.e.thread(tid);
    for(long lane=tid-offset;lane<units;lane+=(long)gridDim.x*blockDim.x-offset) job(lane);
    if(!job.e.ok) atomicExch(failed,1);
}
using factor_stack_guard=mpfr_stack_guard;
template<class T> struct factor_buffer {
    T *p=nullptr;
    ~factor_buffer() { cudaFree(p); }
    int allocate(size_t n) { return cudaMalloc(&p,std::max(size_t(1),n)*sizeof(T)); }
    int upload(const T *h,size_t n) { return cudaMemcpy(p,h,n*sizeof(T),cudaMemcpyHostToDevice); }
    int download(T *h,size_t n) { return cudaMemcpy(h,p,n*sizeof(T),cudaMemcpyDeviceToHost); }
};
#else
template<class T> struct factor_buffer {
    std::vector<T> data; T *p=nullptr;
    int allocate(size_t n) { data.resize(std::max(size_t(1),n)); p=data.data(); return 0; }
    int upload(const T *h,size_t n) { std::copy_n(h,n,p); return 0; }
    int download(T *h,size_t n) { std::copy_n(p,n,h); return 0; }
};
#endif
template<class E> class factor_context {
  public:
    using M=typename E::M;
    E engine;
    M constants;
    factor_buffer<long> status,candidates,pivots;
    factor_buffer<int> failed;
    int error=0; bool far=false;
#ifndef MPLAPACK_MPFR_CUDA_HOST_EMULATION
    cudaStream_t streams[2]{}; cudaEvent_t ready{},done{};
    cudaStream_t stream() const { return streams[far?1:0]; }
#endif
    factor_context(E e,M c,long m,long n):engine(e),constants(c) {
        error=status.allocate(2);
        if(!error) error=candidates.allocate((m+127)/128);
        if(!error) error=pivots.allocate(std::min(m,n));
        if(!error) error=failed.allocate(1);
        int z=0; if(!error) error=failed.upload(&z,1);
#ifndef MPLAPACK_MPFR_CUDA_HOST_EMULATION
        if(!error) error=cudaStreamCreateWithFlags(&streams[0],cudaStreamNonBlocking);
        if(!error) error=cudaStreamCreateWithFlags(&streams[1],cudaStreamNonBlocking);
        if(!error) error=cudaEventCreateWithFlags(&ready,cudaEventDisableTiming);
        if(!error) error=cudaEventCreateWithFlags(&done,cudaEventDisableTiming);
#endif
    }
    ~factor_context() {
        finish();
#ifndef MPLAPACK_MPFR_CUDA_HOST_EMULATION
        if(ready) cudaEventDestroy(ready); if(done) cudaEventDestroy(done);
        for(auto s:streams) if(s) cudaStreamDestroy(s);
#endif
    }
    void finish() {
#ifndef MPLAPACK_MPFR_CUDA_HOST_EMULATION
        for(auto s:streams) if(s) { int r=cudaStreamSynchronize(s); if(!error) error=r; }
#endif
    }
    void join() {
#ifndef MPLAPACK_MPFR_CUDA_HOST_EMULATION
        if(!error) error=cudaEventRecord(done,streams[1]);
        if(!error) error=cudaStreamWaitEvent(streams[0],done,0);
#endif
    }
    void release() {
#ifndef MPLAPACK_MPFR_CUDA_HOST_EMULATION
        if(!error) error=cudaEventRecord(ready,streams[0]);
        if(!error) error=cudaStreamWaitEvent(streams[1],ready,0);
#endif
    }
    factor_job<E> job(int op,M a) {
        factor_job<E> j{}; j.e=engine; j.op=op; j.a=a; j.constants=constants;
        j.status=status.p; j.candidates=candidates.p; return j;
    }
    template<class Job> void execute(Job j,long units) {
        if(error || !units) return;
#ifndef MPLAPACK_MPFR_CUDA_HOST_EMULATION
        long blocks=(units+factor_threads-1)/factor_threads;
        // Only cu_mpfr needs a bounded pool (its per-thread arena/registers).
        // Fixed updates expose every output lane to the GPU, as GEMM does.
        if(E::runtime) blocks=std::min(factor_blocks,blocks);
        long offset=E::runtime && far?factor_pool:0;
        factor_kernel<<<(unsigned)(blocks+offset/factor_threads),factor_threads,0,stream()>>>(j,units,failed.p,offset);
        error=cudaGetLastError();
#else
        j.e.thread(0);
        for(long i=0;i<units;++i) j(i);
        if(!j.e.ok) *failed.p=1;
#endif
    }
    void launch(factor_job<E> j,long units) {
        if(j.op==1 || (j.op==0 && j.s.op==L3_SYRK)) {
            factor_update_job<typename E::Update> u{};
            u.e=j.e.update(); u.a=j.a; u.b=j.b; u.c=j.c; u.constants=j.constants;
            u.s=j.s; u.g=j.g; u.op=j.op; execute(u,units);
        } else execute(j,units);
    }
    long feedback() {
        long result=0;
#ifndef MPLAPACK_MPFR_CUDA_HOST_EMULATION
        if(!error) error=cudaMemcpyAsync(&result,status.p,sizeof(long),cudaMemcpyDeviceToHost,stream());
        if(!error) error=cudaStreamSynchronize(stream());
#else
        result=status.p[0];
#endif
        return result;
    }
    void blas(l3_shape s,M a,M c) {
        if(!s.m && s.op==L3_TRSM) return;
        if(!s.n || (s.op==L3_SYRK && !s.k)) return;
        auto j=job(0,a); j.s=s; j.c=c; launch(j,blas_lanes(s));
    }
    void solve(bool upper,long m,long n,M a,M b,bool chol) {
        l3_shape s{}; s.op=L3_TRSM; s.m=m; s.n=n;
        s.left=chol?upper:1; s.upper=chol?upper:0;
        s.trans=chol; s.nounit=chol; s.alpha_one=1;
        s.lda=a.ld; s.ldc=b.ld; blas(s,a,b);
    }
    void rank(bool upper,long n,long k,M a,M c) {
        l3_shape s{}; s.op=L3_SYRK; s.n=n; s.k=k;
        s.upper=upper; s.trans=upper; s.beta_one=1;
        s.lda=a.ld; s.ldc=c.ld; blas(s,a,c);
    }
    void update(bool upper,long m,long n,long k,M a,M b,M c) {
        if(!m || !n || !k) return;
        auto j=job(1,a); j.b=b; j.c=c;
        j.g={ !upper,upper,m,n,k,a.ld,b.ld,c.ld,0,1 };
        launch(j,m*n);
    }
    void swap(long n,M a,long k1,long k2,long *piv) {
        auto j=job(5,a); j.k1=k1; j.k2=k2; j.piv=piv; launch(j,n);
    }
    void adjust(long *piv,long n,long offset) {
        auto j=job(6,constants); j.piv=piv; j.offset=offset; launch(j,n);
    }
    long lu_recursive(long m,long n,M a,long *piv) {
        if(error || !m || !n) return 0;
        if(m==1 || n==1) {
            auto j=job(2,a); j.m=m; j.piv=piv;
            if(m>1) launch(j,(m+127)/128);
            j.op=3; launch(j,1);
            long info=feedback();
            if(m>1 && !info) { j.op=4; launch(j,m-1); }
            return info;
        }
        long n1=std::min(m,n)/2,n2=n-n1;
        long info=lu_recursive(m,n1,a,piv);
        M right=E::shift(a,n1*a.ld), lower=E::shift(a,n1), diag=E::shift(right,n1);
        swap(n2,right,0,n1,piv);
        solve(false,n1,n2,a,right,false);
        // LU uses N,N; Cholesky lower uses N,T.
        auto j=job(1,lower); j.b=right; j.c=diag;
        j.g={1,1,m-n1,n2,n1,a.ld,a.ld,a.ld,0,1}; launch(j,(m-n1)*n2);
        long i=lu_recursive(m-n1,n2,diag,piv+n1);
        if(!info && i) info=i+n1;
        adjust(piv+n1,std::min(m,n)-n1,n1);
        swap(n1,a,n1,std::min(m,n),piv);
        return info;
    }
    long lu(long m,long n,long nb,M a) {
        long mn=std::min(m,n),info=0;
        if(nb<=1 || nb>=mn) return lu_recursive(m,n,a,pivots.p);
        for(long col=0;col<mn && !error;col+=nb) {
            long jb=std::min(nb,mn-col);
            long i=lu_recursive(m-col,jb,E::shift(a,col+col*a.ld),pivots.p+col);
            if(!info && i) info=col+i;
            join(); adjust(pivots.p+col,jb,col);
            swap(col,a,col,col+jb,pivots.p);
            for(int part=0;part<2;++part) {
                long c=col+jb+(part?std::min(nb,n-col-jb):0);
                long cols=part?n-c:std::min(nb,n-c); far=part;
                if(cols) {
                    swap(cols,E::shift(a,c*a.ld),col,col+jb,pivots.p);
                    M right=E::shift(a,col+c*a.ld);
                    solve(false,jb,cols,E::shift(a,col+col*a.ld),right,false);
                    auto j=job(1,E::shift(a,col+jb+col*a.ld)); j.b=right; j.c=E::shift(right,jb);
                    j.g={1,1,m-col-jb,cols,jb,a.ld,a.ld,a.ld,0,1}; launch(j,(m-col-jb)*cols);
                }
                if(!part) release();
            }
            far=false;
        }
        return info;
    }
    long chol_recursive(bool upper,long n,M a) {
        if(error || !n) return 0;
        if(n==1) { launch(job(7,a),1); return feedback(); }
        long n1=n/2,n2=n-n1;
        long info=chol_recursive(upper,n1,a);
        if(error || info) return info;
        M off=E::shift(a,upper?n1*a.ld:n1), diag=E::shift(a,n1+n1*a.ld);
        solve(upper,upper?n1:n2,upper?n2:n1,a,off,true);
        rank(upper,n2,n1,off,diag);
        info=chol_recursive(upper,n2,diag);
        return info?info+n1:0;
    }
    long chol(bool upper,long n,long nb,M a) {
        if(nb<=1 || nb>=n) return chol_recursive(upper,n,a);
        for(long col=0;col<n && !error;col+=nb) {
            long jb=std::min(nb,n-col);
            M diag=E::shift(a,col+col*a.ld);
            rank(upper,jb,col,E::shift(a,upper?col*a.ld:col),diag);
            long info=chol_recursive(upper,jb,diag);
            if(error || info) { finish(); return info?info+col:0; }
            join();
            for(int part=0;part<2;++part) {
                long c=col+jb+(part?std::min(nb,n-col-jb):0);
                long rest=part?n-c:std::min(nb,n-c); far=part;
                if(rest) {
                    M off=E::shift(a,upper?col+c*a.ld:c+col*a.ld);
                    update(upper,upper?jb:rest,upper?rest:jb,col,
                           E::shift(a,upper?col*a.ld:c),E::shift(a,upper?c*a.ld:col),off);
                    solve(upper,upper?jb:rest,upper?rest:jb,diag,off,true);
                }
                if(!part) release();
            }
            far=false;
        }
        return 0;
    }
    long factor(const factor_shape &s,M a) {
        long info=s.chol?chol(s.upper,s.n,s.nb,a):lu(s.m,s.n,s.nb,a);
        finish(); int f=0;
        if(!error) error=failed.download(&f,1);
        if(!error && f) error=1;
        return info;
    }
};
}
}
#endif
