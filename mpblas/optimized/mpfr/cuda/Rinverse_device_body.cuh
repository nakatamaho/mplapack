// MPFR residency/launch adapter for the common inverse reference schedule.
#ifndef MPLAPACK_MPFR_INVERSE_DEVICE_BODY_CUH
#define MPLAPACK_MPFR_INVERSE_DEVICE_BODY_CUH
#include "Rinverse_kernel_cuda.h"
#include "Cfactor_device_body.cuh"
#include <vector>
#include "../../../../mplapack/optimized/common/Rinverse_schedule.h"
namespace mplapack_mpfr_inverse {
using namespace mplapack_mpfr_cuda;
namespace {
// Generated inverse leaves add TRMV/TRMM frames above the factor arithmetic.
// Reserve a conservative extra factor-policy frame; restore on every exit.
size_t inverse_stack_bytes(long p,bool fixed) {
    return 2*(fixed?complex_fixed_stack_bytes(p):complex_runtime_stack_bytes(p));
}
struct inverse_engine {
    struct M {};
    static constexpr bool runtime=true;
    bool ok=true;
    L3_HD void thread(long) { failure()=false; factor_arena_reset(); }
};
template<class J> struct checked_job {
    inverse_engine e; J j;
    L3_HD void operator()(long lane) { j(lane); if(failure()) e.ok=false; }
};
// Block-cooperative checked wrapper for the shared-memory chunk reductions
// (the chunk jobs carry no engine of their own).
template<class J> struct block_checked_job {
    inverse_engine e; J j;
    L3_HD void block_run(int tid,int T,long u,char *sh) { j.block_run(tid,T,u,sh); if(failure()) e.ok=false; }
};
// Goal 13 wavefront: the whole triangular schedule runs in one block with
// per-element flags (see Rtrisolve_jobs.h).  The runtime-precision engine
// keeps one live value per owner thread — the same residency the flat
// per-column schedule uses.
template<class J> struct wave_checked_job {
    inverse_engine e; J j;
    L3_HD void block_run(int tid,int T,int *flags,int *prog) { e.thread(tid); j.block_run(tid,T,flags,prog); if(failure()) e.ok=false; }
};
#ifndef MPLAPACK_MPFR_CUDA_HOST_EMULATION
template<class J> __global__ void k_mpfr_wave(wave_checked_job<J> job,int *failed,int *buf) {
    const long gtid=(long)blockIdx.x*blockDim.x+threadIdx.x;
    const long gT=(long)gridDim.x*blockDim.x;
    job.block_run(gtid,gT,buf,buf+job.j.s.n);
    if(!job.e.ok) atomicExch(failed,1);
}
#endif
struct launcher : factor_context<inverse_engine> {
    // the checked arithmetic keeps the lane-parallel sub-lane schedule
    static constexpr bool chunked_reductions=false;
    launcher():factor_context<inverse_engine>(inverse_engine{},inverse_engine::M{},1,1) {}
    long *status_data() { return status.p; }
    long feedback(long *) { return factor_context<inverse_engine>::feedback(); }
    template<class J> void launch(J j,long lanes) { this->execute(checked_job<J>{inverse_engine{},j},lanes); }
    template<class J> void launch_wave(J j) {
        if(error) return;
        wave_checked_job<J> w{inverse_engine{},j};
#ifndef MPLAPACK_MPFR_CUDA_HOST_EMULATION
        const long n=j.s.n;
        if(n<=0) return;
        const int T=256;
        int perSM=0,sms=0;
        if(cudaOccupancyMaxActiveBlocksPerMultiprocessor(&perSM,(const void *)k_mpfr_wave<J>,T,0)!=cudaSuccess || perSM<=0)
            perSM=1;
        if(cudaDeviceGetAttribute(&sms,cudaDevAttrMultiProcessorCount,0)!=cudaSuccess || sms<=0)
            sms=1;
        // The checked arithmetic arena is indexed by the physical thread id
        // (factor_pool slots); the grid must not exceed it.
        long cap=(long)perSM*sms;
        const long poolcap=mplapack_mpfr_cuda::factor_pool/T;
        if(cap>poolcap) cap=poolcap;
        long blocks=cap;
        factor_buffer<int> fb;
        if(fb.allocate(size_t(n)*2) || cudaMemset(fb.p,0,size_t(n)*2*sizeof(int))) { error=1;return; }
        wave_checked_job<J> w2{inverse_engine{},j};
        int *p=fb.p,*fp=failed.p;
        void *args[3]={&w2,&fp,&p};
        error=cudaLaunchCooperativeKernel((const void *)k_mpfr_wave<J>,dim3(unsigned(blocks)),dim3(T),args,0,0);
        if(!error) error=cudaDeviceSynchronize();
#else
        std::vector<int> f(size_t(j.s.n>0?j.s.n:0)*2,0);
        w.block_run(0,1,f.data(),f.data()+(j.s.n>0?j.s.n:0));
        if(!w.e.ok) *failed.p=1;
#endif
    }
    // Goal 13 (continued) fast reductions: one chunk of the fixed two-level
    // order per block iteration, folded in shared memory (see
    // Rblas1_jobs_common.h); 32 threads = one warp = the sub-lane workers.
    // The grid stays within the pool: the arena is indexed by the physical
    // thread id.  The engine wraps the chunk job like checked_job wraps the
    // lane jobs.
    template<class J> void launch_chunked(J job,long units) {
        if(error) return;
        block_checked_job<J> w{inverse_engine{},job};
        size_t full=J::smem_bytes();
        if(full<=160*1024) {
            error=execute_blocks(w,units,32,full,true);
            return;
        }
        // chunk elements too large for the shared memory: unstaged variant
        // (x/y from global memory, same order)
        size_t sm=J::sublane_smem_bytes();
        if(sm>160*1024) { error=1;return; }
        error=execute_blocks(w,units,32,sm,false);
    }
    // Returns 0 when the one-kernel shared-memory tree was launched; nonzero
    // asks the schedule for the per-stage global-memory fallback.
    template<class J> int launch_chunk_tree(J job,long units) {
        if(error) return 1;
        const size_t need=J::smem_bytes(units);
        if(need>160*1024) return 1; // the chunk states do not fit
        block_checked_job<J> w{inverse_engine{},job};
        return execute_blocks(w,1,256,need,true);
    }
    void finish() {
        factor_context<inverse_engine>::finish(); int f=0;
        if(!error) error=failed.download(&f,1);
        if(!error && f) error=1;
    }
};
}
template<class Tag,bool Complex> int run(const shape &s,long prec,value<Tag,Complex> *a,value<Tag,Complex> *work,
    const long *piv,long &info,factor_times &times) {
    using V=value<Tag,Complex>;
    using A=mplapack_inverse_gpu::inverse_arithmetic<value<Tag,false>,value<Tag,true>,Complex>;
    using Clock=std::chrono::steady_clock;
    size_t count=size_t(s.n)*s.n,wcount=s.kind==3?size_t(s.n)*s.columns:0;
#ifndef MPLAPACK_MPFR_CUDA_HOST_EMULATION
    int devices=0; if(cudaGetDeviceCount(&devices) || !devices) return 1;
    // Leaves call generated BLAS and MPC division; use the complex policy for
    // both kinds, including cumulative caller frames, not the BLAS minimum.
    constexpr bool fixed=std::is_same<Tag,factor_fixed_tag<256>>::value || std::is_same<Tag,factor_fixed_tag<512>>::value ||
        std::is_same<Tag,factor_fixed_tag<768>>::value || std::is_same<Tag,factor_fixed_tag<1024>>::value || std::is_same<Tag,factor_fixed_tag<2048>>::value;
    mpfr_stack_guard stack(inverse_stack_bytes(prec,fixed));
    if(stack.error) return stack.error;
    factor_buffer<char> arena; factor_buffer<size_t> tops;
    size_t slab=16384*((prec+1023)/1024);
    if(arena.allocate(factor_pool*slab) || tops.allocate(factor_pool) || cudaMemset(tops.p,0,factor_pool*sizeof(size_t))) return 1;
    complex_arena_guard guard(arena.p,slab,tops.p,prec); if(guard.error) return guard.error;
#else
    long previous=mplapack_mpfr_rt_prec_h; mplapack_mpfr_rt_prec_h=prec;
    struct precision_guard { long p; ~precision_guard() { mplapack_mpfr_rt_prec_h=p; } } guard{previous};
#endif
    factor_buffer<V> matrix,workspace; factor_buffer<long> pivots;
    if(matrix.allocate(count) || (wcount && (workspace.allocate(wcount) || pivots.allocate(s.n)))) return 1;
    mplapack_inverse_gpu::inverse_schedule<A,launcher> ctx;
    if(ctx.error) return ctx.error;
    auto start=Clock::now();
    if(matrix.upload(a,count) || (wcount && (workspace.upload(work,wcount) || pivots.upload(piv,s.n)))) return 1;
#ifndef MPLAPACK_MPFR_CUDA_HOST_EMULATION
    if(cudaStreamSynchronize(nullptr)) return 1;
#endif
    auto uploaded=Clock::now(); long result=0;
    if(s.kind==0 || s.kind==2) result=ctx.tri(s.upper,s.unit,s.n,s.nbtri,matrix.p,s.n);
    if(!result && (s.kind==1 || s.kind==2)) ctx.lauum(s.upper,s.n,s.nblau,matrix.p,s.n);
    if(s.kind==3) result=ctx.getri(s.n,s.nbtri,s.nb,s.nbmin,matrix.p,workspace.p,pivots.p,s.n);
    ctx.finish(); if(ctx.error) return ctx.error;
    auto computed=Clock::now();
    if(matrix.download(a,count) || (wcount && !result && workspace.download(work,wcount))) return 1;
#ifndef MPLAPACK_MPFR_CUDA_HOST_EMULATION
    if(guard.restore() || stack.restore()) return 1;
#endif
    info=result;
    times.upload=std::chrono::duration<double>(uploaded-start).count();
    times.factor=std::chrono::duration<double>(computed-uploaded).count();
    times.download=std::chrono::duration<double>(Clock::now()-computed).count();
    return 0;
}
}
#endif
