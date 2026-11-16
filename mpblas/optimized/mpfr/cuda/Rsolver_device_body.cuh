// Continuation on the factor context; no public BLAS calls (mutex already held).
#ifndef MPLAPACK_MPFR_SOLVER_DEVICE_BODY_CUH
#define MPLAPACK_MPFR_SOLVER_DEVICE_BODY_CUH
namespace mplapack_mpfr_cuda {
namespace {
template<class V> static int solver_download(V *host,const V *device,size_t count,size_t offset) {
#ifndef MPLAPACK_MPFR_CUDA_HOST_EMULATION
    return cudaMemcpy(host+offset,device+offset,count*sizeof(V),cudaMemcpyDeviceToHost);
#else
    std::copy_n(device+offset,count,host+offset); return 0;
#endif
}
#ifndef MPLAPACK_MPFR_FACTOR_COMPLEX
// The panel engine deliberately outlines add/sub/mul to bound panel frames.
// Whole-system TRSM must instead use its existing inline update/BLAS policy:
// identical signed, checked arithmetic, without a device call at every term.
template<class E> struct solver_real_trsm_job {
    E e; typename E::M a{},b{},constants{}; l3_shape s{};
    L3_HD void operator()(long lane) {
        blas_mats<E> x{}; x.a=a; x.c=b; x.al=constants; x.be=constants;
        l3_trsm(e,s,x,lane);
    }
};
#endif
template<class E> struct solver_swaps {
    E e; typename E::M a; const long *piv; long n; bool reverse;
    L3_HD void operator()(long j) {
        auto t=e.reg(0),v=e.reg(1);
        for(long k=0;k<n;++k) {
            long i=reverse?n-1-k:k,p=piv[i]-1;
            if(i==p) continue;
            e.set(t,e.at(a,i,j)); e.set(v,e.at(a,p,j));
            e.put(a,i,j,v); e.put(a,p,j,t);
        }
    }
};
template<class E> static long factor_and_solve(factor_context<E> &ctx,const factor_shape &s,typename E::M a,long *piv) {
    if(!s.factorize && !s.chol && !ctx.error) {
        ctx.error=ctx.pivots.upload(piv,s.n);
#ifndef MPLAPACK_MPFR_CUDA_HOST_EMULATION
        // Pageable H2D may return after staging, before DMA completes. The
        // nonblocking solve stream must not read pivots ahead of that copy.
        if(!ctx.error) ctx.error=cudaStreamSynchronize(nullptr);
#endif
    }
    long result=s.factorize?ctx.factor(s,a):0;
    if(s.nrhs && !result && !ctx.error) {
        auto b=E::shift(a,s.m*s.n);
        auto swap=[&](bool reverse) { ctx.execute(solver_swaps<E>{ctx.engine,b,ctx.pivots.p,s.n,reverse},s.nrhs); };
        auto trsm=[&](bool upper,int trans,bool unit) {
            l3_shape t{}; t.op=L3_TRSM; t.left=1; t.upper=upper; t.trans=trans;
            t.nounit=!unit; t.alpha_one=1; t.m=s.n; t.n=s.nrhs; t.lda=a.ld; t.ldc=b.ld;
#ifdef MPLAPACK_MPFR_FACTOR_COMPLEX
            ctx.blas(t,a,b); // Keep the already measured MPC policy unchanged.
#else
            solver_real_trsm_job<typename E::Update> job{};
            job.e=ctx.engine.update(); job.a=a; job.b=b; job.constants=ctx.constants; job.s=t;
            ctx.execute(job,blas_lanes(t));
#endif
        };
        if(s.chol) { trsm(s.upper,s.upper?1:0,false); trsm(s.upper,s.upper?0:1,false); }
        else if(s.trans=='N') { swap(false); trsm(false,0,true); trsm(true,0,false); }
        else { int t=s.trans=='T'?2:1; trsm(true,t,false); trsm(false,t,true); swap(true); }
        ctx.finish(); int failed=0;
        if(!ctx.error) ctx.error=ctx.failed.download(&failed,1);
        if(!ctx.error && failed) ctx.error=1;
    }
    return result;
}
}
}
#endif
