// Included by Rgemm_rt_cuda.cu to reuse its MPFR engine, packs and arena.
#ifndef MPLAPACK_MPFR_FACTOR_RT_BODY_CUH
#define MPLAPACK_MPFR_FACTOR_RT_BODY_CUH
#include "Rfactor_device_body.cuh"
#include "Rsolver_device_body.cuh"
namespace mplapack_mpfr_cuda {
namespace {
struct factor_runtime_engine : l3_rt_engine {
    struct M { dpack P; long ld; int nl; };
    static constexpr bool runtime=true;
    using Update=factor_runtime_engine;
    RT_HD Update update() const { return *this; }
    mutable bool ok=true;
    RT_HD void thread(long tid) { W+=(size_t)tid*8*nl; }
    RT_HD static M shift(M a,long n) { a.P.L+=(size_t)n*a.nl; a.P.E+=n; a.P.K+=n; return a; }
    RT_HD V at(const M &a,long i,long j) const { return l3_rt_engine::at({a.P,a.ld},i,j); }
    RT_HD void put(M &a,long i,long j,const V &v) const { l3_rt_engine::M b{a.P,a.ld}; l3_rt_engine::put(b,i,j,v); }
    RT_HD void checked(V &r) const {
        int k=cu_mpfr_custom_get_kind(r.x);
        if(k==CU_MPFR_ZERO_KIND || k==-CU_MPFR_ZERO_KIND) return;
        if((k!=CU_MPFR_REGULAR_KIND && k!=-CU_MPFR_REGULAR_KIND) ||
           cu_mpfr_custom_get_exp(r.x)<-1073741823L || cu_mpfr_custom_get_exp(r.x)>1073741823L) {
            ok=false; cu_mpfr_set_zero(r.x,1);
        }
    }
    RT_HD void add(V &r,const V &a,const V &b) const { l3_rt_engine::add(r,a,b); checked(r); }
    RT_HD void sub(V &r,const V &a,const V &b) const { l3_rt_engine::sub(r,a,b); checked(r); }
    RT_HD void mul(V &r,const V &a,const V &b) const { l3_rt_engine::mul(r,a,b); checked(r); }
    RT_HD void div(V &r,const V &a,const V &b) const { l3_rt_engine::div(r,a,b); checked(r); }
    RT_HD void absolute(V &r,const V &a) const { arena_reset(); cu_mpfr_abs(r.x,a.x,CU_MPFR_RNDN); }
    RT_HD int compare(const V &a,const V &b) const { return cu_mpfr_cmp(a.x,b.x); }
    RT_HD void root(V &r,const V &a) const { arena_reset(); cu_mpfr_sqrt(r.x,a.x,CU_MPFR_RNDN); checked(r); }
};
#ifndef MPLAPACK_MPFR_CUDA_HOST_EMULATION
struct factor_arena_guard {
    char *base=mpc_cuda_arena_base;
    size_t slab=mpc_cuda_arena_slab;
    size_t *top=mpc_cuda_arena_top;
    factor_arena_guard(char *b,size_t s,size_t *t) {
        mpc_cuda_arena_base=b; mpc_cuda_arena_slab=s; mpc_cuda_arena_top=t;
    }
    ~factor_arena_guard() { mpc_cuda_arena_base=base; mpc_cuda_arena_slab=slab; mpc_cuda_arena_top=top; }
};
#endif
}
int factor_runtime(const factor_shape &s,long prec,rt_pack &a,const rt_pack &safe,long *piv,long &info,factor_times &times) {
    using E=factor_runtime_engine; using Clock=std::chrono::steady_clock;
    int nl=(int)((prec+63)/64);
    size_t count=(size_t)s.m*(s.n+s.nrhs);
    E e; e.prec=prec; e.nl=nl;
#ifdef MPLAPACK_MPFR_CUDA_HOST_EMULATION
    std::vector<cu_mp_limb_t> work(8*nl), limbs(4*nl);
    long exps[4]{}; signed char kinds[4]{};
    dpack mat{(cu_mp_limb_t*)a.limbs,a.exp,a.kind}, c{limbs.data(),exps,kinds};
    e.W=work.data();
#else
    int ndev=0;
    if(cudaGetDeviceCount(&ndev)!=cudaSuccess || !ndev) return 1;
    factor_stack_guard stack(factor_stack_bytes(prec,true));
    if(stack.error) return stack.error;
    const size_t pool=2*factor_pool,slab=16384*((prec+1023)/1024);
    size_t available=0,total=0;
    size_t extra=pool*(slab+sizeof(size_t)+8*nl*sizeof(cu_mp_limb_t))+1048576;
    if(cudaMemGetInfo(&available,&total)!=cudaSuccess || available<extra ||
       count>(available-extra)/(nl*sizeof(cu_mp_limb_t)+sizeof(long)+1)) return 1;
    device_pack matrix,constants;
    factor_buffer<cu_mp_limb_t> work;
    factor_buffer<char> arena; factor_buffer<size_t> top;
    if(!matrix.allocate(count,nl) || !constants.allocate(4,nl) ||
       work.allocate(pool*8*nl) || arena.allocate(pool*slab) || top.allocate(pool)) return 1;
    if(cudaMemset(top.p,0,pool*sizeof(size_t))!=cudaSuccess) return 1;
    factor_arena_guard guard(arena.p,slab,top.p);
    e.W=work.p;
    dpack mat=matrix.d,c=constants.d;
#endif
    // Four constants: one, minus one, sfmin, panel reciprocal. Initialize
    // with the same MPFR engine (no conversion of matrix elements here).
    std::vector<cu_mp_limb_t> cl(4*nl); long ce[4]{}; signed char ck[4]{};
    dpack hc{cl.data(),ce,ck};
#ifdef MPLAPACK_MPFR_CUDA_HOST_EMULATION
    auto v=e.reg(0); e.one(v); store(v.x,nl,hc,0);
    mpfr_neg(v.x,v.x,MPFR_RNDN); store(v.x,nl,hc,1);
#else
    // MPFR layout of +/-1: exponent=1, high significand bit set.
    ck[0]=CU_MPFR_REGULAR_KIND; ck[1]=-CU_MPFR_REGULAR_KIND;
    ce[0]=ce[1]=1; cl[nl-1]=cl[2*nl-1]=cu_mp_limb_t(1)<<63;
#endif
    hc.K[2]=safe.kind[0]; hc.E[2]=safe.exp[0]; hc.K[3]=CU_MPFR_ZERO_KIND;
    for(int i=0;i<nl;++i) hc.L[2*nl+i]=safe.limbs[i];
#ifdef MPLAPACK_MPFR_CUDA_HOST_EMULATION
    std::copy_n(hc.L,4*nl,c.L); std::copy_n(hc.E,4,c.E); std::copy_n(hc.K,4,c.K);
#else
    rt_pack packed{(unsigned long long*)hc.L,hc.E,hc.K};
    if(!constants.upload(packed,4,nl)) return 1;
#endif
    factor_context<E> ctx(e,E::M{c,4,nl},s.m,s.n);
    if(ctx.error) return ctx.error;
    auto start=Clock::now();
#ifndef MPLAPACK_MPFR_CUDA_HOST_EMULATION
    if(!matrix.upload(a,count,nl) || cudaStreamSynchronize(nullptr)!=cudaSuccess) return 1;
#endif
    auto uploaded=Clock::now();
    long result=s.nrhs?factor_and_solve(ctx,s,E::M{mat,s.m,nl},piv):ctx.factor(s,E::M{mat,s.m,nl});
    if(ctx.error) return ctx.error;
    auto factored=Clock::now();
#ifndef MPLAPACK_MPFR_CUDA_HOST_EMULATION
    if(s.nrhs && !s.factorize) {
        size_t offset=size_t(s.m)*s.n, nc=size_t(s.m)*s.nrhs;
        if(solver_download(a.limbs,reinterpret_cast<const unsigned long long*>(mat.L),nc*nl,offset*nl) ||
           solver_download(a.exp,mat.E,nc,offset) || solver_download(a.kind,mat.K,nc,offset)) return 1;
    } else if(!matrix.download(a,count,nl)) return 1;
#endif
    if(s.factorize && !s.chol && ctx.pivots.download(piv,std::min(s.m,s.n))) return 1;
    auto downloaded=Clock::now();
#ifndef MPLAPACK_MPFR_CUDA_HOST_EMULATION
    if(stack.restore()) return 1;
#endif
    info=result;
    times.upload=std::chrono::duration<double>(uploaded-start).count();
    times.factor=std::chrono::duration<double>(factored-uploaded).count();
    times.download=std::chrono::duration<double>(downloaded-factored).count();
    return 0;
}
}
#endif
