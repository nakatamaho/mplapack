// Included inside mplapack_mpfr_inverse in Rinverse_dynamic_cuda.cu.
// Dynamic real BLAS1 uses the same owning scalar pool and packed views.
namespace {
MPLAPACK_DEV_HD bool operator<(const dynamic_real &a,const dynamic_real &b) {
#ifdef MPLAPACK_RT_ACTIVE
    MPLAPACK_RT_DECL(x);MPLAPACK_RT_DECL(y);a.view(x);b.view(y);
    return MPLAPACK_RT(mpfr_cmp)(x,y)<0;
#else
    return false;
#endif
}
MPLAPACK_DEV_HD bool operator>(const dynamic_real &a,const dynamic_real &b) { return b<a; }
MPLAPACK_DEV_HD dynamic_real abs(const dynamic_real &a) { return a.kind<0?-a:a; }
MPLAPACK_MPFR_DEV_OP dynamic_real blas1_dynamic_sqrt(const dynamic_real &a) {
    dynamic_real r;
#ifdef MPLAPACK_RT_ACTIVE
    if(!a.safe_operand() || !r.storage()) { failure()=true;return r; }
    MPLAPACK_RT_DECL(x);MPLAPACK_RT_DECL(z);a.view(x);r.view(z);factor_arena_reset();
    MPLAPACK_RT(mpfr_sqrt)(z,x,MPLAPACK_RTU(MPFR_RNDN));MPLAPACK_RT_TAKE(r,z);
    if(!r.valid()) failure()=true;
#endif
    return r;
}
struct dynamic_blas1_job {
    dynamic_engine e;
    mplapack_blas1::operation op;long n,stride;int phase;
    dynamic_real *a;long *indices;
    L3_HD void operator()(long lane) {
        using namespace mplapack_blas1;
        {
            dynamic_real *x=a,*y=a+n,*re=a+2*n,*scale=a+3*n;
            if(phase==0) {
                auto s=leaf<dynamic_real>(op,x[lane],op==dot?y[lane]:dynamic_real(0.0),lane);
                re[lane]=s.re;scale[lane]=s.scale;indices[lane]=lane;
            } else if(phase==1) {
                long i=lane*stride*2;
                if(i<n-stride) {
                    state<dynamic_real> left,right;
                    left.re=re[i];left.scale=scale[i];left.index=indices[i];
                    right.re=re[i+stride];right.scale=scale[i+stride];right.index=indices[i+stride];
                    auto s=merge(op,left,right);
                    re[i]=s.re;scale[i]=s.scale;indices[i]=s.index;
                }
            } else if(phase==2) {
                if(op==mplapack_blas1::norm) re[0]=scale[0]*blas1_dynamic_sqrt(re[0]);
            } else {
                if(op==axpy) y[lane]=y[lane]+a[4*n]*x[lane];else y[lane]=x[lane];
            }
        }
        if(slots()!=0) failure()=true;
        if(failure()) e.ok=false;
    }
};
}
bool dynamic_blas1_memory_fits(long n,long p) { return dynamic_storage_fits(4*size_t(n)+1,n,p); }
int run_dynamic_blas1(mplapack_blas1::operation op,long n,long prec,rt_pack &packed,long &index) {
    const size_t count=4*size_t(n)+1,nl=(prec+63)/64;
    if(!dynamic_storage_fits(count,n,prec)) return 1;
#ifndef MPLAPACK_MPFR_CUDA_HOST_EMULATION
    int devices=0;if(cudaGetDeviceCount(&devices) || !devices) return 1;
    mpfr_stack_guard stack(inverse_stack_bytes(prec,false));if(stack.error) return stack.error;
    factor_buffer<char> arena;factor_buffer<size_t> tops;
    const size_t pool=factor_pool,slab=16384*((prec+1023)/1024);
    if(arena.allocate(pool*slab) || tops.allocate(pool) || cudaMemset(tops.p,0,pool*sizeof(size_t))) return 1;
    complex_arena_guard guard(arena.p,slab,tops.p,prec);if(guard.error) return guard.error;
#else
    const size_t pool=1;
    long previous=mplapack_mpfr_rt_prec_h;mplapack_mpfr_rt_prec_h=prec;
    struct precision_guard {long p;~precision_guard(){mplapack_mpfr_rt_prec_h=p;}} guard{previous};
#endif
    factor_buffer<mplapack_rt_limb> limbs,values,one;factor_buffer<dynamic_real> matrix;
    factor_buffer<dynamic_cell> descriptions;factor_buffer<long> indices;
    if(limbs.allocate(count*nl) || values.allocate(pool*scalar_slots*nl) || one.allocate(nl) ||
       matrix.allocate(count) || descriptions.allocate(count) || indices.allocate(n)) return 1;
    std::vector<mplapack_rt_limb> literal(nl);literal[nl-1]=mplapack_rt_limb(1)<<63;
    if(one.upload(literal.data(),nl)) return 1;
#ifndef MPLAPACK_MPFR_CUDA_HOST_EMULATION
    dynamic_state_guard state(values.p,one.p);if(state.error) return state.error;
#else
    auto old_values=values_host,old_one=one_host;values_host=values.p;one_host=one.p;
    struct state_guard {mplapack_rt_limb *v,*o;~state_guard(){values_host=v;one_host=o;}} state{old_values,old_one};
#endif
    dynamic_launcher ctx;if(ctx.error) return ctx.error;
    std::vector<dynamic_cell> cells(count);
    for(size_t i=0;i<count;++i) cells[i]={limbs.p+i*nl,packed.exp[i],packed.kind[i],-2};
#ifndef MPLAPACK_MPFR_CUDA_HOST_EMULATION
    if(limbs.upload(reinterpret_cast<const mplapack_rt_limb*>(packed.limbs),count*nl) || descriptions.upload(cells.data(),count)) return 1;
    initialize_views<<<unsigned(std::min(size_t(factor_blocks),(count+factor_threads-1)/factor_threads)),factor_threads>>>(matrix.p,descriptions.p,count);
    if(cudaGetLastError() || cudaStreamSynchronize(nullptr)) return 1;
#else
    std::memcpy(limbs.p,packed.limbs,count*nl*sizeof(mplapack_rt_limb));
    for(size_t i=0;i<count;++i) {matrix.p[i].m=cells[i].m;matrix.p[i].exp=cells[i].exp;matrix.p[i].kind=cells[i].kind;matrix.p[i].slot=-2;}
#endif
    bool element=op==mplapack_blas1::axpy || op==mplapack_blas1::copy;
    dynamic_blas1_job job{dynamic_engine{},op,n,1,element?3:0,matrix.p,indices.p};
    ctx.execute(job,n);
    if(!element) {
        for(long stride=1;stride<n;) {
            job.phase=1;job.stride=stride;ctx.execute(job,(n-1)/stride/2+1);
            if(stride>n/2) break;stride*=2;
        }
        job.phase=2;ctx.execute(job,1);
    }
    ctx.finish();if(ctx.error) return ctx.error;
    long found=0;if(!element && indices.download(&found,1)) return 1;
#ifndef MPLAPACK_MPFR_CUDA_HOST_EMULATION
    if(cudaMemcpy(cells.data(),matrix.p,count*sizeof(dynamic_cell),cudaMemcpyDeviceToHost) ||
       limbs.download(reinterpret_cast<mplapack_rt_limb*>(packed.limbs),count*nl)) return 1;
#else
    for(size_t i=0;i<count;++i) {cells[i].kind=matrix.p[i].kind;cells[i].exp=matrix.p[i].exp;}
    std::memcpy(packed.limbs,limbs.p,count*nl*sizeof(mplapack_rt_limb));
#endif
    for(size_t i=0;i<count;++i) {packed.kind[i]=cells[i].kind;packed.exp[i]=cells[i].exp;}
#ifndef MPLAPACK_MPFR_CUDA_HOST_EMULATION
    if(state.restore() || guard.restore() || stack.restore()) return 1;
#endif
    index=found;return 0;
}
