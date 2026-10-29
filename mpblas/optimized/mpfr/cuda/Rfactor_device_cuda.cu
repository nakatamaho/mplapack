// MPFR-only real device-resident factorizations, fixed-precision engine.
#include "Rfactor_device_body.cuh"
namespace mplapack_mpfr_cuda {
namespace {
template<int PB> struct factor_fixed_update : l3_fixed_engine<PB,true> {
    using M=typename l3_fixed_engine<PB,true>::M;
    L3_HD void thread(long) {}
    L3_HD static M shift(M a,long n) { a.p+=n; return a; }
};
template<int PB> struct factor_fixed_engine : l3_fixed_engine<PB,true> {
    using Base=l3_fixed_engine<PB,true>;
    using V=typename Base::V; using M=typename Base::M;
    static constexpr bool runtime=false;
    using Update=factor_fixed_update<PB>;
    L3_HD Update update() const { return Update{}; }
#if defined(__CUDACC__)
#define FACTOR_OP __host__ __device__ __noinline__
#else
#define FACTOR_OP inline
#endif
    FACTOR_OP void add(V &r,const V &a,const V &b) const { Base::add(r,a,b); }
    FACTOR_OP void sub(V &r,const V &a,const V &b) const { Base::sub(r,a,b); }
    FACTOR_OP void mul(V &r,const V &a,const V &b) const { Base::mul(r,a,b); }
#undef FACTOR_OP
    L3_HD void thread(long) {}
    L3_HD static M shift(M a,long n) { a.p+=n; return a; }
    L3_HD void absolute(V &r,const V &a) const { r=a; r.sign=1; }
    L3_HD int compare(const V &a,const V &b) const { return cu_fp::cu_cmp<PB>(a,b); }
#if defined(__CUDACC__)
    __noinline__
#endif
    L3_HD void root(V &r,const V &a) const {
        if(cu_fp::cu_sqrt_rn<PB>(r,a)) { this->ok=false; r.set_zero(); }
        this->checked(r);
    }
};
}
size_t factor_stack_bytes(long prec,bool runtime) {
    return runtime?runtime_stack_bytes(prec):fixed_stack_bytes(prec);
}
bool factor_memory_fits(long m,long n,long prec,bool fixed) {
#ifndef MPLAPACK_MPFR_CUDA_HOST_EMULATION
    size_t free=0,total=0;
    int device=0; cudaDeviceProp properties{};
    if(cudaGetDevice(&device)!=cudaSuccess || cudaGetDeviceProperties(&properties,device)!=cudaSuccess ||
       cudaMemGetInfo(&free,&total)!=cudaSuccess) return false;
    // Account for the stack reservation before allocating/packing the matrix.
    size_t overhead=size_t(properties.multiProcessorCount)*properties.maxThreadsPerMultiProcessor*factor_stack_bytes(prec,!fixed);
    overhead+=1048576+(size_t(std::min(m,n))+(size_t(m)+127)/128+2)*sizeof(long);
    if(!fixed) overhead+=size_t(2*factor_pool)*(16384*((prec+1023)/1024)+sizeof(size_t)+8*((prec+63)/64)*8);
    if(free<overhead) return false;
    size_t bytes=((prec+63)/64)*8+(fixed?16:9);
    if(size_t(n)>(free-overhead)/size_t(m)/bytes) return false;
    return true;
#else
    (void)m; (void)n; (void)prec; (void)fixed; return true;
#endif
}
template<int PB> int factor_fixed(const factor_shape &s,real_t<PB> *a,const real_t<PB> &safe,long *piv,long &info,factor_times &times) {
    using E=factor_fixed_engine<PB>; using V=real_t<PB>;
    using Clock=std::chrono::steady_clock;
    size_t count=(size_t)s.m*s.n;
#ifndef MPLAPACK_MPFR_CUDA_HOST_EMULATION
    int ndev=0;
    if(cudaGetDeviceCount(&ndev)!=cudaSuccess || !ndev) return 1;
    factor_stack_guard stack(factor_stack_bytes(PB,false));
    if(stack.error) return stack.error;
    size_t available=0,total=0;
    if(cudaMemGetInfo(&available,&total)!=cudaSuccess ||
       count>(available-std::min(available,size_t(1048576)))/sizeof(V)) return 1;
#endif
    factor_buffer<V> matrix,constants;
    if(matrix.allocate(count) || constants.allocate(4)) return 1;
    V c[4]={V::from_double(1),V::from_double(-1),safe,V::from_double(0)};
    if(constants.upload(c,4)) return 1;
    factor_context<E> ctx(E{},typename E::M{constants.p,4},s.m,s.n);
    if(ctx.error) return ctx.error;
    auto start=Clock::now();
    if(matrix.upload(a,count)) return 1;
#ifndef MPLAPACK_MPFR_CUDA_HOST_EMULATION
    if(cudaStreamSynchronize(nullptr)!=cudaSuccess) return 1;
#endif
    auto uploaded=Clock::now();
    long result=ctx.factor(s,typename E::M{matrix.p,s.m});
    if(ctx.error) return ctx.error;
    auto factored=Clock::now();
    if(matrix.download(a,count) || (!s.chol && ctx.pivots.download(piv,std::min(s.m,s.n)))) return 1;
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
#define INSTANTIATE(PB) template int factor_fixed<PB>(const factor_shape &,real_t<PB> *,const real_t<PB> &,long *,long &,factor_times &);
MPLAPACK_MPFR_CUDA_FIXED_PRECISIONS(INSTANTIATE)
#undef INSTANTIATE
}
