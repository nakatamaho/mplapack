// Unbounded-precision real inverse values. Matrix values are views of packed
// limbs; scalar temporaries own slots in a per-physical-thread limb pool.
// MPC scratch may therefore be reclaimed between operations without ever
// invalidating a live reference-template temporary. No device heap is used.
#include "Rinverse_dynamic_cuda.h"
#include "Rinverse_device_body.cuh"
#include <new>
#include <cstring>
namespace mplapack_mpfr_inverse {
namespace {
constexpr int scalar_slots=64;
struct dynamic_cell { mplapack_rt_limb *m; long exp; int kind,slot; };
static mplapack_rt_limb *values_host=nullptr,*one_host=nullptr;
static thread_local unsigned long long slots_host=0;
#ifndef MPLAPACK_MPFR_CUDA_HOST_EMULATION
static __device__ mplapack_rt_limb *values_device,*one_device;
static __device__ unsigned long long slots_device[factor_pool];
#endif
MPLAPACK_DEV_HD unsigned long long &slots() {
#if defined(__CUDA_ARCH__)
    return slots_device[blockIdx.x*blockDim.x+threadIdx.x];
#else
    return slots_host;
#endif
}
MPLAPACK_DEV_HD mplapack_rt_limb *value_limbs(int slot) {
#if defined(__CUDA_ARCH__)
    return values_device+((size_t(blockIdx.x)*blockDim.x+threadIdx.x)*scalar_slots+slot)*((rt_prec()+63)/64);
#else
    return values_host+size_t(slot)*((rt_prec()+63)/64);
#endif
}
MPLAPACK_DEV_HD mplapack_rt_limb *literal_one() {
#if defined(__CUDA_ARCH__)
    return one_device;
#else
    return one_host;
#endif
}
struct view_tag {};
struct dynamic_real {
    mplapack_rt_limb *m=nullptr;
    long exp=0;
    int kind=2,slot=-1; // -1 literal (0, +/-1), -2 matrix view, >=0 scalar slot
    MPLAPACK_DEV_HD dynamic_real() {}
    MPLAPACK_DEV_HD dynamic_real(double d):exp(d==0?0:1),kind(d==0?2:d<0?-3:3) {}
    MPLAPACK_DEV_HD dynamic_real(view_tag,const dynamic_cell &c):m(c.m),exp(c.exp),kind(c.kind),slot(-2) {}
    MPLAPACK_DEV_HD dynamic_real(const dynamic_real &r) { *this=r; }
    MPLAPACK_DEV_HD dynamic_real(dynamic_real &&r):m(r.m),exp(r.exp),kind(r.kind),slot(r.slot) {
        r.m=nullptr;r.kind=2;r.exp=0;r.slot=-1;
    }
    MPLAPACK_DEV_HD ~dynamic_real() { if(slot>=0) slots()&=~(1ULL<<slot); }
    MPLAPACK_DEV_HD bool storage() {
        if(m) return true;
        for(int i=0;i<scalar_slots;++i) if(!(slots()&(1ULL<<i))) {
            slots()|=1ULL<<i;slot=i;m=value_limbs(i);return true;
        }
        failure()=true;return false;
    }
    MPLAPACK_DEV_HD dynamic_real &operator=(const dynamic_real &r) {
        if(this==&r) return *this;
        if(!m && !r.m) { kind=r.kind;exp=r.exp;return *this; }
        if(!storage()) return *this;
        kind=r.kind;exp=r.exp;
        if(kind==3 || kind==-3) {
            const auto *src=r.m?r.m:literal_one();
            for(long i=0;i<(rt_prec()+63)/64;++i) m[i]=src[i];
        }
        return *this;
    }
    MPLAPACK_DEV_HD bool valid() const {
        int k=kind<0?-kind:kind;
        return k==2 || (k==3 && exp>=-1073741823L && exp<=1073741823L);
    }
    MPLAPACK_DEV_HD bool safe_operand() const {
        const long bound=(1073741823L-4L*rt_prec()-128)/4;
        return valid() && (kind==2 || kind==-2 || (exp>=-bound && exp<=bound));
    }
    MPLAPACK_DEV_HD dynamic_real operator-() const { dynamic_real r(*this);r.kind=-r.kind;return r; }
    MPLAPACK_DEV_HD dynamic_real operator+() const { return dynamic_real(*this); }
#ifdef MPLAPACK_RT_ACTIVE
    MPLAPACK_DEV_HD void view(MPLAPACK_RT(mpfr_ptr) x) const {
        MPLAPACK_RT(mpfr_custom_init_set)(x,kind,exp,rt_prec(),m?m:literal_one());
    }
#endif
    MPLAPACK_MPFR_DEV_OP static dynamic_real op(const dynamic_real &a,const dynamic_real &b,int which) {
        dynamic_real r;
        if(!a.safe_operand() || !b.safe_operand() || (which==3 && (b.kind==2 || b.kind==-2))) {
            failure()=true;return r;
        }
#ifdef MPLAPACK_RT_ACTIVE
        if(!r.storage()) return r;
        MPLAPACK_RT_DECL(x);MPLAPACK_RT_DECL(y);MPLAPACK_RT_DECL(z);
        a.view(x);b.view(y);r.view(z);factor_arena_reset();
        if(which==0) MPLAPACK_RT(mpfr_add)(z,x,y,MPLAPACK_RTU(MPFR_RNDN));
        else if(which==1) MPLAPACK_RT(mpfr_sub)(z,x,y,MPLAPACK_RTU(MPFR_RNDN));
        else if(which==2) MPLAPACK_RT(mpfr_mul)(z,x,y,MPLAPACK_RTU(MPFR_RNDN));
        else MPLAPACK_RT(mpfr_div)(z,x,y,MPLAPACK_RTU(MPFR_RNDN));
        MPLAPACK_RT_TAKE(r,z);
        if(!r.valid()) failure()=true;
#else
        failure()=true;
#endif
        return r;
    }
#define DYNAMIC_OP(OP,ASSIGN,N) \
    friend MPLAPACK_DEV_HD dynamic_real operator OP(const dynamic_real &a,const dynamic_real &b) { return op(a,b,N); } \
    MPLAPACK_DEV_HD dynamic_real &operator ASSIGN(const dynamic_real &b) { return *this=op(*this,b,N); }
    DYNAMIC_OP(+,+=,0)
    DYNAMIC_OP(-,-=,1)
    DYNAMIC_OP(*,*=,2)
    DYNAMIC_OP(/,/=,3)
#undef DYNAMIC_OP
    friend MPLAPACK_DEV_HD bool operator==(const dynamic_real &a,const dynamic_real &b) {
        if(!a.valid() || !b.valid()) { failure()=true;return false; }
#ifdef MPLAPACK_RT_ACTIVE
        MPLAPACK_RT_DECL(x);MPLAPACK_RT_DECL(y);a.view(x);b.view(y);
        return MPLAPACK_RT(mpfr_cmp)(x,y)==0;
#else
        return a.kind==b.kind && a.exp==b.exp;
#endif
    }
    friend MPLAPACK_DEV_HD bool operator!=(const dynamic_real &a,const dynamic_real &b) { return !(a==b); }
};
static_assert(sizeof(dynamic_real)==sizeof(dynamic_cell),"runtime view layout");
static_assert(offsetof(dynamic_real,m)==offsetof(dynamic_cell,m) && offsetof(dynamic_real,exp)==offsetof(dynamic_cell,exp) &&
              offsetof(dynamic_real,kind)==offsetof(dynamic_cell,kind) && offsetof(dynamic_real,slot)==offsetof(dynamic_cell,slot),"runtime view fields");
struct dynamic_engine : inverse_engine {
    L3_HD void thread(long tid) { inverse_engine::thread(tid);slots()=0; }
};
using dynamic_arithmetic=mplapack_inverse_gpu::inverse_arithmetic<dynamic_real,dynamic_real,false>;
using dynamic_inverse_job=mplapack_inverse_gpu::inverse_job<dynamic_arithmetic>;
// A kernel argument must be trivially copyable. Transfer only a POD job
// description; construct the owning scalar job after entering the thread.
struct dynamic_job {
    dynamic_engine e;
    dynamic_real *a,*b,*c;long *status,*piv;
    long m,n,k,ld,col;int op,alpha;bool upper,left,unit;char ta,tb;
    dynamic_job(const dynamic_inverse_job &j):a(j.a),b(j.b),c(j.c),status(j.status),piv(j.piv),
        m(j.m),n(j.n),k(j.k),ld(j.ld),col(j.col),op(j.op),alpha(j.alpha.kind==2 || j.alpha.kind==-2?0:j.alpha.kind<0?-1:1),
        upper(j.upper),left(j.left),unit(j.unit),ta(j.ta),tb(j.tb) {}
    L3_HD void operator()(long lane) {
        dynamic_inverse_job j;
        j.a=a;j.b=b;j.c=c;j.status=status;j.piv=piv;j.m=m;j.n=n;j.k=k;j.ld=ld;j.col=col;
        j.op=op;j.alpha=dynamic_real(alpha);j.upper=upper;j.left=left;j.unit=unit;j.ta=ta;j.tb=tb;
        j(lane);
        // No scalar temporary may escape a lane. Views and literal alpha own
        // no slots; an exhausted/leaked pool rejects the entire publication.
        if(slots()!=0) failure()=true;
        if(failure()) e.ok=false;
    }
};
static_assert(std::is_trivially_copyable<dynamic_job>::value,"CUDA job argument must be POD");
struct dynamic_launcher : factor_context<dynamic_engine> {
    dynamic_launcher():factor_context<dynamic_engine>(dynamic_engine{},dynamic_engine::M{},1,1) {}
    long *status_data() { return status.p; }
    long feedback(long *) { return factor_context<dynamic_engine>::feedback(); }
    void launch(const dynamic_inverse_job &j,long lanes) { this->execute(dynamic_job(j),lanes); }
    void finish() {
        factor_context<dynamic_engine>::finish();int f=0;
        if(!error) error=failed.download(&f,1);
        if(!error && f) error=1;
    }
};
#ifndef MPLAPACK_MPFR_CUDA_HOST_EMULATION
__global__ void initialize_views(dynamic_real *out,const dynamic_cell *cells,size_t count) {
    for(size_t i=size_t(blockIdx.x)*blockDim.x+threadIdx.x;i<count;i+=size_t(gridDim.x)*blockDim.x)
        new(out+i) dynamic_real(view_tag{},cells[i]);
}
struct dynamic_state_guard {
    mplapack_rt_limb *old_values=nullptr,*old_one=nullptr;
    bool known=false,changed=false;int error=0;
    dynamic_state_guard(mplapack_rt_limb *v,mplapack_rt_limb *one) {
        error=cudaMemcpyFromSymbol(&old_values,values_device,sizeof(v));
        if(!error) error=cudaMemcpyFromSymbol(&old_one,one_device,sizeof(one));
        if(error) return;
        known=changed=true;
        error=cudaMemcpyToSymbol(values_device,&v,sizeof(v));
        if(!error) error=cudaMemcpyToSymbol(one_device,&one,sizeof(one));
    }
    int restore() {
        if(!known || !changed) return 0;
        int a=cudaMemcpyToSymbol(values_device,&old_values,sizeof(old_values));
        int b=cudaMemcpyToSymbol(one_device,&old_one,sizeof(old_one));
        if(!a && !b) changed=false;
        return a?a:b;
    }
    ~dynamic_state_guard() { (void)restore(); }
};
#endif
}
bool dynamic_real_memory_fits(const shape &s,long p) {
#ifndef MPLAPACK_MPFR_CUDA_HOST_EMULATION
    size_t free=0,total=0,old=0;int device=0;cudaDeviceProp prop{};
    if(cudaGetDevice(&device) || cudaGetDeviceProperties(&prop,device) || cudaMemGetInfo(&free,&total) ||
       cudaDeviceGetLimit(&old,cudaLimitStackSize)) return false;
    size_t nl=(p+63)/64,count=size_t(s.n)*(s.n+(s.kind==3?s.columns:0));
    size_t overhead=size_t(prop.multiProcessorCount)*prop.maxThreadsPerMultiProcessor*std::max(old,inverse_stack_bytes(p,false));
    overhead+=size_t(factor_pool)*(16384*((p+1023)/1024)+sizeof(size_t)+scalar_slots*nl*sizeof(mplapack_rt_limb));
    overhead+=1048576+size_t(s.n)*sizeof(long)+nl*sizeof(mplapack_rt_limb);
    return free>=overhead && count<=(free-overhead)/(nl*sizeof(mplapack_rt_limb)+sizeof(dynamic_real)+sizeof(dynamic_cell));
#else
    (void)s;(void)p;return true;
#endif
}
int run_dynamic_real(const shape &s,long prec,rt_pack &packed,const long *piv,long &info,factor_times &times) {
    using Clock=std::chrono::steady_clock;
    const size_t acount=size_t(s.n)*s.n,wcount=s.kind==3?size_t(s.n)*s.columns:0,count=acount+wcount,nl=(prec+63)/64;
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
    factor_buffer<mplapack_rt_limb> limbs,values,one;
    factor_buffer<dynamic_real> matrix;
    factor_buffer<dynamic_cell> descriptions;
    factor_buffer<long> pivots;
    if(limbs.allocate(count*nl) || values.allocate(pool*scalar_slots*nl) || one.allocate(nl) ||
       matrix.allocate(count) || descriptions.allocate(count) || (wcount && pivots.allocate(s.n))) return 1;
    std::vector<mplapack_rt_limb> literal(nl,0);literal[nl-1]=mplapack_rt_limb(1)<<63;
    if(one.upload(literal.data(),nl)) return 1;
#ifndef MPLAPACK_MPFR_CUDA_HOST_EMULATION
    dynamic_state_guard state(values.p,one.p);if(state.error) return state.error;
#else
    auto old_values=values_host,old_one=one_host;values_host=values.p;one_host=one.p;
    struct state_guard {
        mplapack_rt_limb *v,*o;
        ~state_guard() {values_host=v;one_host=o;}
    } state{old_values,old_one};
#endif
    mplapack_inverse_gpu::inverse_schedule<dynamic_arithmetic,dynamic_launcher> ctx;if(ctx.error) return ctx.error;
    std::vector<dynamic_cell> cells(count);
    for(size_t i=0;i<count;++i) cells[i]={limbs.p+i*nl,packed.exp[i],packed.kind[i],-2};
    auto start=Clock::now();
    if(wcount && pivots.upload(piv,s.n)) return 1;
#ifndef MPLAPACK_MPFR_CUDA_HOST_EMULATION
    if(limbs.upload(reinterpret_cast<const mplapack_rt_limb*>(packed.limbs),count*nl) || descriptions.upload(cells.data(),count)) return 1;
    initialize_views<<<unsigned(std::min(size_t(factor_blocks),(count+factor_threads-1)/factor_threads)),factor_threads>>>(matrix.p,descriptions.p,count);
    if(cudaGetLastError() || cudaStreamSynchronize(nullptr)) return 1;
#else
    std::memcpy(limbs.p,packed.limbs,count*nl*sizeof(mplapack_rt_limb));
    for(size_t i=0;i<count;++i) {matrix.p[i].m=cells[i].m;matrix.p[i].exp=cells[i].exp;matrix.p[i].kind=cells[i].kind;matrix.p[i].slot=-2;}
#endif
    auto uploaded=Clock::now();long result=0;
    if(s.kind==0 || s.kind==2) result=ctx.tri(s.upper,s.unit,s.n,s.nbtri,matrix.p,s.n);
    if(!result && (s.kind==1 || s.kind==2)) ctx.lauum(s.upper,s.n,s.nblau,matrix.p,s.n);
    if(s.kind==3) result=ctx.getri(s.n,s.nbtri,s.nb,s.nbmin,matrix.p,matrix.p+acount,pivots.p,s.n);
    ctx.finish();if(ctx.error) return ctx.error;
    auto computed=Clock::now();
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
    info=result;
    times.upload=std::chrono::duration<double>(uploaded-start).count();
    times.factor=std::chrono::duration<double>(computed-uploaded).count();
    times.download=std::chrono::duration<double>(Clock::now()-computed).count();
    return 0;
}
}
