// Included by Rinverse_bridge_cuda.cpp: private checked packing, mutex and policy.
#include "../../common/Rtrisolve_policy.h"
#include "Rtrisolve_kernel_cuda.h"
namespace {
template<class Tag,bool Complex,class H> bool tri_packed(const mplapack_trisolve::shape &s,long p,bool fixed,H *a,H *x,long inc) {
    using V=value<Tag,Complex>;
    size_t total=size_t(s.count())+2*size_t(s.n)+1;
    if(total>std::numeric_limits<size_t>::max()/sizeof(V) || !memory_fits(total,sizeof(V),s.n,p,fixed)) return false;
    std::vector<V> matrix(s.count()),v(s.n);
    for(long j=0;j<s.n;++j) {
        for(long i=s.first(j);i<s.last(j);++i) { long at=s.index(i,j);if(!get(a[at],matrix[at],p)) return false; }
        if(!s.unit) { long at=s.index(j,j);if(!get(a[at],matrix[at],p)) return false; }
    }
    long start=inc<0?-(s.n-1)*inc:0;
    for(long i=0;i<s.n;++i) if(!get(x[start+i*inc],v[i],p)) return false;
    int error=run_tri<Tag,Complex>(s,p,matrix.data(),v.data());
    if(error) {
        if(env_long("MPLAPACK_MPFR_CUDA_VERBOSE",0)) std::fprintf(stderr,"MPFR triangular solve p=%ld CUDA error %d\n",p,error);
        return false;
    }
    for(const auto &z:v) {
        if constexpr(Complex) { if(!fits(z.v.re_) || !fits(z.v.im_)) return false; }
        else if(!fits(z.v)) return false;
    }
    for(long i=0;i<s.n;++i) put(v[i],x[start+i*inc],p);
    return true;
}
bool tri_dynamic(const mplapack_trisolve::shape &s,long p,mpfr_class *a,mpfr_class *x,long inc) {
    if(!dynamic_tri_memory_fits(s,p)) return false;
    size_t acount=s.count(),count=acount+2*size_t(s.n)+1;
    int nl=(p+63)/64;rt_storage stage;stage.resize(count,nl);
    std::fill(stage.kind.begin(),stage.kind.end(),RT_ZERO_KIND);
    for(long j=0;j<s.n;++j) {
        for(long i=s.first(j);i<s.last(j);++i) { long at=s.index(i,j);
            if(!rt_pack_value(a[at].get_mpfr_t(),p,nl,stage,at,RT_DEVICE_EMAX)) return false; }
        if(!s.unit) { long at=s.index(j,j);if(!rt_pack_value(a[at].get_mpfr_t(),p,nl,stage,at,RT_DEVICE_EMAX)) return false; }
    }
    long start=inc<0?-(s.n-1)*inc:0;
    for(long i=0;i<s.n;++i) if(!rt_pack_value(x[start+i*inc].get_mpfr_t(),p,nl,stage,acount+i,RT_DEVICE_EMAX)) return false;
    if(run_dynamic_tri(s,p,stage.pack)) return false;
    for(long i=0;i<s.n;++i) if(!rt_fits(stage,acount+i)) return false;
    for(long i=0;i<s.n;++i) rt_unpack(stage,acount+i,p,nl,x[start+i*inc].get_mpfr_t());
    return true;
}
template<bool Complex,class H> bool tri_entry(int storage,const char *u,const char *t,const char *d,
    long n,long k,H *a,long lda,H *x,long inc) {
    mplapack_trisolve::shape s;
    long p=mpfrxx::default_precision_bits();
    if(p<MPFR_PREC_MIN || p>std::numeric_limits<int>::max()-1024L || (Complex && p>4096) || mpfr_get_default_rounding_mode()!=MPFR_RNDN ||
       mpfr_get_emax()!=RT_DEVICE_EMAX || mpfr_get_emin()!=-RT_DEVICE_EMAX) return false;
    if(Complex && (mpfrxx::default_mpc_real_precision_bits()!=p || mpfrxx::default_mpc_imag_precision_bits()!=p ||
       mpfrxx::default_mpc_rounding_mode()!=MPC_RNDNN)) return false;
    if(!a || !x || !mplapack_trisolve::validate(s,storage,u,t,d,n,k,lda,inc,sizeof(H),2*(size_t((p+63)/64)*8+64))) return false;
    if(!env_long("MPLAPACK_MPFR_CUDA",1) || !env_long("MPLAPACK_MPFR_CUDA_L12",0) ||
       double(n)*(storage==1?k+1:n)<env_long("MPLAPACK_MPFR_CUDA_MIN_MNK",32768)) return false;
    try {
        inverse_lock lock(gpu_mutex());
        if(p>4096) {
            if(!env_long("MPLAPACK_MPFR_CUDA_RUNTIME",1)) return false;
            if constexpr(!Complex) return tri_dynamic(s,p,a,x,inc);
            return false;
        }
        if(!env_long("MPLAPACK_MPFR_CUDA_FORCE_RUNTIME",0)) switch(p) {
#define TRI_FIXED(PB) case PB: return tri_packed<factor_fixed_tag<PB>,Complex>(s,p,true,a,x,inc);
            MPLAPACK_MPFR_CUDA_FIXED_PRECISIONS(TRI_FIXED)
#undef TRI_FIXED
        }
        if(!env_long("MPLAPACK_MPFR_CUDA_RUNTIME",1)) return false;
#define TRI_RT(NL) if(p<=64*NL) return tri_packed<factor_runtime_tag<NL>,Complex>(s,p,false,a,x,inc);
        MPLAPACK_MPFR_CUDA_RT_TIERS(TRI_RT)
#undef TRI_RT
        return false;
    } catch(const std::bad_alloc &) { return false; } catch(const std::length_error &) { return false; }
}
}
#include "../../common/Rtrisolve_entries.h"
