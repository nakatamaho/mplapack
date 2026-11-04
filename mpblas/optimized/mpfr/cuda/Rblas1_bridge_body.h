// Included by Rinverse_bridge_cuda.cpp to share its private packing adapters.
#include <mpblas.h>
#include "../../common/Rblas_cuda_entry_common.h"
#include "../../common/Rblas1_policy.h"
#include "Rblas1_kernel_cuda.h"
namespace {
template<class Tag,bool Complex,class H> bool blas1_packed(mplapack_blas1::operation op,long n,long p,bool fixed,
    H *x,long incx,H *y,long incy,const mpfr_class &alpha,mpfr_class &r,mpc_class &c,long &index) {
    using namespace mplapack_blas1;
    using V=value<Tag,Complex>;using R=value<Tag,false>;using S=state<R>;
    const bool product=op==dot || op==dotc || op==dotu,element=op==axpy || op==copy;
    // Account for both operands and the reduction state (including complex norm's
    // two real leaves), plus the existing arena, pool and inherited stack budget.
    size_t count=Complex && op==mplapack_blas1::norm?2*size_t(n):size_t(n);
    if(!memory_fits(count,sizeof(S)+2*sizeof(V),n,p,fixed)) return false;
    std::vector<V> a(n),b(product || element?n:0);
    long ix=incx<0?-(n-1)*incx:0,iy=(product || element) && incy<0?-(n-1)*incy:0;
    R al;if(!get(alpha,al,p)) return false;
    for(long i=0;i<n;++i) {
        if(!get(x[ix+i*incx],a[i],p)) return false;
        if(product || op==axpy) { if(!get(y[iy+i*incy],b[i],p)) return false; }
        else if(op==copy && !precision(y[iy+i*incy],p)) return false;
    }
    S s;
    int error=run_blas1<Tag,Complex>(op,n,p,a.data(),b.data(),al,s);
    if(error) {
        if(env_long("MPLAPACK_MPFR_CUDA_VERBOSE",0)) std::fprintf(stderr,"MPFR BLAS1 op=%d p=%ld CUDA error %d\n",int(op),p,error);
        return false;
    }
    if(element) {
        // The checked launch has validated all intermediate arithmetic; all
        // downloads/restoration completed before any caller output is changed.
        for(long i=0;i<n;++i) put(b[i],y[iy+i*incy],p);
    } else {
        if(!fits(s.re.v) || !fits(s.im.v)) return false;
        put(s.re,r,p);value<Tag,true> z(s.re,s.im);put(z,c,p);index=s.index;
    }
    return true;
}
bool blas1_dynamic(mplapack_blas1::operation op,long n,long p,mpfr_class *x,long incx,mpfr_class *y,long incy,
    const mpfr_class &alpha,mpfr_class &r,mpc_class &c,long &index) {
    using namespace mplapack_blas1;
    if(!dynamic_blas1_memory_fits(n,p)) return false;
    const bool product=op==dot,element=op==axpy || op==copy;
    const size_t count=4*size_t(n)+1;int nl=(p+63)/64;
    rt_storage stage;stage.resize(count,nl);
    std::fill(stage.kind.begin(),stage.kind.end(),RT_ZERO_KIND);
    long ix=incx<0?-(n-1)*incx:0,iy=(product || element) && incy<0?-(n-1)*incy:0;
    if(!rt_pack_value(alpha.get_mpfr_t(),p,nl,stage,4*size_t(n),RT_DEVICE_EMAX)) return false;
    for(long i=0;i<n;++i) {
        if(!rt_pack_value(x[ix+i*incx].get_mpfr_t(),p,nl,stage,i,RT_DEVICE_EMAX)) return false;
        if(product || op==axpy) { if(!rt_pack_value(y[iy+i*incy].get_mpfr_t(),p,nl,stage,n+i,RT_DEVICE_EMAX)) return false; }
        else if(op==copy && !precision(y[iy+i*incy],p)) return false;
    }
    long found=0;int error=run_dynamic_blas1(op,n,p,stage.pack,found);
    if(error) {
        if(env_long("MPLAPACK_MPFR_CUDA_VERBOSE",0)) std::fprintf(stderr,"MPFR dynamic BLAS1 op=%d p=%ld CUDA error %d\n",int(op),p,error);
        return false;
    }
    if(element) {
        for(long i=0;i<n;++i) if(!rt_fits(stage,n+i)) return false;
        for(long i=0;i<n;++i) rt_unpack(stage,n+i,p,nl,y[iy+i*incy].get_mpfr_t());
    } else {
        if(!rt_fits(stage,2*size_t(n))) return false;
        rt_unpack(stage,2*size_t(n),p,nl,r.get_mpfr_t());c=mpc_class(r);index=found;
    }
    return true;
}
template<bool Complex,class H> bool blas1_entry(mplapack_blas1::operation op,long n,H *x,long incx,H *y,long incy,
    mpfr_class &r,mpc_class &c,long &index,const mpfr_class &alpha=mpfr_class(0.0)) {
    using namespace mplapack_blas1;
    const bool product=op==dot || op==dotc || op==dotu,element=op==axpy || op==copy;
    if(op!=iamax && !element && !fast_enabled("MPLAPACK_MPFR_MPBLAS_FAST")) return false;
    if(!stride_ok(n,incx) || !x || ((op==asum || op==iamax) && incx<0) ||
       ((product || element) && (!y || !stride_ok(n,incy))) || (Complex && op==mplapack_blas1::norm && n>std::numeric_limits<long>::max()/2) ||
       (op==axpy && alpha==mpfr_class(0.0))) return false;
    if(!env_long("MPLAPACK_MPFR_CUDA",1) || !env_long("MPLAPACK_MPFR_CUDA_L12",0) ||
       n<env_long("MPLAPACK_MPFR_CUDA_MIN_MNK",32768)) return false;
    long p=mpfrxx::default_precision_bits();
    if(p<MPFR_PREC_MIN || p>std::numeric_limits<int>::max()-1024L || (Complex && p>4096) || mpfr_get_default_rounding_mode()!=MPFR_RNDN ||
       mpfr_get_emax()!=RT_DEVICE_EMAX || mpfr_get_emin()!=-RT_DEVICE_EMAX) return false;
    if(Complex && (mpfrxx::default_mpc_real_precision_bits()!=p || mpfrxx::default_mpc_imag_precision_bits()!=p ||
       mpfrxx::default_mpc_rounding_mode()!=MPC_RNDNN)) return false;
    if(size_t(n)>(std::numeric_limits<size_t>::max()/4-1)/size_t((p+63)/64*8+64)) return false;
    try {
        inverse_lock lock(gpu_mutex());
        if(p>4096) {
            if(!env_long("MPLAPACK_MPFR_CUDA_RUNTIME",1)) return false;
            if constexpr(!Complex) return blas1_dynamic(op,n,p,x,incx,y,incy,alpha,r,c,index);
            return false;
        }
        if(!env_long("MPLAPACK_MPFR_CUDA_FORCE_RUNTIME",0)) switch(p) {
#define BLAS1_FIXED(PB) case PB: return blas1_packed<factor_fixed_tag<PB>,Complex>(op,n,p,true,x,incx,y,incy,alpha,r,c,index);
            MPLAPACK_MPFR_CUDA_FIXED_PRECISIONS(BLAS1_FIXED)
#undef BLAS1_FIXED
        }
        if(!env_long("MPLAPACK_MPFR_CUDA_RUNTIME",1)) return false;
#define BLAS1_RT(NL) if(p<=64*NL) return blas1_packed<factor_runtime_tag<NL>,Complex>(op,n,p,false,x,incx,y,incy,alpha,r,c,index);
        MPLAPACK_MPFR_CUDA_RT_TIERS(BLAS1_RT)
#undef BLAS1_RT
        return false;
    } catch(const std::bad_alloc &) { return false; } catch(const std::length_error &) { return false; }
}
}
#include "../../common/Rblas1_entries.h"
bool Raxpy_mpfr_cuda(mplapackint const n,mpfr_class const da,mpfr_class *dx,mplapackint const incx,mpfr_class *dy,mplapackint const incy) {
    mpfr_class r;mpc_class c;long index=0;
    return blas1_entry<false>(mplapack_blas1::axpy,n,dx,incx,dy,incy,r,c,index,da);
}
bool Rcopy_mpfr_cuda(mplapackint const n,mpfr_class *dx,mplapackint const incx,mpfr_class *dy,mplapackint const incy) {
    mpfr_class r;mpc_class c;long index=0;
    return blas1_entry<false>(mplapack_blas1::copy,n,dx,incx,dy,incy,r,c,index);
}
