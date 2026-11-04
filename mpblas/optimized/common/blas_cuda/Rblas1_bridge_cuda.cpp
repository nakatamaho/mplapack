// New standalone entries only; existing public dispatch defaults are unchanged.
#include <mpblas.h>
#include "../Rblas_cuda_entry_common.h"
#include "../Rblas_dev_types.h"
#include "../Rblas1_policy.h"
#include <type_traits>
#include <vector>
#include <cstring>
#include <stdexcept>
#if defined(MPLAPACK_GEMM_CUDA_HOST_EMULATION)
#include "../Rblas1_cuda_common.cuh"
#else
#include "../Rblas1_reduction_common.h"
namespace mplapack_blas1 {
template<class Ops,bool Complex> int run_bf(operation,long,
    const typename std::conditional<Complex,typename mplapack_blas_dev::cgpu_types<Ops>::DC,typename mplapack_blas_dev::cgpu_types<Ops>::DR>::type *,
    const typename std::conditional<Complex,typename mplapack_blas_dev::cgpu_types<Ops>::DC,typename mplapack_blas_dev::cgpu_types<Ops>::DR>::type *,
    state<typename mplapack_blas_dev::cgpu_types<Ops>::DR> &);
}
#endif
#if defined(MPLAPACK_BUILD_WITH_DD)
using OPS=mplapack_gemm::dd_ops;
#define PREFIX "MPLAPACK_DD"
#elif defined(MPLAPACK_BUILD_WITH_TD)
using OPS=mplapack_gemm::td_ops;
#define PREFIX "MPLAPACK_TD"
#elif defined(MPLAPACK_BUILD_WITH_QD)
using OPS=mplapack_gemm::qd_ops;
#define PREFIX "MPLAPACK_QD"
#endif
#if defined(PREFIX)
namespace {
using DR=mplapack_blas_dev::cgpu_types<OPS>::DR;
using DC=mplapack_blas_dev::cgpu_types<OPS>::DC;
bool finite(const DR &x) { for(double d:x.v.x) if(!std::isfinite(d)) return false;return true; }
bool finite(const DC &x) { return finite(x.real()) && finite(x.imag()); }
template<bool Complex,class H> bool blas1_entry(mplapack_blas1::operation op,long n,H *x,long incx,H *y,long incy,REAL &r,COMPLEX &c,long &index) {
    using namespace mplapack_blas1;
    using X=typename std::conditional<Complex,DC,DR>::type;
    const bool product=op==dot || op==dotc || op==dotu;
    if(op!=iamax && !fast_enabled(PREFIX "_MPBLAS_FAST")) return false;
    if(!stride_ok(n,incx) || !x || ((op==asum || op==iamax) && incx<0) ||
       (product && (!y || !stride_ok(n,incy))) || (Complex && op==mplapack_blas1::norm && n>std::numeric_limits<long>::max()/2)) return false;
    if(!mplapack_gemm::env_long(PREFIX "_BLAS_CUDA",1) || !mplapack_gemm::env_long(PREFIX "_BLAS_CUDA_L12",0) ||
       n<mplapack_gemm::env_long(PREFIX "_BLAS_CUDA_MIN",32768) || mplapack_gemm::gpu_device_count<OPS>()<=0) return false;
    try {
        std::vector<X> a(n),b(product?n:0);
        long ix=incx<0?-(n-1)*incx:0,iy=product && incy<0?-(n-1)*incy:0;
        for(long i=0;i<n;++i) {
            std::memcpy(&a[i],&x[ix+i*incx],sizeof(X));if(!finite(a[i])) return false;
            if(product) { std::memcpy(&b[i],&y[iy+i*incy],sizeof(X));if(!finite(b[i])) return false; }
        }
        state<DR> s;
        int error=run_bf<OPS,Complex>(op,n,a.data(),b.data(),s);
        if(error) {
            if(mplapack_gemm::env_long(PREFIX "_CUDA_VERBOSE",0)) std::fprintf(stderr,PREFIX " BLAS1 CUDA error %d\n",error);
            return false;
        }
        if(!finite(s.re) || !finite(s.im)) return false;
        std::memcpy(&r,&s.re,sizeof r);
        REAL im;std::memcpy(&im,&s.im,sizeof im);c=COMPLEX(r,im);index=s.index;
        return true;
    } catch(const std::bad_alloc &) { return false; } catch(const std::length_error &) { return false; }
}
}
#include "../Rblas1_entries.h"
#undef PREFIX
#endif
