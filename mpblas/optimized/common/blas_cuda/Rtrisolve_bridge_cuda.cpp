// Standalone opt-in solves; no CPU fallback or caller writes on false.
#include <mpblas.h>
#include "../Rblas_cuda_entry_common.h"
#include "../Rblas_dev_types.h"
#include "../Rtrisolve_policy.h"
#include <type_traits>
#include <vector>
#include <cstring>
#include <stdexcept>
#if defined(MPLAPACK_GEMM_CUDA_HOST_EMULATION)
#include "../Rtrisolve_cuda_common.cuh"
#else
namespace mplapack_trisolve {
template<class Ops,bool Complex> int run_bf(const shape &,
    const typename std::conditional<Complex,typename mplapack_blas_dev::cgpu_types<Ops>::DC,typename mplapack_blas_dev::cgpu_types<Ops>::DR>::type *,
    typename std::conditional<Complex,typename mplapack_blas_dev::cgpu_types<Ops>::DC,typename mplapack_blas_dev::cgpu_types<Ops>::DR>::type *);
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
namespace {
// Flatten only this policy shim so the shared inline env helper is not
// incidentally added to the library's public ABI by this translation unit.
__attribute__((flatten)) long tri_env(const char *name,long value) { return mplapack_gemm::env_long(name,value); }
using DR=mplapack_blas_dev::cgpu_types<OPS>::DR;
using DC=mplapack_blas_dev::cgpu_types<OPS>::DC;
bool finite(const DR &x) { for(double d:x.v.x) if(!std::isfinite(d)) return false;return true; }
bool finite(const DC &x) { return finite(x.real()) && finite(x.imag()); }
template<bool Complex,class H> bool tri_entry(int storage,const char *u,const char *t,const char *d,
    long n,long k,H *a,long lda,H *x,long inc) {
    using namespace mplapack_trisolve;
    using V=typename std::conditional<Complex,DC,DR>::type;
    shape s;
    if(!a || !x || !validate(s,storage,u,t,d,n,k,lda,inc,sizeof(H),sizeof(V))) return false;
    if(!tri_env(PREFIX "_BLAS_CUDA",1) || !tri_env(PREFIX "_BLAS_CUDA_L12",0) ||
       double(n)*(storage==1?k+1:n)<tri_env(PREFIX "_BLAS_CUDA_MIN",32768) ||
       mplapack_gemm::gpu_device_count<OPS>()<=0) return false;
    try {
        std::vector<V> matrix(s.count()),v(n);
        for(long j=0;j<n;++j) {
            for(long i=s.first(j);i<s.last(j);++i) {
                long at=s.index(i,j);std::memcpy(&matrix[at],&a[at],sizeof(V));if(!finite(matrix[at])) return false;
            }
            if(!s.unit) { long at=s.index(j,j);std::memcpy(&matrix[at],&a[at],sizeof(V));if(!finite(matrix[at])) return false; }
        }
        long start=inc<0?-(n-1)*inc:0;
        for(long i=0;i<n;++i) { std::memcpy(&v[i],&x[start+i*inc],sizeof(V));if(!finite(v[i])) return false; }
        int error=run_bf<OPS,Complex>(s,matrix.data(),v.data());
        if(error) {
            if(tri_env(PREFIX "_CUDA_VERBOSE",0)) std::fprintf(stderr,PREFIX " triangular solve CUDA error %d\n",error);
            return false;
        }
        for(const auto &z:v) if(!finite(z)) return false;
        for(long i=0;i<n;++i) std::memcpy(&x[start+i*inc],&v[i],sizeof(V));
        return true;
    } catch(const std::bad_alloc &) { return false; } catch(const std::length_error &) { return false; }
}
}
#include "../Rtrisolve_entries.h"
#undef PREFIX
