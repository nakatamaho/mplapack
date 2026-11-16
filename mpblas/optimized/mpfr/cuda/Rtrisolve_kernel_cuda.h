#ifndef MPLAPACK_MPFR_TRISOLVE_KERNEL_CUDA_H
#define MPLAPACK_MPFR_TRISOLVE_KERNEL_CUDA_H
#include "Rinverse_kernel_cuda.h"
#include "../../common/Rtrisolve_jobs.h"
namespace mplapack_mpfr_inverse __attribute__((visibility("hidden"))) {
bool dynamic_tri_memory_fits(const mplapack_trisolve::shape &,long);
int run_dynamic_tri(const mplapack_trisolve::shape &,long,mplapack_mpfr_cuda::rt_pack &);
template<class Tag,bool Complex> int run_tri(const mplapack_trisolve::shape &,long,value<Tag,Complex> *,value<Tag,Complex> *);
}
#endif
