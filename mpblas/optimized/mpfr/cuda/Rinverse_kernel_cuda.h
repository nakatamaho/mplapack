// Internal MPFR inverse interface. No host MPFR classes in device code.
#ifndef MPLAPACK_MPFR_INVERSE_KERNEL_CUDA_H
#define MPLAPACK_MPFR_INVERSE_KERNEL_CUDA_H
#include "Rinverse_types.h"
#include "Rfactor_kernel_cuda.h"
namespace mplapack_mpfr_inverse {
struct shape {
    int kind=0; bool upper=false,unit=false;
    long n=0,nbtri=1,nblau=1,nb=1,nbmin=2,columns=1;
};
template<class Tag,bool Complex> __attribute__((visibility("hidden")))
int run(const shape &,long,value<Tag,Complex> *,value<Tag,Complex> *,const long *,long &,mplapack_mpfr_cuda::factor_times &);
__attribute__((visibility("hidden"))) bool memory_fits(size_t count,size_t bytes,long n,long prec,bool fixed);
}
#endif
