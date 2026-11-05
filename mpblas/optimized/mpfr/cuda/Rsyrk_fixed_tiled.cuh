// Fixed MPFR adapter to the same ordered SYRK/SYR2K tiles as BF.
#ifndef MPLAPACK_RSYRK_FIXED_TILED_CUH
#define MPLAPACK_RSYRK_FIXED_TILED_CUH
#include "blas3_kernel_cuda.h"
#include "../../common/Rblas_cuda_kernels_common.cuh"
namespace mplapack_mpfr_cuda {
template<int PB> struct syrk_fixed_ops {
    using F=real_t<PB>;
    static L3_HD F from_d(double x) { return F::from_double(x); }
    static L3_HD F add(const F &a,const F &b) { return a+b; }
    static L3_HD F sub(const F &a,const F &b) { return a-b; }
    static L3_HD F mul(const F &a,const F &b) { return a*b; }
    static L3_HD bool eq(const F &a,const F &b) { return cu_fp::cu_cmp<PB>(a,b)==0; }
    // Only needed to instantiate the common host dispatcher; the caller
    // selects this adapter exclusively for SYRK/SYR2K, never for TRSM.
    static L3_HD F div(const F &a,const F &b) { F r; cu_fp::cu_fdiv_rn<PB>(r,a,b); return r; }
};
static inline bool fixed_syrk_tiled_shape(const l3_shape &s) {
    return (s.op==mplapack_blas::L3_SYRK || s.op==mplapack_blas::L3_SYR2K) && s.k>0;
}
template<int PB> static int fixed_syrk_tiled(const l3_shape &s,real_t<PB> *al,real_t<PB> *be,real_t<PB> *a,real_t<PB> *b,real_t<PB> *c) {
    using O=syrk_fixed_ops<PB>;
    mplapack_blas::gpu_lane_args<O> g{};
    g.s=s; g.al=al; g.be=be; g.a=a; g.b=b; g.c=c;
    g.lda=s.lda; g.ldb=s.ldb; g.ldc=s.ldc;
#ifdef MPLAPACK_MPFR_CUDA_HOST_EMULATION
    return mplapack_blas::blas_device<O>(g);
#else
    constexpr int T=mplapack_blas::syrk_shared<O>::BLAS_TILE;
    const unsigned nb=(unsigned)((s.n+T-1)/T);
    if(nb) mplapack_blas::k_syrk_tiled<O><<<dim3(nb,nb),dim3(T,T)>>>(g);
    return cudaGetLastError();
#endif
}
}
#endif
