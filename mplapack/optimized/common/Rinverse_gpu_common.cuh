// Copyright (c) 2026 Nakata, Maho. BSD-2-Clause (see Rgetrf.cpp).
// BF adapter for the shared reference-order inverse scheduler.
#ifndef MPLAPACK_INVERSE_GPU_COMMON_CUH
#define MPLAPACK_INVERSE_GPU_COMMON_CUH
#include "Rpotrf_gpu_common.cuh"
#include "Rinverse_schedule.h"
namespace mplapack_inverse_gpu {
class bf_inverse_launcher : public mplapack_lapack_gpu::chol_schedule {
    mplapack_lapack_gpu::buffer<long> status;
  public:
    static constexpr bool tiled_mm=true;
    bf_inverse_launcher() { if(!error) error=status.allocate(1); }
    ~bf_inverse_launcher() { finish(); }
    long *status_data() { return status.data(); }
};
template<class Ops,bool Complex> using arithmetic=inverse_arithmetic<
    mplapack_blas_dev::dreal<Ops>,mplapack_blas_dev::dcomplex<mplapack_blas_dev::dreal<Ops>>,Complex>;
template<class Ops,bool Complex> using job=inverse_job<arithmetic<Ops,Complex>>;
template<class Ops,bool Complex> using context=inverse_schedule<arithmetic<Ops,Complex>,bf_inverse_launcher>;
}
#endif
