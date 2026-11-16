#ifndef MPLAPACK_SOLVER_CPU_FALLBACK_H
#define MPLAPACK_SOLVER_CPU_FALLBACK_H
// A rejected solver must use CPU BLAS, not re-enter the CUDA dispatchers.
// In particular legacy MPFR fixed complex BLAS normalizes computed zeros;
// using it here loses the reference solver's signed-zero results.
namespace mplapack_solver_cpu_detail {
template<class> struct boolean_result;
template<class... A> struct boolean_result<void(A...)> { using type=bool(A...); };
}
#define MPLAPACK_SOLVER_CPU_BLAS(R) \
extern mplapack_solver_cpu_detail::boolean_result<decltype(R)>::type R##_omp; \
extern decltype(R) R##_reference_cpu; \
namespace { template<class... A> void solver_##R##_cpu(A... a) { \
    if(!R##_omp(a...)) R##_reference_cpu(a...); \
} }
MPLAPACK_SOLVER_CPU_BLAS(Rtrsm)
MPLAPACK_SOLVER_CPU_BLAS(Ctrsm)
MPLAPACK_SOLVER_CPU_BLAS(Cgemm)
MPLAPACK_SOLVER_CPU_BLAS(Rsyrk)
MPLAPACK_SOLVER_CPU_BLAS(Cherk)
MPLAPACK_SOLVER_CPU_BLAS(Rswap)
MPLAPACK_SOLVER_CPU_BLAS(Cswap)
MPLAPACK_SOLVER_CPU_BLAS(Rscal)
MPLAPACK_SOLVER_CPU_BLAS(Cscal)
#undef MPLAPACK_SOLVER_CPU_BLAS
#endif
