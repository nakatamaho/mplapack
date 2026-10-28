/* Copyright (c) 2026 Nakata, Maho. BSD-2-Clause (see Rgetrf.cpp). */
#define Rgetrf2 Rgetrf2_reference_cpu
#include "../../reference/Rgetrf2.cpp"
#undef Rgetrf2
#include "Rgetrf_common.h"

namespace {
struct real_lu_ops {
    INTEGER amax(INTEGER n, REAL *a) const { return iRamax(n, a, 1); }
    REAL mag(const REAL &a) const { return abs(a); }
    void scal(INTEGER n, const REAL &s, REAL *a) const { Rscal(n, s, a, 1); }
    void swap(INTEGER n, REAL *a, INTEGER ld, INTEGER first, INTEGER last, INTEGER *p) const
    { Rlaswp(n, a, ld, first, last, p, 1); }
    void solve(INTEGER m, INTEGER n, REAL *a, INTEGER ld, REAL *b) const
    { Rtrsm("L", "L", "N", "U", m, n, REAL(1.0), a, ld, b, ld); }
    void update(INTEGER m, INTEGER n, INTEGER k, REAL *a, INTEGER ld, REAL *b, REAL *c) const
    { Rgemm("N", "N", m, n, k, -REAL(1.0), a, ld, b, ld, REAL(1.0), c, ld); }
};
}
void Rgetrf2(INTEGER const m, INTEGER const n, REAL *a, INTEGER const lda, INTEGER *ipiv, INTEGER &info)
{
    if (m <= 0 || n <= 0 || lda < std::max(INTEGER(1), m)) {
        Rgetrf2_reference_cpu(m, n, a, lda, ipiv, info);
        return;
    }
    // Serial leaves need no scheduler. Reuse the unchanged recursion rather
    // than paying task-policy/tiling overhead at every small recursion node.
    if (m < 4096 && n <= 64) {
        Rgetrf2_reference_cpu(m, n, a, lda, ipiv, info);
        return;
    }
    const bool tasks = mplapack_getrf::tasks_allowed();
#ifdef _OPENMP
    if (tasks && omp_get_level() == 0 && omp_get_max_threads() > 1 && (m >= 4096 || std::min(m, n) >= 128)) {
        const mplapack_blas::blas_thread_context ctx;
#pragma omp parallel
        {
            const mplapack_blas::blas_thread_context::scope sc(ctx);
#pragma omp single
            mplapack_getrf::recursive(m, n, a, lda, ipiv, info, true, real_lu_ops{});
        }
        return;
    }
#endif
    mplapack_getrf::recursive(m, n, a, lda, ipiv, info,
#ifdef _OPENMP
                             tasks && mplapack_lapack::lu_team_active,
#else
                             false,
#endif
                             real_lu_ops{});
}
