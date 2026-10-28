/* Copyright (c) 2026 Nakata, Maho. BSD-2-Clause (see Rgetrf.cpp). */
#define Cgetrf2 Cgetrf2_reference_cpu
#include "../../reference/Cgetrf2.cpp"
#undef Cgetrf2
#include "Rgetrf_common.h"

namespace {
struct complex_lu_ops {
    INTEGER amax(INTEGER n, COMPLEX *a) const { return iCamax(n, a, 1); }
    REAL mag(const COMPLEX &a) const { return abs(a.real()) + abs(a.imag()); }
    void scal(INTEGER n, const COMPLEX &s, COMPLEX *a) const { Cscal(n, s, a, 1); }
    void swap(INTEGER n, COMPLEX *a, INTEGER ld, INTEGER first, INTEGER last, INTEGER *p) const
    { Claswp(n, a, ld, first, last, p, 1); }
    void solve(INTEGER m, INTEGER n, COMPLEX *a, INTEGER ld, COMPLEX *b) const
    { Ctrsm("L", "L", "N", "U", m, n, COMPLEX(1.0), a, ld, b, ld); }
    void update(INTEGER m, INTEGER n, INTEGER k, COMPLEX *a, INTEGER ld, COMPLEX *b, COMPLEX *c) const
    { Cgemm("N", "N", m, n, k, -COMPLEX(1.0), a, ld, b, ld, COMPLEX(1.0), c, ld); }
};
}
void Cgetrf2(INTEGER const m, INTEGER const n, COMPLEX *a, INTEGER const lda, INTEGER *ipiv, INTEGER &info)
{
    if (m <= 0 || n <= 0 || lda < std::max(INTEGER(1), m)) {
        Cgetrf2_reference_cpu(m, n, a, lda, ipiv, info);
        return;
    }
    // The unchanged serial recursion is the leaf for ordinary small panels.
    if (m < 4096 && n <= 64) {
        Cgetrf2_reference_cpu(m, n, a, lda, ipiv, info);
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
            mplapack_getrf::recursive(m, n, a, lda, ipiv, info, true, complex_lu_ops{});
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
                             complex_lu_ops{});
}
