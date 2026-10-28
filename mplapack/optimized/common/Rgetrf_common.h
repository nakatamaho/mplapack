/* Dedicated LU scheduling; arithmetic and pivot semantics follow GETRF2.
 * Copyright (c) 2026 Nakata, Maho. BSD-2-Clause (see Rgetrf.cpp).
 */
#ifndef MPLAPACK_GETRF_COMMON_H
#define MPLAPACK_GETRF_COMMON_H
#include "Rlookahead_common.h"

namespace mplapack_getrf {
inline bool tasks_allowed()
{
#ifdef _OPENMP
#ifdef MPLAPACK_ACCELERATOR_LIBRARY
    return mplapack_lapack::la_env("MPLAPACK_LAPACK_LOOKAHEAD", 0) >= 2;
#else
    return mplapack_lapack::la_env("MPLAPACK_LAPACK_LOOKAHEAD", 1) != 0;
#endif
#else
    return false;
#endif
}
inline INTEGER width()
{
    return std::max(long(1), mplapack_lapack::la_env("MPLAPACK_GETRF_TILE", 32));
}

// A taskgroup, not a nested parallel region: usable inside a look-ahead
// panel task. Each worker inherits the caller's MPFR/GMP arithmetic context.
// Small panels stay serial to avoid a task per tiny column.
template <class F> void rows(INTEGER n, bool tasks, F f)
{
    if (!tasks || n < 4096) {
        f(0, n);
        return;
    }
#pragma omp taskgroup
    for (INTEGER p = 0; p < n; p += 1024) {
        const INTEGER end = std::min(n, p + 1024);
#pragma omp task firstprivate(p, end) shared(f)
        f(p, end);
    }
}

template <class T, class Ops>
void recursive(INTEGER m, INTEGER n, T *a, INTEGER ld, INTEGER *piv,
               INTEGER &info, bool tasks, const Ops &op)
{
    info = 0;
    if (m == 1) {
        piv[0] = 1;
        if (a[0] == T(0.0)) info = 1;
        return;
    }
    if (n == 1) {
        // Independent contiguous row blocks, reduced in increasing order:
        // strict > retains the first index on ties, just like iRamax/iCamax.
        // Keep the serial search for NaNs (their ordering is not total).
        const INTEGER blocks = tasks && m >= 4096 ? (m + 1023) / 1024 : 1;
        INTEGER pivot;
        if (blocks == 1) pivot = op.amax(m, a) - 1;
        else {
            std::vector<INTEGER> candidates((size_t)blocks);
            rows(m, tasks, [&](INTEGER lo, INTEGER hi) {
                candidates[lo / 1024] = lo + op.amax(hi - lo, a + lo) - 1;
            });
            pivot = candidates[0];
            bool nan = false;
            for (INTEGER b = 0; b < blocks; ++b) {
                const REAL v = op.mag(a[candidates[b]]);
                nan = nan || v != v;
                if (v > op.mag(a[pivot])) pivot = candidates[b];
            }
            if (nan) pivot = op.amax(m, a) - 1;
        }
        piv[0] = pivot + 1;
        if (a[pivot] == T(0.0)) { info = 1; return; }
        if (pivot != 0) { T tmp = a[0]; a[0] = a[pivot]; a[pivot] = tmp; }
        const bool scale = abs(a[0]) >= Rlamch("S");
        // Evaluate the reciprocal once, exactly as the reference SCAL call.
        const T factor = scale ? T(1.0) / a[0] : T(0.0);
        rows(m - 1, tasks, [&](INTEGER lo, INTEGER hi) {
            if (scale) op.scal(hi - lo, factor, a + 1 + lo);
            else for (INTEGER i = lo + 1; i <= hi; ++i) a[i] = a[i] / a[0];
        });
        return;
    }
    const INTEGER n1 = std::min(m, n) / 2, n2 = n - n1;
    INTEGER iinfo = 0;
    recursive(m, n1, a, ld, piv, info, tasks, op);
    auto update = [&](INTEGER lo, INTEGER hi) {
        T *b = a + (n1 + lo) * ld;
        op.swap(hi - lo, b, ld, 1, n1, piv);
        op.solve(n1, hi - lo, a, ld, b);
        op.update(m - n1, hi - lo, n1, a + n1, ld, b, b + n1);
    };
    if (tasks && n2 > width()) {
#pragma omp taskgroup
        for (INTEGER c = 0; c < n2; c += width()) {
            const INTEGER end = std::min(n2, c + width());
#pragma omp task firstprivate(c, end) shared(update)
            update(c, end);
        }
    } else update(0, n2);
    recursive(m - n1, n2, a + n1 + n1 * ld, ld, piv + n1, iinfo, tasks, op);
    if (info == 0 && iinfo > 0) info = iinfo + n1;
    for (INTEGER i = n1; i < std::min(m, n); ++i) piv[i] += n1;
    op.swap(n1, a, ld, n1 + 1, std::min(m, n), piv);
}
} // namespace mplapack_getrf
#endif
