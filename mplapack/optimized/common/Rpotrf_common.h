/*
 * Copyright (c) 2026 Nakata, Maho. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 */

// CPU Cholesky scheduling helpers. The arithmetic belongs to the BLAS and
// the unchanged reference potrf2 leaf; only independent outputs are tiled.
#ifndef MPLAPACK_POTRF_COMMON_H
#define MPLAPACK_POTRF_COMMON_H

#include "Rlookahead_common.h"

namespace mplapack_potrf {

// CPU and accelerator builds share the single accelerator-library switch.
// Concurrent GPU BLAS is safe after b626121f, but only explicitly requested.
inline bool tasks_enabled()
{
#if defined(_OPENMP)
    return mplapack_lapack::la_lookahead(false, 0, 1);
#else
    return false;
#endif
}

inline INTEGER tile_size()
{
#if defined(MPLAPACK_BUILD_WITH_MPFR) || defined(MPLAPACK_BUILD_WITH_GMP)
    const long dflt = 8;
#else
    const long dflt = 16;
#endif
    return std::max(long(1), mplapack_lapack::la_env("MPLAPACK_POTRF_TILE", dflt));
}

inline INTEGER leaf_size()
{
#if defined(MPLAPACK_BUILD_WITH_MPFR) || defined(MPLAPACK_BUILD_WITH_GMP)
    const long dflt = 32;
#else
    const long dflt = 64;
#endif
    return std::max(long(1), mplapack_lapack::la_env("MPLAPACK_POTRF_LEAF", dflt));
}

// Same binary recursion as potrf2 in default mode. In fast mode the split
// is aligned to our own leaf size, rather than to iMlaenv's blocked loop.
// taskgroup is local to a recursion node, not a team barrier: panel tiles
// can begin immediately after their diagonal block finishes. All workers
// retain the caller's precision/rounding context (installed by the team).
template <class T, class Ops>
void recursive(bool upper, INTEGER n, T *a, INTEGER ld, INTEGER &info,
               INTEGER tile, bool tasks, bool fast, const Ops &op)
{
    info = 0;
    if (n <= tile) {
        op.leaf(upper, n, a, ld, info);
        return;
    }
    INTEGER n1 = n / 2;
    if (fast && n1 >= tile)
        n1 = (n1 / tile) * tile;
    const INTEGER n2 = n - n1;
    recursive(upper, n1, a, ld, info, tile, tasks, fast, op);
    if (info != 0)
        return;
    T *panel = a + (upper ? n1 * ld : n1);
    T *diag = a + n1 + n1 * ld;
    if (tasks && n2 > tile) {
#pragma omp taskgroup
        for (INTEGER p = 0; p < n2; p += tile) {
            const INTEGER count = std::min(tile, n2 - p);
#pragma omp task firstprivate(p, count, n1, panel, a, upper, ld) shared(op)
            op.solve(upper, n1, count, a, panel + (upper ? p * ld : p), ld);
        }
        // No arithmetic splitting in k: each output gets exactly the full
        // reference syrk/herk reduction, including its alpha/beta handling.
#pragma omp taskgroup
        for (INTEGER p = 0; p < n2; p += tile) {
            const INTEGER end = std::min(n2, p + tile);
#pragma omp task firstprivate(p, end, n1, n2, panel, diag, upper, ld) shared(op)
            op.rank_columns(upper, n2, n1, panel, diag, ld, p, end);
        }
    } else {
        op.solve(upper, n1, n2, a, panel, ld);
        op.rank(upper, n2, n1, panel, diag, ld);
    }
    recursive(upper, n2, diag, ld, info, tile, tasks, fast, op);
    if (info != 0)
        info += n1;
}

} // namespace mplapack_potrf
#endif
