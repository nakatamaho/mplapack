/*
 * One-pass chunked iamax on the CPU (factorization goal 13).
 *
 * The reference scan keeps the first index of the exact maximum (strict >
 * comparisons over per-element exact scores).  Maximum with first-index
 * tie-break is associative and commutative for exact scores, so contiguous
 * chunks scanned serially and combined in chunk order return the reference
 * index bit for bit, for every thread count.  The chunk length is fixed by n
 * alone (4096 elements); the chunks run under OpenMP when allowed and
 * sequentially otherwise, with the same result either way.
 */

#ifndef _MPLAPACK_RIAMAX_CHUNKS_COMMON_H_
#define _MPLAPACK_RIAMAX_CHUNKS_COMMON_H_

#include <algorithm>
#include <vector>
#include "Rblas_cpu_common.h"

namespace mplapack_blas {

constexpr long iamax_chunk_length = 4096;

template <class Score> bool iamax_chunks(long n, Score score_of, INTEGER &result)
{
    if (n < 2)
        return false;
    const long chunks = (n + iamax_chunk_length - 1) / iamax_chunk_length;
    std::vector<REAL> cmax((size_t)chunks);
    std::vector<INTEGER> cidx((size_t)chunks);
    volatile bool bad = false;
    const blas_thread_context ctx;
#ifdef _OPENMP
#pragma omp parallel if (chunks >= 2 && blas_env_long("MPLAPACK_BLAS_OMP", 1) != 0)
#endif
    {
        const blas_thread_context::scope sc(ctx);
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
        for (long q = 0; q < chunks; q++) {
            if (bad)
                continue;
            const long i0 = q * iamax_chunk_length;
            const long i1 = std::min(n, i0 + iamax_chunk_length);
            REAL dmax = score_of(i0);
            long best = i0;
            if (dmax != dmax) {
                bad = true;
                continue;
            }
            for (long i = i0 + 1; i < i1; i++) {
                const REAL v = score_of(i);
                if (v != v) {
                    bad = true;
                    break;
                }
                if (v > dmax) {
                    dmax = v;
                    best = i;
                }
            }
            cmax[(size_t)q] = dmax;
            cidx[(size_t)q] = best;
        }
    }
    if (bad)
        return false;
    // Chunk order is index order, so this combines exactly like the reference
    // scan over [0, n): strict > keeps the first maximum.
    REAL dmax = cmax[0];
    INTEGER best = cidx[0];
    for (long q = 1; q < chunks; q++) {
        if (cmax[(size_t)q] > dmax) {
            dmax = cmax[(size_t)q];
            best = cidx[(size_t)q];
        }
    }
    result = best + 1;
    return true;
}

} // namespace mplapack_blas

#endif
