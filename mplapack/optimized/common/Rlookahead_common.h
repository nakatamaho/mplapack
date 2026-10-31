/*
 * Copyright (c) 2026
 *      Nakata, Maho
 *      All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 *
 */

/*
 * Look-ahead for the blocked factorizations of the optimized libraries
 * (Rgetrf.cpp, Rpotrf.cpp, Cgetrf.cpp, Cpotrf.cpp).
 *
 * The steps of the reference blocked loop are split into OpenMP tasks on
 * tiles of nb columns (rows) with dependencies, so that the next panel is
 * factored as soon as its tile is updated, while the other tiles are still
 * being updated.  Every tile still receives the operations of the reference
 * routine in their order, on the same data, so the results are the
 * reference results bit for bit.
 *
 *   MPLAPACK_LAPACK_LOOKAHEAD=0   the reference blocked loop
 *   MPLAPACK_LAPACK_LOOKAHEAD=2   tasks also for small LU and for Cholesky
 *   MPLAPACK_LAPACK_FAST=1        dedicated LU panel width / Cholesky
 *                                 recursion (see README.md); rounding
 *                                 differs from the reference blocked loop
 */

#ifndef _MPLAPACK_LOOKAHEAD_COMMON_H_
#define _MPLAPACK_LOOKAHEAD_COMMON_H_

#include "../../../mpblas/optimized/common/Rblas_cpu_common.h"

// These CUDA flavors now bind every CPU factor/panel BLAS call to CPU
// implementations; their look-ahead policy can match the CPU library.
#if defined(MPLAPACK_ACCELERATOR_LIBRARY) && (defined(MPLAPACK_BUILD_WITH_DD) || defined(MPLAPACK_BUILD_WITH_TD) || defined(MPLAPACK_BUILD_WITH_QD) || defined(MPLAPACK_BUILD_WITH_MPFR))
#define MPLAPACK_FACTOR_CPU_BLAS
#endif

namespace mplapack_lapack {

// Recursive panels may share a team only when every worker has installed
// the caller's arithmetic context (not an arbitrary external OpenMP team).
inline thread_local bool lu_team_active = false;
struct lu_team_scope {
    bool previous = lu_team_active;
    lu_team_scope() { lu_team_active = true; }
    ~lu_team_scope() { lu_team_active = previous; }
};

inline long la_env(const char *name, long dflt) { return mplapack_blas::blas_env_long(name, dflt); }

inline bool la_fast() { return la_env("MPLAPACK_LAPACK_FAST", 0) != 0; }

// Whether the look-ahead (task) form is used: two threads or more, not
// already inside a parallel region, and
//   MPLAPACK_LAPACK_LOOKAHEAD=1 (default): LU with min(m,n) >= 12*nb;
//     also real LU >= 2*nb with more than eight requested threads. It then
//     shares one bounded team instead of repeatedly synchronizing full
//     BLAS teams. Cholesky never (the left-looking order of the reference
//     puts the diagonal blocks on a serial path: 0.66-1.17x)
//   MPLAPACK_LAPACK_LOOKAHEAD=2: always;  0: never
// Accelerator flavors still using per-call GPU BLAS default to no tasks:
// transfers/serialized GPU calls do not overlap. CUDA CPU-factor bindings
// instead use the CPU library policy (no transfers in these tasks).
inline bool la_lookahead(bool lu, long mn, long nb, bool small_real_lu = false)
{
#ifdef _OPENMP
#if defined(MPLAPACK_ACCELERATOR_LIBRARY) && !defined(MPLAPACK_FACTOR_CPU_BLAS)
    const long mode = la_env("MPLAPACK_LAPACK_LOOKAHEAD", 0);
#else
    const long mode = la_env("MPLAPACK_LAPACK_LOOKAHEAD", 1);
#endif
    if (mode == 0 || omp_get_max_threads() < 2 || omp_get_level() != 0)
        return false;
#if defined(MPLAPACK_ACCELERATOR_LIBRARY) && !defined(MPLAPACK_FACTOR_CPU_BLAS)
    return mode >= 2;
#else
    return mode >= 2 || (lu && (mn >= 12 * nb || (small_real_lu && omp_get_max_threads() > 8 && mn >= 2 * nb)));
#endif
#else
    (void)lu, (void)mn, (void)nb, (void)small_real_lu;
    return false;
#endif
}

// The blocked LU of the reference Rgetrf/Cgetrf (block size nb) as OpenMP
// tasks on the column tiles [t*nb, (t+1)*nb): for every step s (j = s*nb),
//   panel(j, jb, iinfo)          Rgetrf2 of the rows [j, m), tile s;
//   left(j, jb, c0, c1)          Rlaswp of step s on tile t < s;
//   step(j, jb, c0, c1)          Rlaswp, Rtrsm, Rgemm of step s on tile t > s;
// each task depends on the tiles it reads (in) and writes (inout), so every
// tile goes through the operations of the reference in their order, and
// the next panel can start as soon as its own tile is updated (look-ahead).
// The tasks run their BLAS calls on one thread.  ipiv holds 1-based
// indices, as in the reference; the panel task adds j to its pivots.
template <class FP, class FL, class FS> void la_getrf(long m, long n, long nb, INTEGER *ipiv, INTEGER &info, FP panel, FL left, FS step)
{
    const long mn = std::min(m, n), steps = (mn + nb - 1) / nb, tiles = (n + nb - 1) / nb;
    std::vector<INTEGER> iinfo((size_t)steps, 0);
    std::vector<char> dep((size_t)tiles);
    char *d = dep.data();
    const mplapack_blas::blas_thread_context ctx;
#ifdef _OPENMP
    // At most one update per column tile is ready at a time. Tall panels
    // can additionally expose independent 1024-row tasks. Do not wake a
    // full machine's team for a handful of small tiles (especially costly
    // with libgomp's active waits on a loaded host). No ICV is modified.
    const int threads = mn < 12 * nb
        ? int(std::min(long(omp_get_max_threads()), std::max(tiles, (m + 1023) / 1024)))
        : omp_get_max_threads();
#pragma omp parallel num_threads(threads)
#endif
    {
        const mplapack_blas::blas_thread_context::scope sc(ctx);
        const lu_team_scope team;
#pragma omp single
        for (long s = 0; s < steps; s++) {
            const long j = s * nb, jb = std::min(nb, mn - j);
#pragma omp task depend(inout : d[s])
            {
                panel(j, jb, iinfo[(size_t)s]);
                for (long i = j; i < std::min(m, j + jb); i++)
                    ipiv[i] += j;
            }
            for (long t = 0; t < s; t++) {
#pragma omp task depend(in : d[s]) depend(inout : d[t])
                left(j, jb, t * nb, (t + 1) * nb);
            }
            // the columns of tile s after a last panel narrower than nb
            if (j + jb < std::min(n, (s + 1) * nb)) {
#pragma omp task depend(inout : d[s])
                step(j, jb, j + jb, std::min(n, (s + 1) * nb));
            }
            for (long t = s + 1; t < tiles; t++) {
#pragma omp task depend(in : d[s]) depend(inout : d[t])
                step(j, jb, t * nb, std::min(n, (t + 1) * nb));
            }
        }
    }
    info = 0;
    for (long s = 0; s < steps; s++)
        if (info == 0 && iinfo[(size_t)s] > 0)
            info = iinfo[(size_t)s] + s * nb;
}

// The blocked Cholesky of the reference Rpotrf/Cpotrf as tasks on the tiles
// of nb columns (upper) or rows (lower): diag(j, jb, iinfo) (Rsyrk/Cherk and
// Rpotrf2 of tile s), step(j, jb, p0, p1) (Rgemm and Rtrsm of step s on tile
// t > s).  The reference stops at the first diagonal block that is not
// positive definite: the tasks of that step and the later ones do nothing.
template <class FD, class FS> void la_potrf(long n, long nb, INTEGER &info, FD diag, FS step)
{
    const long tiles = (n + nb - 1) / nb;
    std::vector<char> dep((size_t)tiles);
    char *d = dep.data();
    long failed = tiles; // the first step whose diagonal block failed
    INTEGER finfo = 0;
    const mplapack_blas::blas_thread_context ctx;
#pragma omp parallel
    {
        const mplapack_blas::blas_thread_context::scope sc(ctx);
#pragma omp single
        for (long s = 0; s < tiles; s++) {
            const long j = s * nb, jb = std::min(nb, n - j);
#pragma omp task depend(inout : d[s]) shared(failed, finfo)
            {
                long f;
#pragma omp atomic read
                f = failed;
                if (f > s) {
                    INTEGER iinfo = 0;
                    diag(j, jb, iinfo);
                    if (iinfo != 0) {
                        finfo = iinfo + j;
#pragma omp atomic write
                        failed = s;
                    }
                }
            }
            for (long t = s + 1; t < tiles; t++) {
#pragma omp task depend(in : d[s]) depend(inout : d[t]) shared(failed)
                {
                    long f;
#pragma omp atomic read
                    f = failed;
                    if (f > s)
                        step(j, jb, t * nb, std::min(n, (t + 1) * nb));
                }
            }
        }
    }
    info = failed < tiles ? finfo : 0;
}

} // namespace mplapack_lapack

#endif
