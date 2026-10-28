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

// Rpotrf of the optimized libraries: the reference blocked algorithm (same
// block size, Rsyrk and Rpotrf2 on the diagonal blocks, Rgemm and Rtrsm
// on the others) with look-ahead (Rlookahead_common.h), with the reference
// result bit for bit. The diagonal recursion has task-tiled solves and
// reductions; FAST=1 uses a separately blocked whole-matrix recursion.

#define Rpotrf Rpotrf_reference_cpu
#include "../../reference/Rpotrf.cpp"
#undef Rpotrf

#include "Rpotrf_common.h"

namespace {
struct CholeskyOps {
    void leaf(bool upper, INTEGER n, REAL *a, INTEGER ld, INTEGER &info) const
    { Rpotrf2(upper ? "U" : "L", n, a, ld, info); }
    void solve(bool upper, INTEGER k, INTEGER n, REAL *a, REAL *b, INTEGER ld) const
    {
        Rtrsm(upper ? "L" : "R", upper ? "U" : "L", "T", "N",
              upper ? k : n, upper ? n : k, REAL(1), a, ld, b, ld);
    }
    void rank(bool upper, INTEGER n, INTEGER k, REAL *b, REAL *c, INTEGER ld) const
    { Rsyrk(upper ? "U" : "L", upper ? "T" : "N", n, k, REAL(-1), b, ld, REAL(1), c, ld); }
    void rank_columns(bool upper, INTEGER n, INTEGER k, REAL *b, REAL *c, INTEGER ld, INTEGER p, INTEGER end) const
    {
        using namespace mplapack_blas;
        blas_call<REAL> call;
        make_syrk<REAL>(call, [](const char *x, const char *y) { return Mlsame(x, y); },
                        false, upper ? "U" : "L", upper ? "T" : "N", n, k,
                        REAL(-1), b, ld, (REAL *)nullptr, 1, REAL(1), c, ld);
        op_engine<REAL> engine;
        blas_mats<op_engine<REAL>> mats;
        mats.al.p = &call.alpha; mats.al.ld = 1;
        mats.be.p = &call.beta; mats.be.ld = 1;
        mats.a = blas_cpu_mat(call.a); mats.c = blas_cpu_mat(call.c);
        for (long lane = long(p) * n; lane < long(end) * n; ++lane)
            blas_lane(engine, call.s, mats, lane);
    }
};
void potrf_diagonal(const char *uplo, INTEGER n, REAL *a, INTEGER ld, INTEGER &info, bool tasks = false)
{
    if (!tasks) {
        Rpotrf2(uplo, n, a, ld, info);
        return;
    }
    const CholeskyOps op;
    mplapack_potrf::recursive(Mlsame(uplo, "U"), n, a, ld, info,
                             mplapack_potrf::tile_size(), tasks, false, op);
}
// Reuse the generated blocked loop, without changing its oracle above.
#define Rpotrf potrf_blocked
#define Rpotrf2 potrf_diagonal
#include "../../reference/Rpotrf.cpp"
#undef Rpotrf2
#undef Rpotrf
} // namespace

#ifdef MPLAPACK_POTRF_CPU_DISPATCH
static void potrf_cpu(const char *uplo, INTEGER const n, REAL *a, INTEGER const lda, INTEGER &info)
#else
void Rpotrf(const char *uplo, INTEGER const n, REAL *a, INTEGER const lda, INTEGER &info)
#endif
{
    const bool upper = Mlsame(uplo, "U");
    if ((!upper && !Mlsame(uplo, "L")) || n <= 0 || lda < std::max((INTEGER)1, n)) {
        Rpotrf_reference_cpu(uplo, n, a, lda, info);
        return;
    }
    if (mplapack_lapack::la_fast()) {
        info = 0;
        const CholeskyOps op;
        mplapack_potrf::recursive(upper, n, a, lda, info,
                                 mplapack_potrf::leaf_size(), false, true, op);
        return;
    }
    const INTEGER nb = iMlaenv(1, "Rpotrf", uplo, n, -1, -1, -1);
    if (nb <= 1 || nb >= n || !mplapack_potrf::tasks_enabled()) {
        potrf_blocked(uplo, n, a, lda, info);
        return;
    }
    const REAL one = 1.0;
    mplapack_lapack::la_potrf(
        n, nb, info,
        [&](long j, long jb, INTEGER &iinfo) {
            if (upper)
                Rsyrk("Upper", "Transpose", jb, j, -one, &a[j * lda], lda, one, &a[j + j * lda], lda);
            else
                Rsyrk("Lower", "No transpose", jb, j, -one, &a[j], lda, one, &a[j + j * lda], lda);
            potrf_diagonal(uplo, jb, &a[j + j * lda], lda, iinfo, true);
        },
        [&](long j, long jb, long p0, long p1) {
            if (p1 <= p0)
                return;
            if (upper) {
                Rgemm("Transpose", "No transpose", jb, p1 - p0, j, -one, &a[j * lda], lda, &a[p0 * lda], lda, one, &a[j + p0 * lda], lda);
                Rtrsm("Left", "Upper", "Transpose", "Non-unit", jb, p1 - p0, one, &a[j + j * lda], lda, &a[j + p0 * lda], lda);
            } else {
                Rgemm("No transpose", "Transpose", p1 - p0, jb, j, -one, &a[p0], lda, &a[j], lda, one, &a[p0 + j * lda], lda);
                Rtrsm("Right", "Lower", "Transpose", "Non-unit", p1 - p0, jb, one, &a[j + j * lda], lda, &a[p0 + j * lda], lda);
            }
        });
}
