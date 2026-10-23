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
// result bit for bit; with MPLAPACK_LAPACK_FAST=1 the recursive Rpotrf2 on
// the whole matrix (another reference algorithm: faster, different rounding).

#define Rpotrf Rpotrf_reference_cpu
#include "../../reference/Rpotrf.cpp"
#undef Rpotrf

#include "Rlookahead_common.h"

void Rpotrf(const char *uplo, INTEGER const n, REAL *a, INTEGER const lda, INTEGER &info)
{
    const bool upper = Mlsame(uplo, "U");
    if ((!upper && !Mlsame(uplo, "L")) || n <= 0 || lda < std::max((INTEGER)1, n)) {
        Rpotrf_reference_cpu(uplo, n, a, lda, info);
        return;
    }
    if (mplapack_lapack::la_fast()) {
        info = 0;
        Rpotrf2(uplo, n, a, lda, info);
        return;
    }
    const INTEGER nb = iMlaenv(1, "Rpotrf", uplo, n, -1, -1, -1);
    if (nb <= 1 || nb >= n || !mplapack_lapack::la_lookahead(false, n, nb)) {
        Rpotrf_reference_cpu(uplo, n, a, lda, info);
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
            Rpotrf2(uplo, jb, &a[j + j * lda], lda, iinfo);
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
