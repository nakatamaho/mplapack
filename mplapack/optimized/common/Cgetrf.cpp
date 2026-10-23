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

// Cgetrf of the optimized libraries: as Rgetrf.cpp, with Cgetrf2, Claswp,
// Ctrsm and Cgemm.

#define Cgetrf Cgetrf_reference_cpu
#include "../../reference/Cgetrf.cpp"
#undef Cgetrf

#include "Rlookahead_common.h"

void Cgetrf(INTEGER const m, INTEGER const n, COMPLEX *a, INTEGER const lda, INTEGER *ipiv, INTEGER &info)
{
    if (m <= 0 || n <= 0 || lda < std::max((INTEGER)1, m)) {
        Cgetrf_reference_cpu(m, n, a, lda, ipiv, info);
        return;
    }
    if (mplapack_lapack::la_fast()) {
        info = 0;
        Cgetrf2(m, n, a, lda, ipiv, info);
        return;
    }
    const INTEGER nb = iMlaenv(1, "Cgetrf", " ", m, n, -1, -1);
    if (nb <= 1 || nb >= std::min(m, n) || !mplapack_lapack::la_lookahead(true, std::min(m, n), nb)) {
        Cgetrf_reference_cpu(m, n, a, lda, ipiv, info);
        return;
    }
    const COMPLEX one = COMPLEX(1.0, 0.0);
    mplapack_lapack::la_getrf(
        m, n, nb, ipiv, info, [&](long j, long jb, INTEGER &iinfo) { Cgetrf2(m - j, jb, &a[j + j * lda], lda, &ipiv[j], iinfo); },
        [&](long j, long jb, long c0, long c1) { Claswp(c1 - c0, &a[c0 * lda], lda, j + 1, j + jb, ipiv, 1); },
        [&](long j, long jb, long c0, long c1) {
            Claswp(c1 - c0, &a[c0 * lda], lda, j + 1, j + jb, ipiv, 1);
            Ctrsm("Left", "Lower", "No transpose", "Unit", jb, c1 - c0, one, &a[j + j * lda], lda, &a[j + c0 * lda], lda);
            if (j + jb < m)
                Cgemm("No transpose", "No transpose", m - j - jb, c1 - c0, jb, -one, &a[(j + jb) + j * lda], lda, &a[j + c0 * lda], lda, one, &a[(j + jb) + c0 * lda], lda);
        });
}
