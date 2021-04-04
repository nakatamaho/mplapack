/*
 * Copyright (c) 2021
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

#include <mpblas.h>
#include <mplapack.h>

void Rgttrs(const char *trans, INTEGER const &n, INTEGER const &nrhs, REAL *dl, REAL *d, REAL *du, REAL *du2, INTEGER *ipiv, REAL *b, INTEGER const &ldb, INTEGER &info) {
    //
    //  -- LAPACK computational routine --
    //  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    //  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    //
    //     .. Scalar Arguments ..
    //     ..
    //     .. Array Arguments ..
    //     ..
    //
    //  =====================================================================
    //
    //     .. Local Scalars ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    info = 0;
    bool notran = (trans == "N" || trans == "n");
    if (!notran && !(trans == "T" || trans == "t") && !(trans == "C" || trans == "c")) {
        info = -1;
    } else if (n < 0) {
        info = -2;
    } else if (nrhs < 0) {
        info = -3;
    } else if (ldb < max(n, 1)) {
        info = -10;
    }
    if (info != 0) {
        Mxerbla("Rgttrs", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    if (n == 0 || nrhs == 0) {
        return;
    }
    //
    //     Decode TRANS
    //
    INTEGER itrans = 0;
    if (notran) {
        itrans = 0;
    } else {
        itrans = 1;
    }
    //
    //     Determine the number of right-hand sides to solve at a time.
    //
    INTEGER nb = 0;
    if (nrhs == 1) {
        nb = 1;
    } else {
        nb = max((INTEGER)1, iMlaenv[("Rgttrs" - 1) * ldiMlaenv]);
    }
    //
    INTEGER j = 0;
    INTEGER jb = 0;
    if (nb >= nrhs) {
        Rgtts2(itrans, n, nrhs, dl, d, du, du2, ipiv, b, ldb);
    } else {
        for (j = 1; j <= nrhs; j = j + nb) {
            jb = min(nrhs - j + 1, nb);
            Rgtts2(itrans, n, jb, dl, d, du, du2, ipiv, b[(j - 1) * ldb], ldb);
        }
    }
    //
    //     End of Rgttrs
    //
}
