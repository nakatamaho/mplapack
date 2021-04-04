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

void Rgghrd(const char *compq, const char *compz, INTEGER const &n, INTEGER const &ilo, INTEGER const &ihi, REAL *a, INTEGER const &lda, REAL *b, INTEGER const &ldb, REAL *q, INTEGER const &ldq, REAL *z, INTEGER const &ldz, INTEGER &info) {
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
    //     .. Parameters ..
    //     ..
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
    //     Decode COMPQ
    //
    bool ilq = false;
    INTEGER icompq = 0;
    if (Mlsame(compq, "N")) {
        ilq = false;
        icompq = 1;
    } else if (Mlsame(compq, "V")) {
        ilq = true;
        icompq = 2;
    } else if (Mlsame(compq, "I")) {
        ilq = true;
        icompq = 3;
    } else {
        icompq = 0;
    }
    //
    //     Decode COMPZ
    //
    bool ilz = false;
    INTEGER icompz = 0;
    if (Mlsame(compz, "N")) {
        ilz = false;
        icompz = 1;
    } else if (Mlsame(compz, "V")) {
        ilz = true;
        icompz = 2;
    } else if (Mlsame(compz, "I")) {
        ilz = true;
        icompz = 3;
    } else {
        icompz = 0;
    }
    //
    //     Test the input parameters.
    //
    info = 0;
    if (icompq <= 0) {
        info = -1;
    } else if (icompz <= 0) {
        info = -2;
    } else if (n < 0) {
        info = -3;
    } else if (ilo < 1) {
        info = -4;
    } else if (ihi > n || ihi < ilo - 1) {
        info = -5;
    } else if (lda < max((INTEGER)1, n)) {
        info = -7;
    } else if (ldb < max((INTEGER)1, n)) {
        info = -9;
    } else if ((ilq && ldq < n) || ldq < 1) {
        info = -11;
    } else if ((ilz && ldz < n) || ldz < 1) {
        info = -13;
    }
    if (info != 0) {
        Mxerbla("Rgghrd", -info);
        return;
    }
    //
    //     Initialize Q and Z if desired.
    //
    const REAL zero = 0.0;
    const REAL one = 1.0;
    if (icompq == 3) {
        Rlaset("Full", n, n, zero, one, q, ldq);
    }
    if (icompz == 3) {
        Rlaset("Full", n, n, zero, one, z, ldz);
    }
    //
    //     Quick return if possible
    //
    if (n <= 1) {
        return;
    }
    //
    //     Zero out lower triangle of B
    //
    INTEGER jcol = 0;
    INTEGER jrow = 0;
    for (jcol = 1; jcol <= n - 1; jcol = jcol + 1) {
        for (jrow = jcol + 1; jrow <= n; jrow = jrow + 1) {
            b[(jrow - 1) + (jcol - 1) * ldb] = zero;
        }
    }
    //
    //     Reduce A and B
    //
    REAL temp = 0.0;
    REAL c = 0.0;
    REAL s = 0.0;
    for (jcol = ilo; jcol <= ihi - 2; jcol = jcol + 1) {
        //
        for (jrow = ihi; jrow >= jcol + 2; jrow = jrow - 1) {
            //
            //           Step 1: rotate rows JROW-1, JROW to kill A(JROW,JCOL)
            //
            temp = a[((jrow - 1) - 1) + (jcol - 1) * lda];
            Rlartg(temp, a[(jrow - 1) + (jcol - 1) * lda], c, s, a[((jrow - 1) - 1) + (jcol - 1) * lda]);
            a[(jrow - 1) + (jcol - 1) * lda] = zero;
            Rrot(n - jcol, a[((jrow - 1) - 1) + ((jcol + 1) - 1) * lda], lda, a[(jrow - 1) + ((jcol + 1) - 1) * lda], lda, c, s);
            Rrot(n + 2 - jrow, b[((jrow - 1) - 1) + ((jrow - 1) - 1) * ldb], ldb, b[(jrow - 1) + ((jrow - 1) - 1) * ldb], ldb, c, s);
            if (ilq) {
                Rrot(n, q[((jrow - 1) - 1) * ldq], 1, q[(jrow - 1) * ldq], 1, c, s);
            }
            //
            //           Step 2: rotate columns JROW, JROW-1 to kill B(JROW,JROW-1)
            //
            temp = b[(jrow - 1) + (jrow - 1) * ldb];
            Rlartg(temp, b[(jrow - 1) + ((jrow - 1) - 1) * ldb], c, s, b[(jrow - 1) + (jrow - 1) * ldb]);
            b[(jrow - 1) + ((jrow - 1) - 1) * ldb] = zero;
            Rrot(ihi, a[(jrow - 1) * lda], 1, a[((jrow - 1) - 1) * lda], 1, c, s);
            Rrot(jrow - 1, b[(jrow - 1) * ldb], 1, b[((jrow - 1) - 1) * ldb], 1, c, s);
            if (ilz) {
                Rrot(n, z[(jrow - 1) * ldz], 1, z[((jrow - 1) - 1) * ldz], 1, c, s);
            }
        }
    }
    //
    //     End of Rgghrd
    //
}
