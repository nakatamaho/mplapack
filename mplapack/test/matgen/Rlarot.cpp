/*
 * Copyright (c) 2008-2021
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

#include <mplapack_matgen.h>

void Rlarot(bool const lrows, bool const lleft, bool const lright, INTEGER const nl, REAL const c, REAL const s, REAL *a, INTEGER const lda, REAL &xleft, REAL &xright) {
    //
    //  -- LAPACK auxiliary routine --
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
    //     .. Local Arrays ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Set up indices, arrays for ends
    //
    INTEGER iinc = 0;
    INTEGER inext = 0;
    if (lrows) {
        iinc = lda;
        inext = 1;
    } else {
        iinc = 1;
        inext = lda;
    }
    //
    INTEGER nt = 0;
    INTEGER ix = 0;
    INTEGER iy = 0;
    REAL xt[2];
    REAL yt[2];
    if (lleft) {
        nt = 1;
        ix = 1 + iinc;
        iy = 2 + lda;
        xt[1 - 1] = a[1 - 1];
        yt[1 - 1] = xleft;
    } else {
        nt = 0;
        ix = 1;
        iy = 1 + inext;
    }
    //
    INTEGER iyt = 0;
    if (lright) {
        iyt = 1 + inext + (nl - 1) * iinc;
        nt++;
        xt[nt - 1] = xright;
        yt[nt - 1] = a[iyt - 1];
    }
    //
    //     Check for errors
    //
    if (nl < nt) {
        Mxerbla("Rlarot", 4);
        return;
    }
    if (lda <= 0 || (!lrows && lda < nl - nt)) {
        Mxerbla("Rlarot", 8);
        return;
    }
    //
    //     Rotate
    //
    Rrot(nl - nt, &a[ix - 1], iinc, &a[iy - 1], iinc, c, s);
    Rrot(nt, xt, 1, yt, 1, c, s);
    //
    //     Stuff values back into XLEFT, XRIGHT, etc.
    //
    if (lleft) {
        a[1 - 1] = xt[1 - 1];
        xleft = yt[1 - 1];
    }
    //
    if (lright) {
        xright = xt[nt - 1];
        a[iyt - 1] = yt[nt - 1];
    }
    //
    //     End of Rlarot
    //
}
