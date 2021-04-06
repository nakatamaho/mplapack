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

void Ctgexc(bool const wantq, bool const wantz, INTEGER const n, COMPLEX *a, INTEGER const lda, COMPLEX *b, INTEGER const ldb, COMPLEX *q, INTEGER const ldq, COMPLEX *z, INTEGER const ldz, INTEGER const ifst, INTEGER &ilst, INTEGER &info) {
    INTEGER here = 0;
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
    //     .. External Subroutines ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Decode and test input arguments.
    info = 0;
    if (n < 0) {
        info = -3;
    } else if (lda < max((INTEGER)1, n)) {
        info = -5;
    } else if (ldb < max((INTEGER)1, n)) {
        info = -7;
    } else if (ldq < 1 || wantq && (ldq < max((INTEGER)1, n))) {
        info = -9;
    } else if (ldz < 1 || wantz && (ldz < max((INTEGER)1, n))) {
        info = -11;
    } else if (ifst < 1 || ifst > n) {
        info = -12;
    } else if (ilst < 1 || ilst > n) {
        info = -13;
    }
    if (info != 0) {
        Mxerbla("Ctgexc", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    if (n <= 1) {
        return;
    }
    if (ifst == ilst) {
        return;
    }
    //
    if (ifst < ilst) {
        //
        here = ifst;
    //
    statement_10:
        //
        //        Swap with next one below
        //
        Ctgex2(wantq, wantz, n, a, lda, b, ldb, q, ldq, z, ldz, here, info);
        if (info != 0) {
            ilst = here;
            return;
        }
        here++;
        if (here < ilst) {
            goto statement_10;
        }
        here = here - 1;
    } else {
        here = ifst - 1;
    //
    statement_20:
        //
        //        Swap with next one above
        //
        Ctgex2(wantq, wantz, n, a, lda, b, ldb, q, ldq, z, ldz, here, info);
        if (info != 0) {
            ilst = here;
            return;
        }
        here = here - 1;
        if (here >= ilst) {
            goto statement_20;
        }
        here++;
    }
    ilst = here;
    //
    //     End of Ctgexc
    //
}
