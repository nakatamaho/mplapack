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

void Rgbsv(INTEGER const n, INTEGER const kl, INTEGER const ku, INTEGER const nrhs, REAL *ab, INTEGER const ldab, INTEGER *ipiv, REAL *b, INTEGER const ldb, INTEGER &info) {
    //
    //     Test the input parameters.
    //
    info = 0;
    if (n < 0) {
        info = -1;
    } else if (kl < 0) {
        info = -2;
    } else if (ku < 0) {
        info = -3;
    } else if (nrhs < 0) {
        info = -4;
    } else if (ldab < 2 * kl + ku + 1) {
        info = -6;
    } else if (ldb < max(n, (INTEGER)1)) {
        info = -9;
    }
    if (info != 0) {
        Mxerbla("Rgbsv", -info);
        return;
    }
    //
    //     Compute the LU factorization of the band matrix A.
    //
    Rgbtrf(n, n, kl, ku, ab, ldab, ipiv, info);
    if (info == 0) {
        //
        //        Solve the system A*X = B, overwriting B with X.
        //
        Rgbtrs("No transpose", n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info);
    }
    //
    //     End of Rgbsv
    //
}
