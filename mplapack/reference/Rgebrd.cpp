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

void Rgebrd(INTEGER const m, INTEGER const n, REAL *a, INTEGER const lda, REAL *d, REAL *e, REAL *tauq, REAL *taup, REAL *work, INTEGER const lwork, int &info) {
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
    //     .. External Subroutines ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Test the input parameters
    //
    info = 0;
    INTEGER nb = max((INTEGER)1, iMlaenv(1, "Rgebrd", " ", m, n, -1, -1));
    INTEGER lwkopt = (m + n) * nb;
    work[1 - 1] = REAL(lwkopt);
    bool lquery = (lwork == -1);
    if (m < 0) {
        info = -1;
    } else if (n < 0) {
        info = -2;
    } else if (lda < max((INTEGER)1, m)) {
        info = -4;
    } else if (lwork < max({(INTEGER)1, m, n}) && !lquery) {
        info = -10;
    }
    if (info < 0) {
        Mxerbla("Rgebrd", -info);
        return;
    } else if (lquery) {
        return;
    }
    //
    //     Quick return if possible
    //
    INTEGER minmn = min(m, n);
    if (minmn == 0) {
        work[1 - 1] = 1;
        return;
    }
    //
    INTEGER ws = max(m, n);
    INTEGER ldwrkx = m;
    INTEGER ldwrky = n;
    //
    INTEGER nx = 0;
    INTEGER nbmin = 0;
    if (nb > 1 && nb < minmn) {
        //
        //        Set the crossover poINTEGER NX.
        //
        nx = max(nb, iMlaenv(3, "Rgebrd", " ", m, n, -1, -1));
        //
        //        Determine when to switch from blocked to unblocked code.
        //
        if (nx < minmn) {
            ws = (m + n) * nb;
            if (lwork < ws) {
                //
                //              Not enough work space for the optimal NB, consider using
                //              a smaller block size.
                //
                nbmin = iMlaenv(2, "Rgebrd", " ", m, n, -1, -1);
                if (lwork >= (m + n) * nbmin) {
                    nb = lwork / (m + n);
                } else {
                    nb = 1;
                    nx = minmn;
                }
            }
        }
    } else {
        nx = minmn;
    }
    //
    INTEGER i = 0;
    const REAL one = 1.0;
    INTEGER j = 0;
    for (i = 1; i <= minmn - nx; i = i + nb) {
        //
        //        Reduce rows and columns i:i+nb-1 to bidiagonal form and return
        //        the matrices X and Y which are needed to update the unreduced
        //        part of the matrix
        //
        Rlabrd(m - i + 1, n - i + 1, nb, &a[(i - 1) + (i - 1) * lda], lda, &d[i - 1], &e[i - 1], &tauq[i - 1], &taup[i - 1], work, ldwrkx, &work[(ldwrkx * nb + 1) - 1], ldwrky);
        //
        //        Update the trailing submatrix A(i+nb:m,i+nb:n), using an update
        //        of the form  A := A - V*Y**T - X*U**T
        //
        Rgemm("No transpose", "Transpose", m - i - nb + 1, n - i - nb + 1, nb, -one, &a[((i + nb) - 1) + (i - 1) * lda], lda, &work[(ldwrkx * nb + nb + 1) - 1], ldwrky, one, &a[((i + nb) - 1) + ((i + nb) - 1) * lda], lda);
        Rgemm("No transpose", "No transpose", m - i - nb + 1, n - i - nb + 1, nb, -one, &work[(nb + 1) - 1], ldwrkx, &a[(i - 1) + ((i + nb) - 1) * lda], lda, one, &a[((i + nb) - 1) + ((i + nb) - 1) * lda], lda);
        //
        //        Copy diagonal and off-diagonal elements of B back into A
        //
        if (m >= n) {
            for (j = i; j <= i + nb - 1; j = j + 1) {
                a[(j - 1) + (j - 1) * lda] = d[j - 1];
                a[(j - 1) + ((j + 1) - 1) * lda] = e[j - 1];
            }
        } else {
            for (j = i; j <= i + nb - 1; j = j + 1) {
                a[(j - 1) + (j - 1) * lda] = d[j - 1];
                a[((j + 1) - 1) + (j - 1) * lda] = e[j - 1];
            }
        }
    }
    //
    //     Use unblocked code to reduce the remainder of the matrix
    //
    INTEGER iinfo = 0;
    Rgebd2(m - i + 1, n - i + 1, &a[(i - 1) + (i - 1) * lda], lda, &d[i - 1], &e[i - 1], &tauq[i - 1], &taup[i - 1], work, iinfo);
    work[1 - 1] = ws;
    //
    //     End of Rgebrd
    //
}
