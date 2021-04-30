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

void Rorgtsqr_row(INTEGER const m, INTEGER const n, INTEGER const mb, INTEGER const nb, REAL *a, INTEGER const lda, REAL *t, INTEGER const ldt, REAL *work, INTEGER const lwork, INTEGER &info) {
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
    //     .. Local Arrays ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Test the input parameters
    //
    info = 0;
    bool lquery = lwork == -1;
    if (m < 0) {
        info = -1;
    } else if (n < 0 || m < n) {
        info = -2;
    } else if (mb <= n) {
        info = -3;
    } else if (nb < 1) {
        info = -4;
    } else if (lda < max((INTEGER)1, m)) {
        info = -6;
    } else if (ldt < max({(INTEGER)1, min(nb, n)})) {
        info = -8;
    } else if (lwork < 1 && !lquery) {
        info = -10;
    }
    //
    INTEGER nblocal = min(nb, n);
    //
    //     Determine the workspace size.
    //
    INTEGER lworkopt = 0;
    if (info == 0) {
        lworkopt = nblocal * max(nblocal, (n - nblocal));
    }
    //
    //     Handle error in the input parameters and handle the workspace query.
    //
    if (info != 0) {
        Mxerbla("Rorgtsqr_row", -info);
        return;
    } else if (lquery) {
        work[1 - 1] = castREAL(lworkopt);
        return;
    }
    //
    //     Quick return if possible
    //
    if (min(m, n) == 0) {
        work[1 - 1] = castREAL(lworkopt);
        return;
    }
    //
    //     (0) Set the upper-triangular part of the matrix A to zero and
    //     its diagonal elements to one.
    //
    const REAL zero = 0.0;
    const REAL one = 1.0;
    Rlaset("U", m, n, zero, one, a, lda);
    //
    //     KB_LAST is the column index of the last column block reflector
    //     in the matrices T and V.
    //
    INTEGER kb_last = ((n - 1) / nblocal) * nblocal + 1;
    //
    //     (1) Bottom-up loop over row blocks of A, except the top row block.
    //     NOTE: If MB>=M, then the loop is never executed.
    //
    INTEGER mb2 = 0;
    INTEGER m_plus_one = 0;
    INTEGER itmp = 0;
    INTEGER ib_bottom = 0;
    INTEGER num_all_row_blocks = 0;
    INTEGER jb_t = 0;
    INTEGER ib = 0;
    INTEGER imb = 0;
    INTEGER kb = 0;
    INTEGER knb = 0;
    if (mb < m) {
        //
        //        MB2 is the row blocking size for the row blocks before the
        //        first top row block in the matrix A. IB is the row index for
        //        the row blocks in the matrix A before the first top row block.
        //        IB_BOTTOM is the row index for the last bottom row block
        //        in the matrix A. JB_T is the column index of the corresponding
        //        column block in the matrix T.
        //
        //        Initialize variables.
        //
        //        NUM_ALL_ROW_BLOCKS is the number of row blocks in the matrix A
        //        including the first row block.
        //
        mb2 = mb - n;
        m_plus_one = m + 1;
        itmp = (m - mb - 1) / mb2;
        ib_bottom = itmp * mb2 + mb + 1;
        num_all_row_blocks = itmp + 2;
        jb_t = num_all_row_blocks * n + 1;
        //
        for (ib = ib_bottom; ib >= mb + 1; ib = ib - mb2) {
            //
            //           Determine the block size IMB for the current row block
            //           in the matrix A.
            //
            imb = min(m_plus_one - ib, mb2);
            //
            //           Determine the column index JB_T for the current column block
            //           in the matrix T.
            //
            jb_t = jb_t - n;
            //
            //           Apply column blocks of H in the row block from right to left.
            //
            //           KB is the column index of the current column block reflector
            //           in the matrices T and V.
            //
            for (kb = kb_last; kb >= 1; kb = kb - nblocal) {
                //
                //              Determine the size of the current column block KNB in
                //              the matrices T and V.
                //
                knb = min(nblocal, n - kb + 1);
                //
                Rlarfb_gett("I", imb, n - kb + 1, knb, &t[((jb_t + kb - 1) - 1) * ldt], ldt, &a[(kb - 1) + (kb - 1) * lda], lda, &a[(ib - 1) + (kb - 1) * lda], lda, work, knb);
                //
            }
            //
        }
        //
    }
    //
    //     (2) Top row block of A.
    //     NOTE: If MB>=M, then we have only one row block of A of size M
    //     and we work on the entire matrix A.
    //
    INTEGER mb1 = min(mb, m);
    //
    //     Apply column blocks of H in the top row block from right to left.
    //
    //     KB is the column index of the current block reflector in
    //     the matrices T and V.
    //
    REAL dummy[1];
    for (kb = kb_last; kb >= 1; kb = kb - nblocal) {
        //
        //        Determine the size of the current column block KNB in
        //        the matrices T and V.
        //
        knb = min(nblocal, n - kb + 1);
        //
        if (mb1 - kb - knb + 1 == 0) {
            //
            //           In SLARFB_GETT parameters, when M=0, then the matrix B
            //           does not exist, hence we need to pass a dummy array
            //           reference DUMMY(1,1) to B with LDDUMMY=1.
            //
            Rlarfb_gett("N", 0, n - kb + 1, knb, &t[(kb - 1) * ldt], ldt, &a[(kb - 1) + (kb - 1) * lda], lda, dummy, 1, work, knb);
        } else {
            Rlarfb_gett("N", mb1 - kb - knb + 1, n - kb + 1, knb, &t[(kb - 1) * ldt], ldt, &a[(kb - 1) + (kb - 1) * lda], lda, &a[((kb + knb) - 1) + (kb - 1) * lda], lda, work, knb);
            //
        }
        //
    }
    //
    work[1 - 1] = castREAL(lworkopt);
    //
    //     End of Rorgtsqr_row
    //
}
