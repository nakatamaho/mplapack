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

void Rlasrt(const char *id, INTEGER const n, REAL *d, INTEGER &info) {
    INTEGER dir = 0;
    INTEGER stkpnt = 0;
    INTEGER stacklen = 32;
    INTEGER stack[2 * stacklen];
    INTEGER start = 0;
    INTEGER endd = 0;
    const INTEGER select = 20;
    INTEGER i = 0;
    INTEGER j = 0;
    REAL dmnmx = 0.0;
    REAL d1 = 0.0;
    REAL d2 = 0.0;
    REAL d3 = 0.0;
    REAL tmp = 0.0;
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
    //     .. External Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Test the input parameters.
    //
    info = 0;
    dir = -1;
    if (Mlsame(id, "D")) {
        dir = 0;
    } else if (Mlsame(id, "I")) {
        dir = 1;
    }
    if (dir == -1) {
        info = -1;
    } else if (n < 0) {
        info = -2;
    }
    if (info != 0) {
        Mxerbla("Rlasrt", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    if (n <= 1) {
        return;
    }
    //
    stkpnt = 1;
    stack[(1 - 1) + (1 - 1) * stacklen] = 1;
    stack[(2 - 1) + (1 - 1) * stacklen] = n;
statement_10:
    start = stack[(1 - 1) + (stkpnt - 1) * stacklen];
    endd = stack[(2 - 1) + (stkpnt - 1) * stacklen];
    stkpnt = stkpnt - 1;
    if (endd - start <= select && endd - start > 0) {
        //
        //        Do Insertion sort on D( START:ENDD )
        //
        if (dir == 0) {
            //
            //           Sort INTEGERo decreasing order
            //
            for (i = start + 1; i <= endd; i = i + 1) {
                for (j = i; j >= start + 1; j = j - 1) {
                    if (d[j - 1] > d[(j - 1) - 1]) {
                        dmnmx = d[j - 1];
                        d[j - 1] = d[(j - 1) - 1];
                        d[(j - 1) - 1] = dmnmx;
                    } else {
                        goto statement_30;
                    }
                }
            statement_30:;
            }
            //
        } else {
            //
            //           Sort INTEGERo increasing order
            //
            for (i = start + 1; i <= endd; i = i + 1) {
                for (j = i; j >= start + 1; j = j - 1) {
                    if (d[j - 1] < d[(j - 1) - 1]) {
                        dmnmx = d[j - 1];
                        d[j - 1] = d[(j - 1) - 1];
                        d[(j - 1) - 1] = dmnmx;
                    } else {
                        goto statement_50;
                    }
                }
            statement_50:;
            }
            //
        }
        //
    } else if (endd - start > select) {
        //
        //        Partition D( START:ENDD ) and stack parts, largest one first
        //
        //        Choose partition entry as median of 3
        //
        d1 = d[start - 1];
        d2 = d[endd - 1];
        i = (start + endd) / 2;
        d3 = d[i - 1];
        if (d1 < d2) {
            if (d3 < d1) {
                dmnmx = d1;
            } else if (d3 < d2) {
                dmnmx = d3;
            } else {
                dmnmx = d2;
            }
        } else {
            if (d3 < d2) {
                dmnmx = d2;
            } else if (d3 < d1) {
                dmnmx = d3;
            } else {
                dmnmx = d1;
            }
        }
        //
        if (dir == 0) {
            //
            //           Sort INTEGERo decreasing order
            //
            i = start - 1;
            j = endd + 1;
        statement_60:
        statement_70:
            j = j - 1;
            if (d[j - 1] < dmnmx) {
                goto statement_70;
            }
        statement_80:
            i++;
            if (d[i - 1] > dmnmx) {
                goto statement_80;
            }
            if (i < j) {
                tmp = d[i - 1];
                d[i - 1] = d[j - 1];
                d[j - 1] = tmp;
                goto statement_60;
            }
            if (j - start > endd - j - 1) {
                stkpnt++;
                stack[(1 - 1) + (stkpnt - 1) * stacklen] = start;
                stack[(2 - 1) + (stkpnt - 1) * stacklen] = j;
                stkpnt++;
                stack[(1 - 1) + (stkpnt - 1) * stacklen] = j + 1;
                stack[(2 - 1) + (stkpnt - 1) * stacklen] = endd;
            } else {
                stkpnt++;
                stack[(1 - 1) + (stkpnt - 1) * stacklen] = j + 1;
                stack[(2 - 1) + (stkpnt - 1) * stacklen] = endd;
                stkpnt++;
                stack[(1 - 1) + (stkpnt - 1) * stacklen] = start;
                stack[(2 - 1) + (stkpnt - 1) * stacklen] = j;
            }
        } else {
            //
            //           Sort into increasing order
            //
            i = start - 1;
            j = endd + 1;
        statement_90:
        statement_100:
            j = j - 1;
            if (d[j - 1] > dmnmx) {
                goto statement_100;
            }
        statement_110:
            i++;
            if (d[i - 1] < dmnmx) {
                goto statement_110;
            }
            if (i < j) {
                tmp = d[i - 1];
                d[i - 1] = d[j - 1];
                d[j - 1] = tmp;
                goto statement_90;
            }
            if (j - start > endd - j - 1) {
                stkpnt++;
                stack[(1 - 1) + (stkpnt - 1) * stacklen] = start;
                stack[(2 - 1) + (stkpnt - 1) * stacklen] = j;
                stkpnt++;
                stack[(1 - 1) + (stkpnt - 1) * stacklen] = j + 1;
                stack[(2 - 1) + (stkpnt - 1) * stacklen] = endd;
            } else {
                stkpnt++;
                stack[(1 - 1) + (stkpnt - 1) * stacklen] = j + 1;
                stack[(2 - 1) + (stkpnt - 1) * stacklen] = endd;
                stkpnt++;
                stack[(1 - 1) + (stkpnt - 1) * stacklen] = start;
                stack[(2 - 1) + (stkpnt - 1) * stacklen] = j;
            }
        }
    }
    if (stkpnt > 0) {
        goto statement_10;
    }
    //
    //     End of Rlasrt
    //
}
