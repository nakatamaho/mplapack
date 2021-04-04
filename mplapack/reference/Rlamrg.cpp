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

void Rlamrg(INTEGER const &n1, INTEGER const &n2, REAL *a, INTEGER const &dtrd1, INTEGER const &dtrd2, arr_ref<INTEGER> index) {
    INTEGER n1sv = 0;
    INTEGER n2sv = 0;
    INTEGER ind1 = 0;
    INTEGER ind2 = 0;
    INTEGER i = 0;
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
    //     .. Executable Statements ..
    //
    n1sv = n1;
    n2sv = n2;
    if (dtrd1 > 0) {
        ind1 = 1;
    } else {
        ind1 = n1;
    }
    if (dtrd2 > 0) {
        ind2 = 1 + n1;
    } else {
        ind2 = n1 + n2;
    }
    i = 1;
//     while ( (N1SV > 0) & (N2SV > 0) )
statement_10:
    if (n1sv > 0 && n2sv > 0) {
        if (a[ind1 - 1] <= a[ind2 - 1]) {
            index[i - 1] = ind1;
            i++;
            ind1 += dtrd1;
            n1sv = n1sv - 1;
        } else {
            index[i - 1] = ind2;
            i++;
            ind2 += dtrd2;
            n2sv = n2sv - 1;
        }
        goto statement_10;
    }
    //     end while
    if (n1sv == 0) {
        for (n1sv = 1; n1sv <= n2sv; n1sv = n1sv + 1) {
            index[i - 1] = ind2;
            i++;
            ind2 += dtrd2;
        }
    } else {
        //     N2SV .EQ. 0
        for (n2sv = 1; n2sv <= n1sv; n2sv = n2sv + 1) {
            index[i - 1] = ind1;
            i++;
            ind1 += dtrd1;
        }
    }
    //
    //     End of Rlamrg
    //
}
