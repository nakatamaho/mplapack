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

INTEGER
iMlazlr(INTEGER const m, INTEGER const n, COMPLEX *a, INTEGER const lda) {
    INTEGER return_value = 0;
    //     Quick test for the common case where one corner is non-zero.
    const COMPLEX zero = COMPLEX(0.0, 0.0);
    INTEGER j = 0;
    INTEGER i = 0;
    if (m == 0) {
        return_value = m;
    } else if (a[(m - 1)] != zero || a[(m - 1) + (n - 1) * lda] != zero) {
        return_value = m;
    } else {
        //     Scan up each column tracking the last zero row seen.
        return_value = 0;
        for (j = 1; j <= n; j = j + 1) {
            i = m;
            while ((a[(max(i, (INTEGER)1) - 1) + (j - 1) * lda] == zero) && (i >= 1)) {
                i = i - 1;
            }
            return_value = max(return_value, i);
        }
    }
    return return_value;
}
