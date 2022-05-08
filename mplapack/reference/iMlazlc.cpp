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
iMlazlc(INTEGER const m, INTEGER const n, COMPLEX *a, INTEGER const lda) {
    INTEGER return_value = 0;
    //     Quick test for the common case where one corner is non-zero.
    const COMPLEX zero = COMPLEX(0.0, 0.0);
    INTEGER i = 0;
    if (n == 0) {
        return_value = n;
    } else if (a[(n - 1) * lda] != zero || a[(m - 1) + (n - 1) * lda] != zero) {
        return_value = n;
    } else {
        //     Now scan each column from the end, returning with the first non-zero.
        for (return_value = n; return_value >= 1; return_value = return_value - 1) {
            for (i = 1; i <= m; i = i + 1) {
                if (a[(i - 1) + (return_value - 1) * lda] != zero) {
                    return return_value;
                }
            }
        }
    }
    return return_value;
}
