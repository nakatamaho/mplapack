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

REAL Rlaran(INTEGER *iseed) {
    REAL return_value = 0.0;
    const INTEGER m4 = 2549;
    INTEGER it4 = 0;
    const INTEGER ipw2 = 4096;
    INTEGER it3 = 0;
    const INTEGER m3 = 2508;
    INTEGER it2 = 0;
    const INTEGER m2 = 322;
    INTEGER it1 = 0;
    const INTEGER m1 = 494;
    const REAL one = 1.0;
    const REAL r = one / ipw2;
    REAL rndout = 0.0;
//
//  -- LAPACK auxiliary routine --
//  -- LAPACK is a software package provided by Univ. of Tennessee,    --
//  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
//
//     .. Array Arguments ..
//     ..
//
//  =====================================================================
//
//     .. Parameters ..
//     ..
//     .. Local Scalars ..
//     ..
//     .. Intrinsic Functions ..
//     ..
//     .. Executable Statements ..
statement_10:
    //
    //     multiply the seed by the multiplier modulo 2**48
    //
    it4 = iseed[4 - 1] * m4;
    it3 = it4 / ipw2;
    it4 = it4 - ipw2 * it3;
    it3 += iseed[3 - 1] * m4 + iseed[4 - 1] * m3;
    it2 = it3 / ipw2;
    it3 = it3 - ipw2 * it2;
    it2 += iseed[2 - 1] * m4 + iseed[3 - 1] * m3 + iseed[4 - 1] * m2;
    it1 = it2 / ipw2;
    it2 = it2 - ipw2 * it1;
    it1 += iseed[1 - 1] * m4 + iseed[2 - 1] * m3 + iseed[3 - 1] * m2 + iseed[4 - 1] * m1;
    it1 = mod(it1, ipw2);
    //
    //     return updated seed
    //
    iseed[1 - 1] = it1;
    iseed[2 - 1] = it2;
    iseed[3 - 1] = it3;
    iseed[4 - 1] = it4;
    //
    //     convert 48-bit integer to a real number in the interval (0,1)
    //
    rndout = r * (it1.real() + r * (it2.real() + r * (it3.real() + r * (it4.real()))));
    //
    if (rndout == 1.0) {
        //        If a real number has n bits of precision, and the first
        //        n bits of the 48-bit integer above happen to be all 1 (which
        //        will occur about once every 2**n calls), then Rlaran will
        //        be rounded to exactly 1.0.
        //        Since Rlaran is not supposed to return exactly 0.0 or 1.0
        //        (and some callers of Rlaran, such as CLARND, depend on that),
        //        the statistically correct thing to do in this situation is
        //        simply to iterate again.
        //        N.B. the case Rlaran = 0.0 should not be possible.
        //
        goto statement_10;
    }
    //
    return_value = rndout;
    return return_value;
    //
    //     End of Rlaran
    //
}
